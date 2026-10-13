"""Incremental JSON grammar using the standard decoder for scalar values."""

import json
import re


_CHUNK_SIZE = 64 * 1024
_DECODER = json.JSONDecoder()
_NUMBER = re.compile(r"[0-9eE+.-]+")
_LITERALS = {"n": "null", "t": "true", "f": "false", "N": "NaN", "I": "Infinity"}


class _Tokens:
    def __init__(self, handle, decoder=None):
        self.handle, self.buffer, self.position, self.offset = handle, "", 0, 0
        self.line, self.column, self.eof = 1, 1, False
        self.decoder = decoder

    def refill(self):
        consumed = self.buffer[:self.position]
        lines = consumed.count("\n")
        self.line += lines
        self.column = len(consumed.rsplit("\n", 1)[-1]) + 1 if lines else self.column + len(consumed)
        self.offset += self.position
        self.buffer, self.position = self.buffer[self.position:], 0
        chunk = self.handle.read(_CHUNK_SIZE)
        if not isinstance(chunk, str):
            raise TypeError("Incremental JSON input must be a text stream")
        self.buffer += chunk
        self.eof = not chunk

    def error(self, message, position=None):
        error = json.JSONDecodeError(message, self.buffer, self.position if position is None else position)
        if error.lineno == 1:
            error.colno += self.column - 1
        error.lineno += self.line - 1
        error.pos += self.offset
        error.args = ("%s: line %s column %s (char %s)" %
                      (message, error.lineno, error.colno, error.pos),)
        return error

    def scalar(self):
        first = self.buffer[self.position]
        if first == '"':
            while True:
                try:
                    value, end = _DECODER.raw_decode(self.buffer, self.position)
                except json.JSONDecodeError as error:
                    incomplete = (error.msg.startswith("Unterminated string") or
                                  error.msg == "Invalid \\uXXXX escape" and error.pos + 5 >= len(self.buffer))
                    if incomplete and not self.eof:
                        self.refill()
                        continue
                    raise self.error(error.msg, error.pos) from None
                self.position = end
                return value
        if first == "-" and self.position + 1 == len(self.buffer) and not self.eof:
            self.refill()
        literal = "-Infinity" if self.buffer.startswith("-I", self.position) else _LITERALS.get(first)
        if literal:
            while len(self.buffer) - self.position < len(literal) and not self.eof:
                self.refill()
            if not self.buffer.startswith(literal, self.position):
                raise self.error("Expecting value")
            self.position += len(literal)
            return _DECODER.raw_decode(literal)[0]
        if first not in "-0123456789":
            raise self.error("Expecting value")
        while True:
            match = _NUMBER.match(self.buffer, self.position)
            if match.end() < len(self.buffer) or self.eof:
                break
            self.refill()
        token = match.group()
        try:
            value, end = _DECODER.raw_decode(token)
        except json.JSONDecodeError as error:
            raise self.error(error.msg, self.position + error.pos) from None
        if end != len(token):
            raise self.error("Invalid number", self.position + end)
        self.position = match.end()
        return value

    def __iter__(self):
        while True:
            while self.position < len(self.buffer) and self.buffer[self.position] in " \t\r\n":
                self.position += 1
            if self.position == len(self.buffer):
                if self.eof:
                    return
                self.refill()
                continue
            char = self.buffer[self.position]
            if char in "[{" and self.decoder is not None:
                # Decode complete subtrees already inside this bounded buffer
                # in C. Never read more merely to finish a container here.
                try:
                    value, end = self.decoder.raw_decode(self.buffer, self.position)
                except json.JSONDecodeError:
                    pass  # The iterative grammar handles partial containers.
                else:
                    self.position = end
                    yield "scalar", value
                    continue
            if char in "{}[]:,":
                self.position += 1
                yield char, None
            else:
                yield "scalar", self.scalar()


def json_events(handle, *, decoder=None):
    """Yield validated scalar/container events without a whole-document buffer."""
    tokens = _Tokens(handle, decoder)
    stack = [["root", "value"]]
    for token, value in tokens:
        kind, state = stack[-1]
        closing = "}" if kind == "map" else "]"
        if state in ("first_key", "key"):
            if state == "first_key" and token == "}":
                stack.pop()
                yield "end_map", None
            elif token == "scalar" and isinstance(value, str):
                stack[-1][1] = "colon"
                yield "map_key", value
            else:
                raise tokens.error("Expecting property name enclosed in double quotes")
        elif state == "colon":
            if token != ":":
                raise tokens.error("Expecting ':' delimiter")
            stack[-1][1] = "value"
        elif state == "separator":
            if token == closing:
                stack.pop()
                yield "end_" + kind, None
            elif token == ",":
                stack[-1][1] = "key" if kind == "map" else "value"
            else:
                raise tokens.error("Expecting ',' delimiter")
        elif state == "done":
            raise tokens.error("Extra data")
        elif state == "first_value" and token == "]":
            stack.pop()
            yield "end_array", None
        else:
            stack[-1][1] = "done" if kind == "root" else "separator"
            if token in ("{", "["):
                child = "map" if token == "{" else "array"
                stack.append([child, "first_key" if child == "map" else "first_value"])
                yield "start_" + child, None
            elif token == "scalar":
                yield token, value
            else:
                raise tokens.error("Expecting value")
    if stack != [["root", "done"]]:
        raise tokens.error("Unexpected end of JSON input")


def read_json_stream(handle, *, object_hook=None):
    """Build one result, retaining only bounded parser buffers and key interning.

    The returned tree and a single oversized scalar remain output-sized.
    """
    stack, keys, common_keys = [], [], {}

    def intern(key):
        if len(common_keys) < 1024:
            return common_keys.setdefault(key, key)
        return common_keys.get(key, key)

    def bulk_hook(row):
        return {intern(key): value for key, value in row.items()}

    # Speculative partial decoding must never call a caller's hook twice.
    decoder = json.JSONDecoder(object_hook=bulk_hook) if object_hook is None else None
    for event, value in json_events(handle, decoder=decoder):
        if event == "map_key":
            keys[-1] = intern(value)
        elif event.startswith("start_"):
            stack.append({} if event == "start_map" else [])
            keys.append(None)
        else:
            if event.startswith("end_"):
                value = stack.pop()
                if event == "end_map" and object_hook is not None:
                    value = object_hook(value)
                keys.pop()
            if not stack:
                result = value
            elif isinstance(stack[-1], dict):
                stack[-1][keys[-1]] = value
            else:
                stack[-1].append(value)
    return result
