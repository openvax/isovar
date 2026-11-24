# Copyright (c) 2019. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from importlib import import_module

_LAZY_EXPORTS = {
    "make_isovar_arg_parser": (".main_args", "make_isovar_arg_parser"),
    "run_isovar_from_parsed_args": (".main_args", "run_isovar_from_parsed_args"),
    "protein_sequence_creator_from_args": (".protein_sequence_args", "protein_sequence_creator_from_args"),
    "read_collector_from_args": (".rna_args", "read_collector_from_args"),
}

__all__ = list(_LAZY_EXPORTS)


def __getattr__(name):
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name, __name__), attr_name)
        globals()[name] = value
        return value
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(set(globals()) | set(__all__))
