# Copyright (c) 2016. Mount Sinai School of Medicine
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

from __future__ import print_function, division, absolute_import
from itertools import chain

class CompactObject(object):
    """
    Base class for objects which define their fields
    via __slots__ to decrease memory footporint and
    speed up field access.
    """
    __slots__ = []

    def _values(self):
        return tuple(
            getattr(self, name)
            for name in self._fields())

    @classmethod
    def _fields(cls):
        """
        Get all field names of this object by traversing
        its inheritance hierarchy in reverse order and
        concatenating the names in __slots__.
        """
        inherited_class_order = reversed(cls.__mro__)
        return tuple(
            chain.from_iterable(
                getattr(cls, '__slots__', [])
                for cls in inherited_class_order))

    def __str__(self):
        raise NotImplementedError("__str__ must be implemented")

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self._values())

    def __eq__(self, other):
        return (
            self.__class__ is other.__class__ and
            self._values() == other._values())
