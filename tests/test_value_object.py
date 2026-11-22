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


from isovar.value_object import ValueObject


def test_no_fields_unless_specified():
    v = ValueObject()
    assert v._fields == ()
    assert v._values == ()


def test_default_string_repr():
    v = ValueObject()
    assert str(v) == "ValueObject()"
    assert repr(v) == "ValueObject()"


class DerivedWithoutInit(ValueObject):
    __slots__ = ["a", "b"]


def test_default_init():
    obj = DerivedWithoutInit(a=1, b=2)
    assert obj.a == 1
    assert obj.b == 2


class DerivedWithInit(ValueObject):
    __slots__ = ["a", "b"]

    def __init__(self, a, b):
        self.a = a
        self.b = b


def test_equality_checks_class():
    # two objects of different classes should not be equal
    # even if their fields are the same
    x = DerivedWithInit(a=1, b=2)
    y = DerivedWithoutInit(a=1, b=2)

    assert hash(x) == hash(y)
    assert x != y, "Expected %s != %s" % (x, y)


def test_derived_string_repr():
    x = DerivedWithInit(a=1, b=2)
    assert str(x) == "DerivedWithInit(a=1, b=2)"
    assert repr(x) == "DerivedWithInit(a=1, b=2)"
