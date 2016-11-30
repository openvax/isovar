from isovar.value_object import ValueObject
from nose.tools import eq_

def test_no_fields_unless_specified():
    v = ValueObject()
    eq_(v._fields, ())
    eq_(v._values, ())

def test_default_string_repr():
    v = ValueObject()
    eq_(str(v), "ValueObject()")
    eq_(repr(v), "ValueObject()")

class DerivedWithoutInit(ValueObject):
    __slots__ = ["a", "b"]

def test_default_init():
    obj = DerivedWithoutInit(a=1, b=2)
    eq_(obj.a, 1)
    eq_(obj.b, 2)

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

    eq_(hash(x), hash(y))
    assert x != y, "Expected %s != %s" % (x, y)

def test_derived_string_repr():
    x = DerivedWithInit(a=1, b=2)
    eq_(str(x), "DerivedWithInit(a=1, b=2)")
    eq_(repr(x), "DerivedWithInit(a=1, b=2)")
