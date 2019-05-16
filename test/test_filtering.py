from isovar.filtering import apply_filters

class MockIsovarResult(object):
    """
    Mock object which has whatever properties we give it but can only
    be cloned with new filter_values passed into to clone_with_updates by
    apply_filters.
    """
    def __init__(self, filter_values={}, **kwargs):
        self._init_args = kwargs
        self.filter_values = filter_values
        for (k, v) in kwargs.items():
            setattr(self, k, v)

    def clone_with_updates(self, filter_values):
        return MockIsovarResult(filter_values=filter_values, **self._init_args)

def test_apply_filters_min_pass():
    obj = MockIsovarResult(x=1)
    new_obj = apply_filters(obj, filter_thresholds={"min_x": 1})
    assert new_obj.filter_values["min_x"]

def test_apply_filters_min_fail():
    obj = MockIsovarResult(x=1)
    new_obj = apply_filters(obj, filter_thresholds={"min_x": 2})
    assert not new_obj.filter_values["min_x"]

def test_apply_filters_max_pass():
    obj = MockIsovarResult(x=1)
    new_obj = apply_filters(obj, filter_thresholds={"max_x": 1})
    assert new_obj.filter_values["max_x"]

def test_apply_filters_max_fail():
    obj = MockIsovarResult(x=2)
    new_obj = apply_filters(obj, filter_thresholds={"max_x": 1})
    assert not new_obj.filter_values["max_x"]

def test_apply_filters_bool_pass():
    obj = MockIsovarResult(x=True)
    new_obj = apply_filters(obj, filter_flags=["x"])
    assert new_obj.filter_values["x"]

def test_apply_filters_bool_fail():
    obj = MockIsovarResult(x=False)
    new_obj = apply_filters(obj, filter_flags=["x"])
    assert not new_obj.filter_values["x"]

def test_apply_filters_negated_bool_pass():
    obj = MockIsovarResult(x=False)
    new_obj = apply_filters(obj, filter_flags=["not_x"])
    assert new_obj.filter_values["not_x"]

def test_apply_filters_negated_bool_fail():
    obj = MockIsovarResult(x=True)
    new_obj = apply_filters(obj, filter_flags=["not_x"])
    assert not new_obj.filter_values["not_x"]