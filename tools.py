"""
"""

class Counter:
    """
    Return integers in order at each call.
    """

    def __init__(self, initial=0):
        """
        Parameters
        ----------
        initial : int
            Starting point of the counter. (default: 0)
        """

        assert type(initial) == int
        self.counter = initial - 1

    def __call__(self):

        self.counter += 1
        return self.counter

class MultiIntKeyDict(object):
    """
    Dictionary-like object in which multiple integer keys are associated to the
    same value. (see https://stackoverflow.com/questions/11105115/)
    """

    def __init__(self, **kwargs):
        """
        Initialises with (possibly empty) dictionary.
        """

        self._keys = {}	# secondary dictionary which associates multiple keys to an unique key of self._data
        self._data = {}	# values with unique identification

        for k, v in kwargs.items():
            self[k] = v

    def __getitem__(self, key):

        return self._data[self._keys[key]]

    def __setitem__(self, k, v):

        try:

            self._data[self._keys[k]] = v

        except KeyError:

            # check that key in an integer or a tuple of integers
            try:
                assert isinstance(k, int)
                keys = (k,)
            except AssertionError:
                assert isinstance(k, tuple)
                keys = k
                for key in keys:
                    assert isinstance(key, int)

            # unique indices for new value
            try:
                uniqueIndex = max(self._data) + 1
            except ValueError:
                uniqueIndex = 0

            # associate all keys to value
            self._data[uniqueIndex] = v
            for key in keys:
                self._keys[key] = uniqueIndex

    def __delitem__(self, key):

        uniqueIndex = self._keys[key]
        del self._keys[key]
        if not(uniqueIndex in self._keys.values()):
            del self._data[uniqueIndex]

    def __contains__(self, key):

        return self._keys.__contains__(key)

    def __iter__(self):

        return self._keys.__iter__()

    def __next__(self):

        return self._keys.__next__()

    def values(self):

        return self._data.values()

    def add_key(self, existingKey, *newKeys):
        """
        Associates a new key to an existing key.

        Parameters
        ----------
        existingKey : *
            Key already associated to a value.
        newKeys : *
            New keys to associate to the same value.
        """

        uniqueIndex = self._keys[existingKey]
        for newKey in newKeys:
            self._keys[newKey] = uniqueIndex

