from itertools import islice

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk


def batcher(iterable, batch_size):
    iterator = iter(iterable)
    while batch := list(islice(iterator, batch_size)):
        yield batch