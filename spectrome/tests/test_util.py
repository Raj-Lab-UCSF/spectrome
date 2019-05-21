"""Test suite for the util module."""

import os
from nose.tools import with_setup
from ..utils import path


def f_setup():
    """Creates test file in `./data/` dir."""
    open('../data/test.txt', 'w')


def f_teardown():
    """Removes test file from `./data/` dir."""
    files = os.listdir('../data/')
    if 'test.txt' in files:
        os.remove('../data/test.txt')


@with_setup(setup=f_setup, teardown=f_teardown)
def test_get_file_path():
    """Checks if path.get_file_path(file) returns `./data/file`."""
    f_setup()
    file_path = path.get_file_path('test.txt')
    assert os.path.exists(file_path) is True
    f_teardown()


@with_setup(setup=f_setup, teardown=f_teardown)
def test_get_data_path():
    """Checks if path.get_data_path(file) returns `./data/`."""
    f_setup()
    data_path = path.get_data_path()
    file_path = os.path.join(data_path, 'test.txt')
    assert os.path.exists(file_path) is True
    f_teardown()
