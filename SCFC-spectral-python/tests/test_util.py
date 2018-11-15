"""Test suite for the util module."""

import os
from nose.tools import with_setup
from src import util

def f_setup():
    f = open('./data/test.txt', 'w')


def f_teardown():
    files = os.listdir('./data/')
    if 'test.txt' in files:
        os.remove('./data/test.txt')

def test_get_file_path():
    f_setup()
    FILEPATH = util.get_file_path('test.txt')

    assert os.path.exists(FILEPATH) == True

    # f_teardown()

def test_get_data_path():
    f_setup()
    DATAPATH = util.get_data_path()
    FILEPATH = os.path.join(DATAPATH, 'test.txt')

    assert os.path.exists(FILEPATH) == True



util.get_data_path()

test_get_file_path()
