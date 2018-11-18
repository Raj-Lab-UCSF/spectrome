"""Module with generic useful functions such as to return main dir path."""

import os


def get_file_path(filename):
    """Find filename in the relative directory `../data/` .

    Args:
        filename (str): file we're looking for in the ./data/ directory.

    Returns:
        str: absolute path to file "filename" in ./data/ dir.

    """
    here_dir = os.path.dirname(os.path.realpath('__file__'))
    file_dir = os.path.join(here_dir, '../data', filename)

    return file_dir


def get_data_path():
    """Return absolute path to `/data/`."""
    here_dir = os.path.dirname(os.path.realpath('__file__'))
    data_path = os.path.join(here_dir, '../data/')

    return data_path


def get_absolute_path(relative_path='.'):
    """Return absolute path given `relative_path`.

    Args:
        relative_path (str): path relative to 'here'.

    Returns:
        str: absolute path

    """
    here_dir = os.path.dirname(os.path.realpath('__file__'))
    abs_path = os.path.join(here_dir, relative_path)

    return abs_path

def get_sibling_path(folder):
    '''returns the path of 'folder' on the same level'''
    here_dir    = os.path.dirname(os.path.realpath('__file__'))
    par_dir = os.path.abspath(os.path.join(here_dir, os.pardir))
    sibling_dir = os.path.join(par_dir, folder)
    return sibling_dir
