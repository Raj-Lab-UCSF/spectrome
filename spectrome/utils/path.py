"""Module with generic useful functions such as to return main dir path, and
writing hdf5 files."""

import os
import deepdish as dd

from pathlib import Path


def get_file_path(filename):
    """Find filename in the relative directory `../data/` .

    Args:
        filename (str): file we're looking for in the ./data/ directory.

    Returns:
        str: absolute path to file "filename" in ./data/ dir.

    """
    root_dir = Path(__file__).parent.parent
    file_dir = os.path.join(str(root_dir), "data", filename)

    return file_dir


def get_data_path():
    """Return absolute path to `/data/`."""
    root_path = Path(__file__).parent.parent
    data_path = os.path.join(str(root_path), "data")
    return data_path


def get_absolute_path(relative_path="."):
    """Return absolute path given `relative_path`.

    Args:
        relative_path (str): path relative to 'here'.

    Returns:
        str: absolute path

    """
    here_dir = os.path.dirname(os.path.realpath("__file__"))
    abs_path = os.path.join(str(here_dir), relative_path)

    return abs_path


def get_sibling_path(folder):
    """returns the path of 'folder' on the same level"""
    root_dir = Path(__file__).parent.parent
    sibling_dir = os.path.join(str(root_dir), folder)
    return sibling_dir


def get_root_path():
    root_path = Path(__file__).parent.parent
    return root_path


def save_hdf5(path, dict):
    """Save out a dictionary/numpy array to HDF5 format using deepdish package.

    Args:
        path (type): full path including filename of intended output.
        dict (type): dictionary/numpy array to be saved.

    Returns:
        type: Description of returned object.

    """

    dd.io.save(path, dict)


def read_hdf5(path):
    """Read in dictionary/numpy array from HDF5 format using deepdish package.

    Args:
        path (type): full path including filename of input.

    Returns:
        type: dictionary of data.

    """

    dict = dd.io.load(path)
    return dict


def walk_tree(datapath):
    """Return list of directories in the passed folder.

    Args:
        datapath (type): folder of interest.

    Returns:
        type: list of directories in the passed folder.

    """
    directories = []
    for (path, dirs, files) in os.walk(datapath):
        directories.append(dirs)
    return directories[0]
