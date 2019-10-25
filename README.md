# Spectrome Repository
[Xihe Xie](https://github.com/axiezai) & [Pablo F. Damasceno](https://github.com/pfdamasceno)

`Spectrome` is a combination of the words "spectrum" and "connectome". This package is the collection of codes that constructed the analysis for the publication ["Spectral graph theory of brain oscillations"](https://www.biorxiv.org/content/10.1101/589176v3). 

The spectral graph model (SGM) is a brain structure-function model that simulates brain activity power spectrum given a structural connectome. The model is linear, low-dimensional, and provides an analytical relationship between the brain's structural and functional patterns.

## Citation:
The paper is currently in pre-print stages:
https://doi.org/10.1101/589176

Code citation:

## Set-up:
Set up a conda environment if you do not have all the packages/compatible versions.
The list of dependencies is listed in `environment.yml`. 

Set-up environment using conda: 
`conda env create -f environment.yml`

If conda complains about not finding packages/libraries, make sure `conda-forge` is in the list of channels being searched by `conda`.
You may add `conda-forge` to your list of channels with the command: `conda config --add channels conda-forge`.

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate` or `conda deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`
