# Spectrome Repository
[Xihe Xie](https://github.com/axiezai) & [Pablo F. Damasceno](https://github.com/pfdamasceno)

[![DOI](https://zenodo.org/badge/217634383.svg)](https://zenodo.org/badge/latestdoi/217634383)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Raj-Lab-UCSF/spectrome/master)

`Spectrome` is a combination of the words "spectrum" and "connectome". This package is the collection of codes that constructed the analysis for the publication ["Spectral graph theory of brain oscillations"](https://www.biorxiv.org/content/10.1101/589176v3). 

The spectral graph model (SGM) is a brain structure-function model that simulates brain activity power spectrum given a structural connectome. The model is linear, low-dimensional, and provides an analytical relationship between the brain's structural and functional patterns.

## Citation:
The paper is currently in pre-print stages:
https://doi.org/10.1101/589176

Code citation: Xie, X., Stanley, M. J., & Damasceno, P. F. (2019). Raj-Lab-UCSF/spectrome: Binder Ready Release. Zenodo. https://doi.org/10.5281/ZENODO.3528332

## Abstract:
The relationship between the brain’s structural wiring and the functional patterns of neural activity is of fundamental interest in computational neuroscience. We examine a hierarchical, linear graph spectral model of brain activity at mesoscopic and macroscopic scales. The model formulation yields an elegant closed-form solution for the structure-function problem, specified by the graph spectrum of the structural connectome’s Laplacian, with simple, universal rules of dynamics specified by a minimal set of global parameters. The resulting parsimonious and analytical solution stands in contrast to complex numerical simulations of high dimensional coupled non-linear neural field models. This spectral graph model accurately predicts spatial and spectral features of neural oscillatory activity across the brain and was successful in simultaneously reproducing empirically observed spatial and spectral patterns of alpha-band (8-12 Hz) and beta-band (15-30Hz) activity estimated from source localized scalp magneto-encephalography (MEG). This spectral graph model demonstrates that certain brain oscillations are emergent properties of the graph structure of the structural connectome and provides important insights towards understanding the fundamental relationship between network topology and macroscopic whole-brain dynamics.

## Set-up:

To avoid going through the hustles of setting up a `conda` environment, you can simply click on the `Binder` badge above and run the Jupyter notebooks through a docker image on the cloud.

---
Set up a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) if you do not have all the packages/compatible versions. The list of dependencies is listed in `environment.yml`. 

Set-up environment using conda, detailed instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html): 
`conda env create -f environment.yml`

If conda complains about not finding packages/libraries, make sure `conda-forge` is in the list of channels being searched by `conda`.
You may add `conda-forge` to your list of channels with the command: `conda config --add channels conda-forge`.

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate` or `conda deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`. You may have to restart your terminal to make sure this change takes effect.

## Files:
 - `../spectrome/notebooks`: contains two jupyter notebooks, `run_model` is the basic simulation of frequency spectrums with default parameters for the HCP template connectome. `SGM_ frequency_responses` looks at the eigenvalues and their frequency responses.
 - `../spectrome/data`: contains intermediate data. Includes HCP template connectome and distance matrix, optimized model parameters for the HCP connectome as well as individual subject's connectomes (N = 36).