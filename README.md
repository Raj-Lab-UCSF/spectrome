# Spectrome Repository
[Xihe Xie](https://github.com/axiezai), [Megan J. Stanley](https://github.com/megstanley), & [Pablo F. Damasceno](https://github.com/pfdamasceno)

[![DOI](https://zenodo.org/badge/217634383.svg)](https://zenodo.org/badge/latestdoi/217634383)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Raj-Lab-UCSF/spectrome/master)

`Spectrome` is a combination of the words "spectrum" and "connectome". This package is the collection of codes that constructed the analysis for the publication ["Spectral graph theory of brain oscillations"](https://www.biorxiv.org/content/10.1101/589176v3).

The spectral graph model (SGM) is a brain structure-function model that simulates brain activity power spectrum given a structural connectome. The model is linear, low-dimensional, and provides an analytical relationship between the brain's structural and functional patterns.

## Citation:
The code in this repository is used for the analysis as shown in: 
 - Ashish Raj, Chang Cai, Xihe Xie, Eva Palacios, Julia Owen, Pratik Mukherjee, and Srikantan Nagarajan. “Spectral Graph Theory of Brain Oscillations.” Human Brain Mappin. https://doi.org/10.1002/hbm.24991.

 - Xihe Xie, Chang Cai, Pablo Damasceno, Srikantan Nagarajan, and Ashish Raj. "Emergence of Canonical Functional Networks from the Structural Connectome." NeuroImage. https://doi.org/10.1016/j.neuroimage.2021.118190.


Code citation: Xihe Xie, Megan J. Stanley, & Pablo F. Damasceno. (2019, November 7). Raj-Lab-UCSF/spectrome: Spectral Graph Model of Connectomes (Version 0.15). Zenodo. http://doi.org/10.5281/zenodo.3532497

## Abstract:
The relationship between the brain’s structural wiring and the functional patterns of neural activity is of fundamental interest in computational neuroscience. We examine a hierarchical, linear graph spectral model of brain activity at mesoscopic and macroscopic scales. The model formulation yields an elegant closed-form solution for the structure-function problem, specified by the graph spectrum of the structural connectome’s Laplacian, with simple, universal rules of dynamics specified by a minimal set of global parameters. The resulting parsimonious and analytical solution stands in contrast to complex numerical simulations of high dimensional coupled non-linear neural field models. This spectral graph model accurately predicts spatial and spectral features of neural oscillatory activity across the brain and was successful in simultaneously reproducing empirically observed spatial and spectral patterns of alpha-band (8-12 Hz) and beta-band (15-30Hz) activity estimated from source localized scalp magneto-encephalography (MEG). This spectral graph model demonstrates that certain brain oscillations are emergent properties of the graph structure of the structural connectome and provides important insights towards understanding the fundamental relationship between network topology and macroscopic whole-brain dynamics.

## Set-up:

To avoid going through the hustles of setting up a `conda` environment, you can simply click on the `Binder` badge above and run the Jupyter notebooks through a docker image on the cloud.

---
First clone the environment to your computer, either download this repo as a `.zip` file or `git clone https://github.com/Raj-Lab-UCSF/spectrome.git`.

Set up a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) if you do not have all the packages/compatible versions. The list of dependencies is listed in `environment.yml`.

Set-up environment using conda, detailed instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Or after cloning this repository, go to the repo by typing `cd spectrome` and then typing:
`conda env create -f environment.yml`

If conda complains about not finding packages/libraries, make sure `conda-forge` is in the list of channels being searched by `conda`.
You may add `conda-forge` to your list of channels with the command: `conda config --add channels conda-forge`.

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate` or `conda deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`. You may have to restart your terminal to make sure this change takes effect.

After completing the set-up for conda environment and `spectrome` path, you may go to the `spectrome` folder and type `jupyter notebook` or `jupyter lab` in your terminal to run the Jupyter notebooks.

## Files:
 - `../spectrome/notebooks`: contains three jupyter notebooks, `run_model_example.ipynb` is the basic simulation of frequency spectrums with default parameters for the HCP template connectome. `SGM_frequency_responses` looks at the eigenvalues and their frequency responses, and `reproduce_MEG_spectra.ipynb` shows an example of the spectral graph model recapitulating empirical MEG spectra with optimized parameters.

 - `../spectrome/data`: contains intermediate data.
    - `mean80_fibercount/length.csv`: HCP template connectome and distance matrix.
    - The following `.mat` files come from Matlab structs, `..['output']['param']` are the fields for the optimized parameters, other fields are other outputs from the optimization algorith.
        - `SCFC_opparam_individual.mat`: individual subject's connectivity matrices (N = 36).
        - `SCFC_opparam_HCP.mat`: optimized model parameters for the HCP connectome.
    - `freqMEGdata_8002-101.h5`: one example source localized MEG spectrum. Dictionary with region names attached to each column (keys).
