# Spectral Graph Model Repository

Set up a conda environment if you do not have all the packages/compatible versions.
The list of dependencies is listed in `environment.yml`. Set-up environment using conda: 
`conda env create -f environment.yml`

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`

#### Update Notes For April 2019:
 - Functions added:
    - `forward/get_complex_laplacian.py`: adds complex laplacian to brain class based on inputs of oscillatory frequency and transmission velocity. Gives `brain.laplacian/.norm_eigenmodes/.eigenvalues/.raw_eigenvectors.

    - `forward/eigenmode.py`: functions for getting overlap scores, standardized overlap scores, purity scores dataframes.

    - `utils/functions/`:`minmax_scale_z` normalizes a Pandas dataframe's row to between 0 - 1, almost like a z-score and `highlight_max` - highlights maximum in a dataframe's series
