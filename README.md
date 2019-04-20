# Spectral Graph Model Repository

Set up a conda environment if you do not have all the packages/compatible versions.
The list of dependencies is listed in `environment.yml`. Set-up environment using conda: 
`conda env create -f environment.yml`

The default name of the environment is `spectral`, activate the environment with `source activate spectral`, and deactivate with `source deactivate`.

#### Update Notes For April 2019:
 - Functions added:
    - `forward/get_complex_laplacian.py`: adds complex laplacian to brain class based on inputs of oscillatory frequency and transmission velocity. Gives `brain.laplacian/.norm_eigenmodes/.eigenvalues/.raw_eigenvectors.