# GCMAMBOX

## Building

* requires a fortran compiler, e.g.:
   ```bash
   sudo apt install gfortran
   ```

* and netcdf-fortran (easiest obtained from conda)

   ```bash
   mamba create -n gcmambox_env netcdf-fortran
   mamba activate gcmambox_env

   export NETCDF_LIB=${CONDA_PREFIX}/lib/

   export NETCDF_INCLUDE=${CONDA_PREFIX}/include/
   ```

* build
   ```bash
   mkdir build
   cd build
   cmake ../src
   make
   ```

## Running the dashboard

* using uv to install the python dependencies
   ```bash
   uv run ipython mam_dashboard.py
   ```
