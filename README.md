# GCMAMBOX

gcmambox

## Building

requires
* cmake
* fortran compiler
* netcdf-fortran-devel

### RHEL

* Available from package manager
   ```bash
   dnf in cmake gcc-gfortran netcdf-devel netcdf-fortran netcdf-fortran-devel
   ```

### Ubuntu

* install fortran compiler, via package manager:
   ```bash
   sudo apt install gfortran
   ```
* and netcdf-fortran (via conda)
   ```bash
   mamba create -n gcmambox_env netcdf-fortran
   mamba activate gcmambox_env

   export NETCDF_LIB=${CONDA_PREFIX}/lib/
   export NETCDF_INCLUDE=${CONDA_PREFIX}/include/
   ```

* build
   ```bash
   cd src/
   cmake -S . -B build
   cmake --build build
   ```

## Running the dashboard

* using uv to install the python dependencies
   ```bash
   uv run gunicorn dashboard.mam_dashboard:app
   ```


## Deployment

1. Add new user

   ```bash
   useradd --system --no-create-home --shell /sbin/nologin gcmambox
   ```

2. Enable port on SELinux

   ```bash
   # ensure port is not used
   semanage port -l | grep -w 8003
   # add it
   semanage port -a -t http_port_t -p tcp 8003
   ```

3. Clone
   ```bash
   git clone https://github.com/fsolmon/GCMAMBOX.git --recurse-submodules
   ```

4. Build gcmambox
   ```bash
   cd GCMAMBOX/src
   cmake -S . -B build
   cmake --build build
   ```

5. Setup apach config

   ```apacheconf
   <VirtualHost *:80>

      ServerName gcmambox.iac.ethz.ch
      ServerAdmin iac-linux@env.ethz.ch

      <Proxy *>
         Require all granted
      </Proxy>

      ProxyPass / http://127.0.0.1:8003/
      ProxyPassReverse / http://127.0.0.1:8003/
      RequestHeader set "X-Forwarded-Proto" expr=%{REQUEST_SCHEME}

      UseCanonicalName on
      ProxyPreserveHost on

      ErrorLog "logs/gcmambox.iac.ethz.ch-error.log"
      LogLevel warn
      CustomLog "logs/gcmambox.iac.ethz.ch-access.log" combined

   </VirtualHost>
   ```

6. Certificate (Let's Encrypt)

   ```bash
   certbot # and follow the instructions
   ```

7. Start systemd service

   ```bash
   ln -s /opt/GCMAMBOX/service/gcmambox.service /etc/systemd/system/

   systemctl daemon-reload
   systemctl enable gcmambox.service
   systemctl start gcmambox.service
   systemctl status gcmambox.service
   ```

## Update

 - The code and deps are updated with
   ```bash
   ./deploy.sh
   ```
⚠️ **Warning**: the submodule and executable is currently not autmatically updated!

