#!/bin/bash
set -euo pipefail

APP_DIR=/opt/GCMAMBOX
ENV_DIR=${APP_DIR}/.venv

if [ ! -x "${ENV_DIR}/bin/gunicorn" ]; then
    echo "venv not usable at $ENV_DIR — run deploy first" >&2
    exit 1
fi

# if [ -d "$APP_DIR/.venv" ]; then
#    echo "venv not usable at $APP_DIR/.venv — remove" >&2
#    exit 1
#fi

cd ${APP_DIR}

# not sure if neccesary
source ${ENV_DIR}/bin/activate

# set dir for generated_files
# export CACHE_DIR=/var/cache/caps

export MPLCONFIGDIR=/tmp/mpl

exec "${ENV_DIR}/bin/gunicorn" dashboard.mam_dashboard:app --bind 127.0.0.1:8003 --workers 4 --forwarded-allow-ips 127.0.0.1 --control-socket /tmp/.gunicorn/gunicorn.ctl
