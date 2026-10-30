#!/bin/bash
set -euo pipefail

APP_DIR=/opt/GCMAMBOX
ENV_DIR=${APP_DIR}/.venv

export UV_PROJECT_ENVIRONMENT=${ENV_DIR}
export UV_CACHE_DIR=/var/cache/gcmambox-uv

cd "$APP_DIR"

git fetch upstream --tags
git reset --hard origin/main

uv sync --frozen

restorecon -Rv "${APP_DIR}"
restorecon -Rv "${ENV_DIR}"

systemctl restart gcmambox
systemctl is-active --quiet gcmambox || {
    echo "service failed to start" >&2
    journalctl -u gcmambox -n 30 --no-pager >&2
    exit 1
}
