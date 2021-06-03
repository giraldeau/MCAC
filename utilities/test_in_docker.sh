#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

GIT_COMMIT=$(git log -1 --pretty=%h)
COMPILATOR=ifort
#COMPILATOR=gfortran

bash "${SCRIPT_DIR}/update_docker.sh"
docker build "${SCRIPT_DIR}/.." -t arbalete:mcac-local \
  -f "${SCRIPT_DIR}/../.docker/Dockerfile.local" \
  --build-arg GIT_COMMIT="${GIT_COMMIT}" \
  --build-arg COMPILATOR="${COMPILATOR}" \
  --no-cache
