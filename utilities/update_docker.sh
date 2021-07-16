#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

NAME="gitlab.univ-rouen.fr:4567/coria/archer/arbalete:mcac-gitlab-ci"
VERSION=$(git log -1 --pretty=%h)

docker build "${SCRIPT_DIR}/.." -t "${NAME}-gfortran" \
  -f "${SCRIPT_DIR}/../.docker/Dockerfile.gfortran" \
  --build-arg GIT_COMMIT="${VERSION}"

docker build "${SCRIPT_DIR}/.." -t "${NAME}-intel" \
  -f "${SCRIPT_DIR}/../.docker/Dockerfile.intel" \
  --build-arg GIT_COMMIT="${VERSION}" \
  --build-arg HDF5_MAJOR=1 \
  --build-arg HDF5_MINOR=12 \
  --build-arg HDF5_PATCH=0

docker push "${NAME}-gfortran"
docker push "${NAME}-intel"
