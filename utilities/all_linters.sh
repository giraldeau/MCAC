#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "${SCRIPT_DIR}/.." || exit 1

echo "** C++ **"
clang-format -i src/**/*.cpp include/**/*.hpp
cpplint --recursive src include
cppcheck --enable=all --inconclusive --error-exitcode=2 -I include src

echo "** python **"
isort .
black .
pflake8
mypy pymcac
pylint pymcac

echo "** shell scripts **"
shfmt -i 2 -ci -w -d utilities
shellcheck utilities/*.sh

echo "** cmake **"
cmake-format --check CMakeLists.txt cmake/**/*.cmake
cmake-lint CMakeLists.txt cmake/**/*.cmake

echo "** Dockerfiles **"
hadolint .docker/Dockerfile.*
