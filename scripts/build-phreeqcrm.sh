#!/usr/bin/env bash
set -euo pipefail

version=3.8.6-17100
archive="phreeqcrm-$version.tar.gz"
source_dir=/c/phreeqcrm/source
build_dir=/c/phreeqcrm/build
install_dir="/c/phreeqcrm/phreeqcrm-$version"
mkdir -p "$source_dir"
curl --fail --location --retry 3 \
    "https://github.com/phreeqc-dev/phreeqcrm/releases/download/v3.8.6/$archive" \
    --output "$source_dir/$archive"
echo "d2e532a065b17ab07104c7f445c823aadafab10b263b613b1fec7e5fb50b046a  $source_dir/$archive" \
    | sha256sum --check --strict -
tar -xzf "$source_dir/$archive" -C "$source_dir"

# Retry OpenMP detection if a previous configure ran before libgomp was installed.
cmake -U 'OpenMP*' -S "$source_dir/phreeqcrm-$version" -B "$build_dir" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$(cygpath -m "$install_dir")" \
    -DBUILD_SHARED_LIBS=ON \
    -DPHREEQCRM_BUILD_MPI=OFF \
    -DPHREEQCRM_BUILD_PYTHON=OFF \
    -DPHREEQCRM_FORTRAN_TESTING=OFF \
    -DPHREEQCRM_WITH_YAML_CPP=OFF \
    -DPHREEQCRM_USE_ZLIB=OFF \
    -DPHREEQCRM_DISABLE_OPENMP=OFF \
    -DCMAKE_REQUIRE_FIND_PACKAGE_OpenMP=TRUE
cmake --build "$build_dir" --target PhreeqcRM --parallel "$(nproc)"
cmake --install "$build_dir"
mkdir -p "$install_dir/sources"
cp "$source_dir/$archive" "$install_dir/sources/"
cp "$build_dir/CMakeCache.txt" "$install_dir/BUILD_OPTIONS.txt"
printf '%s\n' \
    "https://github.com/phreeqc-dev/phreeqcrm/releases/download/v3.8.6/$archive" \
    'SHA256=d2e532a065b17ab07104c7f445c823aadafab10b263b613b1fec7e5fb50b046a' \
    > "$install_dir/SOURCES.txt"
