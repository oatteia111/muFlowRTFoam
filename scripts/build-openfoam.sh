#!/usr/bin/env bash
set -euo pipefail

export OF13_ROOT=/c/OF13
export OF13_CLONE="$OF13_ROOT/OpenFOAM-13-Windows"
export OF13_THIRDPARTY="$OF13_ROOT/ThirdParty-13-Windows"

fetch_source() {
    local name=$1 ref=$2
    if [[ ! -e "$OF13_ROOT/$name/.git" ]]; then
        git init "$OF13_ROOT/$name"
        git -C "$OF13_ROOT/$name" remote add origin "https://github.com/fralarco/$name.git"
    fi
    git -C "$OF13_ROOT/$name" fetch --depth=1 origin "$ref"
    git -C "$OF13_ROOT/$name" checkout --detach FETCH_HEAD
    git -C "$OF13_ROOT/$name" rev-parse HEAD
}

# Current upstream master commits, captured 2026-10-19.
fetch_source OpenFOAM-13-Windows f57ca4c8aea92c170cb2e397f88f0dd6cfd4b81f
fetch_source ThirdParty-13-Windows 8fee8283fa502f4835a94dc2291581fc2b053d4d
{
    echo "https://github.com/fralarco/OpenFOAM-13-Windows.git $(git -C "$OF13_CLONE" rev-parse HEAD)"
    echo "https://github.com/fralarco/ThirdParty-13-Windows.git $(git -C "$OF13_THIRDPARTY" rev-parse HEAD)"
} > "$OF13_ROOT/SOURCES.txt"
# Build the serial libraries needed by muFlowRT. Keep the sibling ThirdParty
# checkout for the SDK layout, but this subset does not need Scotch or MPI.
export WM_MPLIB=Dummy FOAM_MPI=dummy
source "$OF13_CLONE/scripts/windows/env.sh"
cd "$WM_PROJECT_DIR"

make -C wmake/src
wmakeLnInclude src/OpenFOAM
src/OSspecific/MSwindows/Allwmake

# Preserve the upstream Windows bootstrap: OpenFOAM and dummy Pstream import
# each other, so a clean build needs a provisional Pstream import library.
if [[ ! -f "$FOAM_LIBBIN/dummy/libPstream.dll.a" || ! -f "$FOAM_LIBBIN/libOpenFOAM.dll" ]]; then
    (
        cd src/Pstream/dummy
        wmakeLnInclude .
        wmake objects
    )
    pstream_object=$(find "platforms/$WM_OPTIONS/obj" -name UPstream.o -print -quit)
    if [[ ! -f "$pstream_object" ]]; then
        echo "Dummy Pstream objects were not created; cannot bootstrap OpenFOAM." >&2
        exit 1
    fi
    mkdir -p "$FOAM_LIBBIN/dummy"
    dlltool --export-all-symbols -D libPstream.dll \
        -l "$FOAM_LIBBIN/dummy/libPstream.dll.a" "$(dirname "$pstream_object")"/*.o
fi

wmake libso src/OpenFOAM
wmake libso src/Pstream/dummy
for library in fileFormats surfMesh triSurface meshTools finiteVolume; do
    wmake libso "src/$library"
done

echo "Built the seven OpenFOAM libraries required by muFlowRT in $FOAM_LIBBIN."
