#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
if [[ ${MSYSTEM:-} != UCRT64 ]]; then
    echo 'Run this script from an MSYS2 UCRT64 terminal.' >&2
    exit 1
fi

export OF13_ROOT="${OF13_ROOT:-/c/OF13}"
export OF13_CLONE="${OF13_CLONE:-$OF13_ROOT/OpenFOAM-13-Windows}"
export PHREEQCRM_ROOT="${PHREEQCRM_ROOT:-/c/phreeqcrm/phreeqcrm-3.8.6-17100}"
if [[ ! -f "$OF13_CLONE/scripts/windows/env.sh" ]]; then
    echo "Missing OpenFOAM environment: $OF13_CLONE/scripts/windows/env.sh" >&2
    echo 'Build or extract the OpenFOAM SDK first, or set OF13_ROOT/OF13_CLONE.' >&2
    exit 1
fi
export WM_MPLIB=Dummy FOAM_MPI=dummy
source "$OF13_CLONE/scripts/windows/env.sh"
export LIB_DIR="$WM_PROJECT_DIR/src"
export LIB_HEADER_DIRS="$WM_PROJECT_DIR/src"

required=(
    "$WM_PROJECT_DIR/etc/controlDict"
    "$WM_DIR/wmake"
    "$WM_DIR/platforms/${WM_ARCH}${WM_COMPILER}/wmkdep.exe"
    "$WM_DIR/platforms/${WM_ARCH}${WM_COMPILER}/dirToString.exe"
    "$LIB_DIR/OpenFOAM/lnInclude/argList.H"
    "$LIB_DIR/finiteVolume/lnInclude/fvMesh.H"
    "$LIB_DIR/meshTools/lnInclude/cellSet.H"
    "$PHREEQCRM_ROOT/include/PhreeqcRM.h"
    "$PHREEQCRM_ROOT/include/IPhreeqc.hpp"
    "$PHREEQCRM_ROOT/include/irm_dll_export.h"
    "$PHREEQCRM_ROOT/lib/libPhreeqcRM.dll.a"
    "$PHREEQCRM_ROOT/bin/libPhreeqcRM.dll"
)
for library in OpenFOAM dummy/Pstream fileFormats surfMesh triSurface meshTools finiteVolume; do
    directory="$FOAM_LIBBIN"
    if [[ $library == dummy/* ]]; then directory+=/dummy; fi
    required+=("$directory/lib${library##*/}.dll" "$directory/lib${library##*/}.dll.a")
done
missing=0
for file in "${required[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "Missing dependency file: $file" >&2
        missing=1
    fi
done
if (( missing )); then
    echo 'Build or restore the dependency SDKs before building muFlowRT.' >&2
    exit 1
fi
for command in g++ make wmake; do
    command -v "$command" >/dev/null || { echo "Missing build tool: $command" >&2; exit 1; }
done
export PATH="$PHREEQCRM_ROOT/bin:$PATH"

# The repository-level Allwmake still references obsolete library directories.
cd "$repo_dir/solvers"
wmake
executable="$FOAM_USER_APPBIN/muFlowRT.exe"
[[ -f "$executable" ]] || { echo "wmake did not produce $executable" >&2; exit 1; }
printf '\nBuilt muFlowRT: %s\n' "$(cygpath -m "$executable")"
