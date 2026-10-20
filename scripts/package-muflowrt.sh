#!/usr/bin/env bash
set -euo pipefail

# CI packaging: the dependency SDKs and the solver must already be built.
: "${RUNNER_TEMP:?RUNNER_TEMP must point to the CI staging directory}"
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
temp=$(cygpath -u "$RUNNER_TEMP")
export OF13_ROOT="${OF13_ROOT:-/c/OF13}"
export OF13_CLONE="${OF13_CLONE:-$OF13_ROOT/OpenFOAM-13-Windows}"
export PHREEQCRM_ROOT="${PHREEQCRM_ROOT:-/c/phreeqcrm/phreeqcrm-3.8.6-17100}"
export WM_MPLIB=Dummy FOAM_MPI=dummy
source "$OF13_CLONE/scripts/windows/env.sh"

name=muflowrt-ucrt64-x86_64
dist="$temp/muflowrt-dist"
# Keep final assets separate from any previous packaging attempt.
mkdir "$dist"
# Fresh staging avoids stale binaries or archive entries on a local retry.
stage=$(mktemp -d "$temp/muflowrt-package.XXXXXX")
bundle="$stage/$name"
mkdir -p "$bundle/bin" "$bundle/OpenFOAM-13" "$bundle/database" "$bundle/licenses"

cp "$FOAM_USER_APPBIN/muFlowRT.exe" "$bundle/bin/"
for library in OpenFOAM fileFormats surfMesh triSurface meshTools finiteVolume; do
    cp "$FOAM_LIBBIN/lib$library.dll" "$bundle/bin/"
done
cp "$FOAM_LIBBIN/dummy/libPstream.dll" "$PHREEQCRM_ROOT/bin/libPhreeqcRM.dll" "$bundle/bin/"
runtime_packages=()
for dll in libgcc_s_seh-1.dll libstdc++-6.dll libwinpthread-1.dll libgomp-1.dll zlib1.dll; do
    cp "/ucrt64/bin/$dll" "$bundle/bin/"
    runtime_packages+=("$(pacman -Qqo "/ucrt64/bin/$dll")")
done

# Inspect every packaged PE file, including runtime DLLs, so transitive imports
# must resolve within this bundle or to known Windows system libraries.
for binary in "$bundle/bin/"*; do
    imports=$(objdump -p "$binary" | sed -n 's/.*DLL Name: //p' | tr -d '\r')
    while IFS= read -r dll; do
        [[ -n "$dll" ]] || continue
        case "${dll,,}" in
            api-ms-win-*|ext-ms-win-*|kernel32.dll|advapi32.dll|ntdll.dll|ucrtbase.dll|msvcrt.dll|user32.dll|ws2_32.dll|shell32.dll|ole32.dll|oleaut32.dll|bcrypt.dll|crypt32.dll|sechost.dll|rpcrt4.dll|gdi32.dll|version.dll|winmm.dll|shlwapi.dll)
                continue ;;
        esac
        # Windows DLL filenames are case-insensitive, as is the runner's volume.
        if [[ ! -f "$bundle/bin/$dll" ]]; then
            echo "Unbundled DLL import: $(basename "$binary") -> $dll" >&2
            exit 1
        fi
    done <<< "$imports"
done

cp -RL "$WM_PROJECT_DIR/etc" "$bundle/OpenFOAM-13/"
cp "$PHREEQCRM_ROOT/database/phreeqc.dat" "$bundle/database/"
cp "$repo_dir/LICENSE" "$bundle/licenses/muFlowRT-LICENSE"
cp "$WM_PROJECT_DIR/COPYING" "$bundle/licenses/OpenFOAM-COPYING"
cp "$WM_PROJECT_DIR/THIRD_PARTY_NOTICES.md" "$bundle/licenses/OpenFOAM-THIRD_PARTY_NOTICES.md"
cp "$PHREEQCRM_ROOT/doc/NOTICE.TXT" "$bundle/licenses/PhreeqcRM-NOTICE.TXT"
for package in "${runtime_packages[@]}"; do
    files=$(pacman -Qlq "$package")
    count=0
    while IFS= read -r file; do
        [[ "$file" == /ucrt64/share/licenses/* && -f "$file" ]] || continue
        relative=${file#/ucrt64/share/licenses/}
        mkdir -p "$bundle/licenses/MSYS2/$(dirname "$relative")"
        cp "$file" "$bundle/licenses/MSYS2/$relative"
        count=$((count + 1))
    done <<< "$files"
    if (( count == 0 )); then
        echo "No installed license files found for runtime package $package" >&2
        exit 1
    fi
done

cp "$repo_dir/scripts/run-muflowrt.ps1" "$bundle/run-muflowrt.ps1"
cat > "$bundle/README.txt" <<'EOF'
Unsigned muFlowRT Windows UCRT64 development build (serial OpenFOAM, OpenMP PhreeqcRM).

Extract the entire ZIP. In an existing PowerShell terminal, check startup with:
    & 'C:\path\to\muflowrt-ucrt64-x86_64\run-muflowrt.ps1' -SolverArguments '-help'
Run a prepared case with:
    & 'C:\path\to\muflowrt-ucrt64-x86_64\run-muflowrt.ps1' -CaseDirectory 'C:\models\my-case'
Use -SolverArguments @('-option', 'value') for additional solver arguments.
If PowerShell blocks the unsigned script, allow it for this terminal only:
    Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass

The launcher forwards solver arguments and exit status, then restores your
environment and working directory afterward.
The case directory must contain its intended phreeqc.dat. A reference copy is
included in database/; copy it into a case only if it is the database you intend
to use. No case data is overwritten automatically.

The launcher uses bundled DLLs and Windows system directories. Normal serial
cases do not require MSYS2. Runtime compilation (#codeStream/#calc), additional
OpenFOAM plugins, MPI, and OpenFOAM command-line utilities require a separate SDK
and toolchain; they are not supplied by this runtime bundle.

BUILD_INFO.txt records sources and toolchains; SHA256SUMS covers the bundle files.
Numerical behavior has not been validated by the build workflow.
EOF
{
    echo "Repository: $GITHUB_REPOSITORY"
    echo "Commit: $GITHUB_SHA"
    echo "Run: $GITHUB_RUN_ID; attempt: $GITHUB_RUN_ATTEMPT"
    echo 'muFlowRT: wmake, UCRT64 x86_64, DP/Int32/Opt, dummy Pstream, -fopenmp'
    cat "$temp/dependency-toolchain.txt"
    for dependency in phreeqcrm openfoam; do
        printf '\nDependency: %s\n' "$dependency"
        cat "$temp/sdk-$dependency/BUILD_INFO.txt"
    done
} > "$bundle/BUILD_INFO.txt"
cp "$temp/dependency-abi.txt" "$bundle/TOOLCHAIN_ABI.txt"
(
    cd "$bundle"
    find . -type f ! -name SHA256SUMS -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum --check --strict SHA256SUMS
)
(
    cd "$stage"
    zip -qr "$stage/$name.zip" "$name"
)

# Stage uniquely named SDK archives together; their metadata is already inside
# each archive and is included in the combined BUILD_INFO above.
shopt -s nullglob
for dependency in phreeqcrm openfoam; do
    archives=("$temp/sdk-$dependency/"*.tar.zst*)
    if (( ${#archives[@]} == 0 )); then
        echo "Missing $dependency SDK release assets" >&2
        exit 1
    fi
    cp "${archives[@]}" "$dist/"
done
cp "$bundle/BUILD_INFO.txt" "$bundle/TOOLCHAIN_ABI.txt" "$dist/"
cp "$repo_dir/.github/DEPENDENCIES.md" "$dist/SDK-README.md"
# Keep every release asset below GitHub's limit, including a future larger ZIP.
if (( $(stat -c %s "$stage/$name.zip") > 1900 * 1024 * 1024 )); then
    split --bytes=1900M --numeric-suffixes=0 --suffix-length=4 \
        "$stage/$name.zip" "$dist/$name.zip.part-"
else
    cp "$stage/$name.zip" "$dist/"
fi
(
    cd "$dist"
    assets=(BUILD_INFO.txt TOOLCHAIN_ABI.txt SDK-README.md *.tar.zst* *.zip*)
    sha256sum "${assets[@]}" > SHA256SUMS
    sha256sum --check --strict SHA256SUMS
)
echo "artifact-name=muflowrt-ucrt64-x86_64-$GITHUB_RUN_ID-$GITHUB_RUN_ATTEMPT" >> "$GITHUB_OUTPUT"
echo "Release assets: $dist"
