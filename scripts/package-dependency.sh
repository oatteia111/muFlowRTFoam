#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

case "$1" in
    phreeqcrm)
        root=/c/phreeqcrm/phreeqcrm-3.8.6-17100
        name=phreeqcrm-3.8.6-17100-ucrt64-x86_64
        ;;
    openfoam)
        root=/c/OF13
        name=openfoam-13-ucrt64-x86_64
        ;;
    *) echo "Unknown dependency: $1" >&2; exit 1 ;;
esac
dist="$(cygpath -u "$RUNNER_TEMP")/dependency-dist"
# Use a fresh directory so retries cannot mix old and new archive parts.
mkdir "$dist"
cp "$repo_dir/.github/DEPENDENCIES.md" "$root/SDK-README.md"
cp "$(cygpath -u "$RUNNER_TEMP")/dependency-abi.txt" "$root/TOOLCHAIN_ABI.txt"
{
    echo "Repository: $GITHUB_REPOSITORY"
    echo "Commit: $GITHUB_SHA"
    echo "Run: $GITHUB_RUN_ID; attempt: $GITHUB_RUN_ATTEMPT"
    cat "$root/SOURCES.txt"
    cat "$(cygpath -u "$RUNNER_TEMP")/dependency-toolchain.txt"
} > "$root/BUILD_INFO.txt"
cp "$root/BUILD_INFO.txt" "$root/SDK-README.md" "$root/TOOLCHAIN_ABI.txt" "$dist/"

# Materialize lnInclude links for extraction on Windows. Stream into parts to
# stay below GitHub's per-asset limit without an extra archive copy on disk.
tar --dereference --hard-dereference \
    --exclude=.git --exclude='*.o' --exclude='*.obj' --exclude='*.dep' --exclude='*.d' \
    --exclude='log.*' --exclude='*.log' \
    -C "$(dirname "$root")" -cf - "$(basename "$root")" \
    | zstd -T0 -3 -c \
    | split --bytes=1900M --numeric-suffixes=0 --suffix-length=4 - "$dist/$name.tar.zst.part-"
if [[ ! -f "$dist/$name.tar.zst.part-0001" ]]; then
    mv "$dist/$name.tar.zst.part-0000" "$dist/$name.tar.zst"
fi
(
    cd "$dist"
    sha256sum BUILD_INFO.txt SDK-README.md TOOLCHAIN_ABI.txt "$name".tar.zst* > SHA256SUMS
)
echo "artifact-name=deps-$name-$GITHUB_RUN_ID-$GITHUB_RUN_ATTEMPT" >> "$GITHUB_OUTPUT"
