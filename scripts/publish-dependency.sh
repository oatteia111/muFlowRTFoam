#!/usr/bin/env bash
set -euo pipefail

tag=${1:?Usage: publish-dependency.sh TAG STAGING_DIRECTORY}
dist=${2:?Usage: publish-dependency.sh TAG STAGING_DIRECTORY}
# Git Bash accepts Windows paths from Actions, but normalize them before using
# Unix utilities. This also keeps the helper usable on a Linux runner.
if command -v cygpath >/dev/null; then dist=$(cygpath -u "$dist"); fi
cd "$dist"
dist=$PWD
sha256sum --check --strict SHA256SUMS

if ! gh release view "$tag" --repo "$GITHUB_REPOSITORY" >/dev/null 2>&1; then
    gh release create "$tag" --repo "$GITHUB_REPOSITORY" \
        --target "$GITHUB_SHA" --title "$tag" --draft --prerelease --latest=false \
        --notes "Unsigned Windows UCRT64 muFlowRT runtime and dependency SDKs. See SDK-README.md and BUILD_INFO.txt in the assets."
fi

# A failed publication can be rerun against the retained build artifact.
# Reuse identical uploads; never replace a different asset with the same name.
existing=$(gh release view "$tag" --repo "$GITHUB_REPOSITORY" --json assets --jq '.assets[].name')
existing_dir=$(mktemp -d)
for asset in "$dist"/*; do
    name=$(basename "$asset")
    if grep -Fxq -- "$name" <<< "$existing"; then
        gh release download "$tag" --repo "$GITHUB_REPOSITORY" \
            --pattern "$name" --dir "$existing_dir" --clobber
        if ! cmp -s "$asset" "$existing_dir/$name"; then
            echo "Refusing to replace conflicting release asset: $name" >&2
            exit 1
        fi
    else
        gh release upload "$tag" "$asset" --repo "$GITHUB_REPOSITORY"
    fi
done
gh release edit "$tag" --repo "$GITHUB_REPOSITORY" --draft=false --prerelease --latest=false
