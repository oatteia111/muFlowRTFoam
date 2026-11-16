#!/bin/sh
cd "${0%/*}/.." || exit 1

version=$(tr -d '\r\n' < VERSION) || exit 1
if [ -z "$version" ]; then
    echo "VERSION must not be empty" >&2
    exit 1
fi
if [ -n "${MUFLOWRT_RELEASE_TAG:-}" ] && [ "$MUFLOWRT_RELEASE_TAG" != "v$version" ]; then
    echo "Release tag $MUFLOWRT_RELEASE_TAG does not match VERSION ($version)" >&2
    exit 1
fi
header=$(sed "s/@MUFLOWRT_VERSION@/$version/g" solvers/version.hpp.in) || exit 1
if [ ! -f solvers/version.hpp ] || [ "$header" != "$(cat solvers/version.hpp)" ]; then
    printf '%s\n' "$header" > solvers/version.hpp || exit 1
fi
