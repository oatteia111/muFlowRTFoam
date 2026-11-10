# gwaterFreakFoam
 gwaterFreakFoam

## Versioning

`VERSION` contains the solver version. `Allwmake`, direct `wmake`, and CMake generate
`version.hpp` from it during the build; do not edit or commit the generated header. Use
`muFlowRT -version` or `muFlowRT --version` to print the embedded version without
loading a model.

For a release, update `VERSION`, merge into `main`, and tag that commit with `v`
followed by the version (for example, `v0.1.0`). Build that tagged commit with
`MUFLOWRT_RELEASE_TAG` set to the tag; both build paths reject a mismatch. In
GitHub Actions, pass the tag name as `MUFLOWRT_RELEASE_TAG` when invoking the build.
After a release, use a development version such as `0.1.1-dev` for subsequent work.
