# gwaterFreakFoam
 gwaterFreakFoam

## Windows builds

See [the Windows build guide](.github/DEPENDENCIES.md) for local builds with
`scripts/build-muflowrt.sh`, the consolidated PhreeqcRM/OpenFOAM/muFlowRT workflow,
and unsigned runtime and SDK releases.

## Chemistry scheduling and threads

The normal chemistry cadence follows `reactionSteps`. At the final transport
endpoint, any elapsed time since the last chemistry calculation is reacted before
observations, budgets, and results are written. Chemistry already scheduled at
that endpoint runs once. Intermediate output times do not change the cadence.

OpenMP builds request `omp_get_max_threads()` threads from PhreeqcRM, honoring
`OMP_NUM_THREADS` in the solver's launch environment. Builds without OpenMP request
one thread. Startup logs the requested count and the count configured by PhreeqcRM.
This is a solver-construction setting, not a PHREEQC `.pqi` option. Arti3D can inherit
the environment setting when launched; no GUI setting is added by this change.
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
