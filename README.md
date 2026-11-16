# gwaterFreakFoam
 gwaterFreakFoam

## Chemistry scheduling and threads

The normal chemistry cadence follows `reactionSteps`. At the final transport
endpoint, any elapsed time since the last chemistry calculation is reacted before
observations, budgets, and results are written. Chemistry already scheduled at
that endpoint runs once. Intermediate output times do not change the cadence.

OpenMP builds request `omp_get_max_threads()` threads from PhreeqcRM, honoring
`OMP_NUM_THREADS` in the solver's launch environment. Builds without OpenMP request
one thread. Startup logs the requested count and the count configured by PhreeqcRM.
This is a solver-construction setting, not a PHREEQC `.pqi` option. 

## Versioning

`VERSION` contains the solver version. `Allwmake`, direct `wmake`, and CMake generate
`version.hpp` from it during the build; do not edit or commit the generated header. Use
`muFlowRT -version` or `muFlowRT --version` to print the embedded version without
loading a model.

The current solver version is `0.2.0`. Update `VERSION` when the solver version
changes. Header generation reads this file independently of Git tags or any
external packaging workflow.

`solvers/version.hpp.in` also declares `MUFLOWRT_TIME_SCHEDULE_VERSION`, currently
`2`. This capability is independent of the solver release version and is included
in the generated header. The solver supports legacy schedule version `1` and the
declared version; cases without `timeScheduleVersion` still default to `1`.
