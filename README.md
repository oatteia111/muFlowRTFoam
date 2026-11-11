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
