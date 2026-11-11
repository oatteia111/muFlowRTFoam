# Time scheduling checks

From the repository root in the existing MSYS2 UCRT64 toolchain, compile the
small checks against the same helper used by the solver:

```sh
g++ -std=c++17 -Wall -Wextra -pedantic solvers/tests/timeScheduleTest.C -o /tmp/muflowrt-time-test.exe
/tmp/muflowrt-time-test.exe
```

These checks cover exact event clipping, startup, single/final output, close
fractional events, a restart's next target, malformed schedules, precise and
legacy boundary files, safe EOF, normal chemistry cadence, and completion of
the final chemistry interval without duplicate calls. They do not compile or
validate OpenFOAM, PHREEQC, equation assembly, restart restoration, or the complete
solver loop.

Version 2 cases retain the four-float binary cell payload and provide one
double-precision timestamp per frame in `<boundary-file>.times`, starting at
zero. The final real frame is held indefinitely. `controlDict` must contain
`timeScheduleVersion 2;` and the exact end time, which must match the final
requested output. Unversioned cases use their original float32 timestamps and
configured end time; existing cases must be re-exported to remove padded ends.

After a normal solver build, use disposable, freshly exported copies of
`flow3d_onlythesetimes` and `flow3d_addzonetimes`. Verify:

- Outputs are `50, 100, ..., 800` and `10, 50, ..., 800`, respectively.
- Both solve to the boundary event at 10. The saved solution at 10 represents
  the preceding interval; the new boundary applies after 10.
- Final fields, `uniform/time`, observations, and budgets end at 800, with no
  801/804 output and no EOF-related collapse in fixed head.
- A deliberately nonconvergent Picard case exits nonzero without writing its
  rejected target. Check the same behavior for an injected PHREEQC failure.
- Nonreactive restart preserves the loaded time/index and selects the next
  pending output. Reactive restart and `reactionSteps=0` fail explicitly.
- Positive reaction subdivisions reach each event; flow-only runs have no
  reaction-interval timestep cap. Negative reaction cadence is unchanged before
  the final step. Chemistry runs at the final endpoint for any remaining elapsed
  time before observations, budgets, and output; a call already due there runs once.
- With `reactionSteps=-5`, a run ending off cadence adds one final chemistry call
  without adding a transport step. A run ending on cadence does not add a duplicate.
- Set `OMP_NUM_THREADS` before launching an OpenMP build and check the startup
  `PHREEQC thread count` message for the configured and requested counts.

Keep existing result directories and bundled binaries untouched during these
checks. Passing this standalone probe does not establish solver correctness.
