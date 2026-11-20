# Tools (python scripts) to post-process simulation results.

See examples in Python notebooks available in mechanics or electronics directoris

- [ElecExamples](../electronics/ElecExamples.ipynb) : post-processing / plots of some of the examples available in electronics
- [MechanicsExamples](../mechanics/MechanicsExamples.ipynb) : post-processing / plots of some of the examples available in mechanics (plots, comparison with references ...)





# Checks of the results without a reference trajectory

- [ResultChecks.hpp](ResultChecks.hpp) (C++) and [result_checks.py](result_checks.py) (Python):
  for the examples whose trajectory depends on an external library (contact detection by
  Bullet: margins, contact manifold, version) or on round-off, a row-by-row comparison
  with a reference file is not meaningful: a contact detected one step later shifts the
  whole trajectory. These examples check instead
  - features of the trajectory, with explicit tolerances (time of the first impact
    +/- a few steps, restitution, final position),
  - physical properties that need no reference (bounded penetration, energy that does not
    increase, no overlap of bodies, bodies at rest at the end),

  and print a table of the checks; the example fails (exit code 1) if one check fails.
  See mechanics/BulletBouncingBox for an example.

  Keep `siconos::algebra::io::compareRefFile` for the examples computed by the kernel and
  numerics only.
