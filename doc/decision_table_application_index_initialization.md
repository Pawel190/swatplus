# Land-use decision-table application-index initialization

## Commit

This document describes commit
`56320338eb9bddef7973594f3eacc0119dd580ea`:

```text
initialize land-use decision-table application indices
```

## Failure

An IntelLLVM (`ifx`) Windows Release executable could fail nondeterministically
with an access violation in `MGT_TRANSPLANT` at `src/mgt_transplant.f90:22`,
called from the land-use decision-table action handler in `src/actions.f90`.
The same executable and inputs could also complete successfully because the
failure depended on undefined memory contents.

The failure is unrelated to GWFLOW. It was exposed while running a model that
also enabled the groundwater-flow module.

## Root cause

`dtbl_lum(i)%act_app` was allocated in `src/dtbl_lum_read.f90` without an
initial value. Its elements were assigned only when an action's application
name matched a corresponding database entry.

When `transplant.plt` was absent and a plant action used a null transplant
pointer, no assignment occurred. The plant action later tested the undefined
integer as an application index. A coincidental positive value caused an
out-of-range access to the transplant database and an access violation.

## Correction

The complete `act_app` array is initialized to zero immediately after
allocation:

```fortran
allocate (dtbl_lum(i)%act_app(dtbl_lum(i)%acts))
dtbl_lum(i)%act_app = 0
```

Zero retains the existing meaning of "no application record". Successful
crosswalks continue to replace zero with the valid one-based database index.

No guard was added to `actions.f90`, and no transplant or GWFLOW calculation
was changed.
