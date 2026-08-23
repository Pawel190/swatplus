# Soil nutrient and water-quality stability repair

## Commit

This document describes commit
`34571d9d5e8e577700325df5d0a356b44638f632`:

```text
restore soil nutrient initialization and stabilize water-quality calculations
```

The repair addresses invalid floating-point operations and other
compiler-dependent behavior caused by missing or incomplete initialization.
The original failure was reproducible with an IntelLLVM (`ifx`) Windows
Release build configured to trap floating-point exceptions. The traceback
terminated in `NUT_NMINRL` at `src/nut_nminrl.f90:278`, called from
`HRU_CONTROL`.

The traceback identified where the invalid arithmetic was detected, not the
root cause. The underlying soil nutrient and carbon state had not been
initialized for every HRU before `nut_nminrl` used it. Uninitialized memory can
contain different values on Windows and Linux, or in different compiler and
optimization configurations. This explains why the same input could appear to
work on Linux while failing in the Windows `ifx` build.

## Source changes

### 1. Restore soil nutrient and carbon initialization

File: `src/soils_init.f90`

`soil_nutcarb_init(isol)` is called after the soil nutrient, carbon, microbial,
manure, and water arrays have been allocated for each HRU. This restores the
required initialization stage before daily nutrient mineralization begins.

This is the primary correction for the floating-point exception reported in
`nut_nminrl.f90`. It restores an existing model initialization routine; it does
not introduce a new mineralization equation.

### 2. Define algal growth for every water-quality branch

File: `src/ch_watqual4.f90`

The local algal growth rate `gra` is initialized to zero before its conditional
calculation. A `case default` branch also assigns zero when `igropt` does not
select one of the implemented growth equations.

Previously, `gra` could retain an undefined value when `algcon >= 5000` or when
`igropt` was outside the handled cases. The repair preserves both implemented
growth equations and defines the physically neutral result for paths on which
no equation is selected.

### 3. Initialize carbon transformation fluxes

File: `src/carbon_module.f90`

All members of `carbon_soil_transformations` now default to zero. These members
represent daily carbon fluxes and CO2 production. A newly created value must
therefore begin with no transformation until a process calculates a flux.

The change prevents undefined values from entering aggregation and arithmetic
operators. It does not modify the carbon transformation formulas.

### 4. Define the optional plant-region count

File: `src/pl_read_regions_cal.f90`

`nspu` is reset to zero before every region record is read. A record containing
only the required two fields therefore retains its documented meaning: the
region applies to all HRUs.

Without the reset, an omitted optional third field could leave `nspu` with a
value from an earlier iteration or with uninitialized stack data, leading to an
incorrect allocation or region selection.

### 5. Avoid plant lookup for an empty plant community

File: `src/ero_cfactor.f90`

The plant-specific `usle_c` lookup is now executed only when
`pcom(j)%npl > 0`. For an empty plant community, the routine retains the ground
cover factor already calculated from residue and cover instead of dereferencing
an invalid or stale plant index.

The plant-present path and its USLE C-factor equation are unchanged.

## Build configuration used for diagnosis

The Windows Release executable was built with IntelLLVM `ifx`. The relevant
Fortran flags were:

```cmake
set(fdialect "/fpp /free /fpe0 /traceback /fp:precise")
set(fdebug   "/warn:all /check:bounds")
set(frelease "/O2")
```

## Verification

The corrected IntelLLVM Windows Release executable was tested against the
reported SWAT+ input scenario.

- 50 consecutive runs completed successfully (`50/50`).