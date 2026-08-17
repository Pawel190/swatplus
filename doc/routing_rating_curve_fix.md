# Routing and rating-curve repair

## Scope

This document describes the hydraulic and routing corrections implemented in:

- `src/sd_rating_curve.f90`
- `src/sd_hydsed_init.f90`

The changes cover channel and flood-plain geometry, rating-curve volumes,
Muskingum parameter initialization, initial Muskingum memory, and numerical
guards. The project was not built or executed as part of this work.

## Parameter conventions

Two different slope conventions are used by the channel data:

- `chss` is the channel side-slope ratio expressed as horizontal change per
  unit vertical change, `z = dx/dy`.
- `fps` is the flood-plain ground slope expressed as vertical change per unit
  horizontal change, `s = dy/dx`.

Consequently, the horizontal expansion on one channel side at water depth `d`
is `chss * d`, whereas the horizontal flood-plain expansion is `d / fps`.

## 1. In-bank channel geometry

Let:

- `b` be the channel-bottom width in metres,
- `d` be the water depth in metres,
- `z = chss` be the channel side-slope ratio.

The bottom width at bankfull geometry is:

```text
b = chw - 2 * chd * chss
```

If this produces a non-positive width, the bottom width is set to half the
bankfull top width and `chss` is recalculated, provided `chd` is positive. The
final width is stored in `ch_rcurv(i)%wid_btm` and copied to
`sd_ch_vel(i)%wid_btm`.

### Previous formulas

```text
wet_perimeter = b + 2 * d * sqrt(1 + 1 / chss^2)
area          = b * d + d / chss
```

These formulas treated `chss` as its inverse. The area expression was also
dimensionally invalid because `b*d` has units of square metres while
`d/chss` has units of metres.

### Corrected formulas

```text
wet_perimeter = b + 2 * d * sqrt(1 + chss^2)
area          = b * d + chss * d^2
top_width     = b + 2 * chss * d
hydraulic_radius = area / wet_perimeter
```

The corrected expressions are used for the rating-curve points at `0.1` and
`1.0` times bankfull depth.

## 2. Overbank and flood-plain geometry

For water depth `d_above` above bankfull level, the additional cross-sectional
area consists of a rectangular layer over the bankfull channel and the
combined triangular flood-plain wedges on both sides:

```text
area_above_bankfull = chw * d_above + d_above^2 / fps
```

The total area and top width are therefore:

```text
area_total = area_bankfull + area_above_bankfull
top_width  = chw + 2 * d_above / fps
```

The wetted perimeter retains the `fps = dy/dx` convention:

```text
wet_perimeter = perimeter_bankfull
              + 2 * d_above * sqrt(1 + 1 / fps^2)
```

### Previous inconsistency

The previous area calculation used the channel-bottom width `b` for the layer
above bankfull. It then calculated flood-plain volume as `top_width*d_above`,
which treats the entire overbank layer as a rectangle at its maximum width.
As a result, rating-curve area and volume described different cross sections.

### Corrected volume and routing-storage partition

The active routing code fills `ch_stor` only to bankfull capacity and assigns
all excess water to `fp_stor`. The rating-curve volumes follow that convention.
With channel length `L` in kilometres:

```text
volume_bankfull   = area_bankfull       * L * 1000
volume_channel    = volume_bankfull
volume_floodplain = area_above_bankfull * L * 1000
volume_total      = volume_channel + volume_floodplain
```

Thus, `volume_floodplain` is a routing-storage category: it contains all water
above bankfull, including the rectangular water layer directly above the
bankfull channel as well as the two flood-plain wedges.

This guarantees:

```text
volume_total = area_total * L * 1000
volume_total = volume_channel + volume_floodplain
```

The partition is subsequently used when initial water and constituent storage
is divided between the channel and flood plain.

## 3. Rating-curve hydraulic guards

The hydraulic radius is calculated only when both area and wetted perimeter
are positive. Otherwise it is set to zero.

Travel time is calculated only for a positive Manning velocity:

```text
travel_time = L / (3.6 * velocity)
```

If velocity is effectively zero, travel time is set to zero to prevent
division by zero or non-finite rating-curve values.

## 4. Muskingum storage-time parameters

Previously, `sd_hydsed_init.f90` duplicated the geometry calculations for
`0.1` and `1.0` bankfull depth. The duplicate contained the same incorrect
`chss` formulas as the old rating curve and could also use a stale local
bottom width after `sd_rating_curve` corrected the geometry.

The duplicate calculations have been removed. Storage times are now derived
from rating-curve travel times. With kinematic-wave celerity equal to `5/3` of
the mean velocity:

```text
storage_time = travel_time / (5/3) = 0.6 * travel_time
```

Therefore:

```text
stor_dis_01bf = 0.6 * rating_curve_point_1%ttime
stor_dis_bf   = 0.6 * rating_curve_point_2%ttime
```

This makes the rating curve and routing coefficients use the same geometry.

## 5. Muskingum weight and coefficient safeguards

The calibration weights `msk_co1` and `msk_co2` are normalized when their sum
is positive. If their sum is effectively zero, defaults of `0.75` and `0.25`
are used locally.

The Muskingum weighting factor is locally bounded to its standard physical
range:

```text
0 <= msk_x <= 0.5
```

Substep calculation divides by the stability limit only when that limit is
positive. Coefficient calculation similarly checks the storage time,
substep count, coefficient denominator, and final coefficient sum.

If a valid set of coefficients cannot be calculated, the safe no-delay
fallback is:

```text
c1 = 1
c2 = 0
c3 = 0
```

This maps current inflow volume directly to outflow volume instead of
producing a division-by-zero or non-finite coefficient.

## 6. Initial Muskingum memory

`sd_ch(i)%in1_vol` and `sd_ch(i)%out1_vol` are Muskingum memory variables and
have units of cubic metres per routing step. They participate directly in the
Muskingum equation implemented in `ch_rtmusk.f90`:

```text
outflow_volume = c1 * current_inflow_volume
               + c2 * previous_inflow_volume
               + c3 * previous_outflow_volume
```

### Previous calculation

```text
initial_volume = flow_rate / seconds_per_step
```

This was dimensionally incorrect:

```text
(m^3/s) / s = m^3/s^2
```

It was also executed before the initial rating curve was interpolated for the
current channel, so the shared `rcurv` value could be unset or belong to a
different channel.

### Corrected calculation

The memory is initialized immediately after `rcurv_interp_dep` for the current
channel:

```text
seconds_per_step = 86400 / nsteps
initial_volume   = flow_rate * seconds_per_step
in1_vol          = initial_volume
out1_vol         = initial_volume
```

The units are now correct:

```text
(m^3/s) * s = m^3
```

For example, a flow of `10 m^3/s` during a one-hour routing step represents:

```text
10 * 3600 = 36000 m^3
```

The correction is operationally important for Muskingum routing
(`bsn_cc%rte == 1`). Variable Storage does not use these memory variables in
its outflow equation, although they are still initialized consistently.

Channels without active storage and configurations with zero routing steps
receive zero initial Muskingum memory.