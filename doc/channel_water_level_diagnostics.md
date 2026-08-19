# Channel water-depth diagnostic

## Purpose

`wat_dep` exposes an existing internal SWAT+ value for comparison with observed river water depth. The change is diagnostic only: it does not change Variable Storage or Muskingum routing, their coefficients, channel storage, inflow, or outflow.

## Output definition

`wat_dep` is appended to the `channel_sdmorph_day`, `channel_sdmorph_mon`, `channel_sdmorph_yr`, and `channel_sdmorph_aa` text/CSV outputs.

| Field | Unit | Definition |
|---|---:|---|
| `wat_dep` | m | Flow depth above the model channel bed, derived from the SWAT+ channel rating curve. Its daily definition depends on `i_fpwet`, as described below. |

For `i_fpwet=0`, `wat_dep` is the time-weighted mean daily depth calculated from `out2%dep` for every routing substep. Dry substeps contribute zero.

For `i_fpwet=1`, the routing path does not provide an equivalent series of routed-outflow substeps. `wat_dep` is therefore the rating-curve depth corresponding to the final mean daily outflow rate, after routing and flow-control calculations are complete:

```text
wat_dep = rating_curve_depth(final_outflow_volume / 86400)
```

Because the rating curve is nonlinear, the depth of mean daily flow is not generally identical to a mean of instantaneous depths. The `i_fpwet=1` value must therefore be interpreted as a daily flow-equivalent depth, not as a reconstructed subdaily mean.

Monthly, yearly, and average-annual values follow the existing `sd_ch_output` aggregation and are arithmetic means of the contributing daily values.

`wat_dep` must not be confused with the existing `depth` column. `depth` is channel bankfull geometry; `wat_dep` is flow-dependent water depth.

## Source-code traceability

1. `src/basin_module.f90`: `rte=0` selects Variable Storage; `rte=1` selects Muskingum.
2. For `i_fpwet=0`, `src/sd_channel_control3.f90` prepares `flo_in`, calls `ch_rtmusk`, and exports `flo_out` and `wat_dep`.
3. `src/ch_rtmusk.f90` implements both routing alternatives used on that path. After routed outflow is calculated, `rcurv_interp_flo` produces `ch_rcurv(jrch)%out2`, including `dep`.
4. For `i_fpwet=1`, `src/sd_channel_control.f90` completes routing and flow controls, then converts only the final mean daily outflow rate to `wat_dep` for output.
5. `src/rcurv_interp_flo.f90` interpolates all rating-curve properties for routing calculations.
6. `src/sd_channel_module.f90` defines `channel_rating_curve_parameters%dep`, output field `sd_ch_output%wat_dep`, and the state-free `rcurv_depth_from_flo` diagnostic function.
7. `src/sd_chanmorph_output.f90` writes the morphology record.
8. `src/header_sd_channel.f90` opens `channel_sdmorph_*` files and writes their headers and units.

For `i_fpwet=0`, no additional call to `rcurv_interp_flo` is made. The already calculated `ch_rcurv(jrch)%out2%dep` is accumulated after each active routing substep. For `i_fpwet=1`, `rcurv_depth_from_flo` reads the established rating curve and returns only its interpolated depth; it does not write `rcurv` or any routing variable.

## Routing equations relevant to the diagnostic

For each routing substep, inflow volume is added to total channel/floodplain storage before outflow is calculated.

Variable Storage:

```text
SC = min(1, scoef_bsn * 2*dt / (2*T_out_previous + dt))
O_t = SC * S_t
```

Muskingum:

```text
O_t = C1*I_t + C2*I_(t-1) + C3*O_(t-1)
```

The Muskingum coefficients are initialized in `src/sd_hydsed_init.f90`. Both methods subsequently calculate:

```text
outflow_rate = O_t / dt
call rcurv_interp_flo(channel, outflow_rate)
wat_dep_day = sum(out2%dep * dt) / sum(dt)
```

Routing therefore determines outflow first; the rating curve maps that outflow to reported depth. `wat_dep` does not feed back into routing.

For `i_fpwet=1`, the diagnostic sequence is:

```text
complete routing, losses, storage and flow controls
final_outflow_rate = final_outflow_volume / 86400
wat_dep = rating_curve_depth(final_outflow_rate)
```

The conversion occurs only when morphology output fields are populated and does not alter the final outflow.

## Comparison with a river gauge

For observations expressed as water depth above a local gauge zero, both series require a common vertical reference. An absolute modeled level can be constructed outside routing as:

```text
H_model = Z_model_bed_at_gauge + wat_dep
```

The modeled channel must correspond spatially to the gauge. For `i_fpwet=0`, daily `wat_dep` is a time-weighted daily mean and should be compared with an observed daily mean calculated over the same calendar day. For `i_fpwet=1`, it is the depth associated with mean daily discharge, so comparisons during strongly varying or event flow must acknowledge that distinction.

At many gauges, discharge is itself inferred from observed stage using a station-specific empirical rating curve. The reverse conversion performed here is methodologically valid, but it uses the SWAT+ model-channel rating curve rather than the gauge rating curve. Agreement in discharge therefore does not guarantee agreement in stage. Applying the station-specific rating curve to modeled discharge outside SWAT+ is an alternative when the objective is to reproduce the gauge-reported stage rather than diagnose the model-channel hydraulics.

## Physical limitations

`wat_dep` is derived from the SWAT+ channel rating curve and routed discharge. It is consistent with the geometry and hydraulic relation used by this implementation, but is not a solution of the full Saint-Venant equations. It does not explicitly represent backwater, downstream control, local structures.

## Non-interference check

Run the same input dataset with baseline and diagnostic builds and compare the pre-existing `flo_in`, `flo_out`, and `flo_stor` fields. They should remain unchanged. This checks that reporting does not interfere with routing; it does not revalidate either algorithm.
