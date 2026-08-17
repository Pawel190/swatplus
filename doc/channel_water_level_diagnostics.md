# Channel water-depth diagnostic

## Purpose

`wat_dep` exposes an existing internal SWAT+ value for comparison with observed river water depth. The change is diagnostic only: it does not change Variable Storage or Muskingum routing, their coefficients, channel storage, inflow, or outflow.

## Output definition

`wat_dep` is appended to the `channel_sdmorph_day`, `channel_sdmorph_mon`, `channel_sdmorph_yr`, and `channel_sdmorph_aa` text/CSV outputs.

| Field | Unit | Definition |
|---|---:|---|
| `wat_dep` | m | Time-weighted mean daily flow depth above the model channel bed, calculated from `out2%dep` for every routing substep. Dry substeps contribute zero. |

Monthly, yearly, and average-annual values follow the existing `sd_ch_output` aggregation and are arithmetic means of the contributing daily values.

`wat_dep` must not be confused with the existing `depth` column. `depth` is channel bankfull geometry; `wat_dep` is flow-dependent water depth.

## Source-code traceability

1. `src/basin_module.f90`: `rte=0` selects Variable Storage; `rte=1` selects Muskingum.
2. `src/sd_channel_control3.f90`: prepares `flo_in`, calls `ch_rtmusk`, and exports `flo_out` and `wat_dep`.
3. `src/ch_rtmusk.f90`: implements both routing alternatives. After routed outflow is calculated, `rcurv_interp_flo` produces `ch_rcurv(jrch)%out2`, including `dep`.
4. `src/rcurv_interp_flo.f90`: interpolates rating-curve properties for the routed outflow rate.
5. `src/sd_channel_module.f90`: defines `channel_rating_curve_parameters%dep` and output field `sd_ch_output%wat_dep`.
6. `src/sd_chanmorph_output.f90`: writes the morphology record.
7. `src/header_sd_channel.f90`: opens `channel_sdmorph_*` files and writes their headers and units.

No additional call to `rcurv_interp_flo` is made. The already calculated `ch_rcurv(jrch)%out2%dep` is accumulated after each active routing substep, preventing changes to shared rating-curve state.

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

## Comparison with a river gauge

For observations expressed as water depth above a local gauge zero, both series require a common vertical reference. An absolute modeled level can be constructed outside routing as:

```text
H_model = Z_model_bed_at_gauge + wat_dep
```

The modeled channel must correspond spatially to the gauge. Daily `wat_dep` is a time-weighted daily mean and should be compared with an observed daily mean calculated over the same civil-day interval and time zone.

## Physical limitations

`wat_dep` is derived from the SWAT+ channel rating curve and routed discharge. It is consistent with the geometry and hydraulic relation used by this implementation, but is not a solution of the full Saint-Venant equations. It does not explicitly represent backwater, downstream control, local structures, or a surveyed gauge datum.

It is therefore suitable as a calibration diagnostic where rating-curve assumptions are appropriate. Calibration against absolute gauge elevation additionally requires a defensible channel-bed elevation and common datum.

## Non-interference check

Run the same input dataset with baseline and diagnostic builds and compare the pre-existing `flo_in`, `flo_out`, and `flo_stor` fields. They should remain unchanged. This checks that reporting does not interfere with routing; it does not revalidate either algorithm.

Record the SWAT+ release and Git commit in publications. The public output reference is <https://swatplus.com/t3m4Docs/output-reference/channel-sediment-nutrient/>.
