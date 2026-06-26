# Sport-relevant heat exposure × socio-economic vulnerability (A × E × V)

A reproducible workflow quantifying, at country level, how the burden of
sport-relevant heat exposure is distributed across socio-economic vulnerability,
under coherent Shared Socioeconomic Pathways.

It couples three layers within an explicit **hazard × exposure × vulnerability** frame:

| Component | Definition | Source |
|-----------|------------|--------|
| **A — Hazard** | Annual-mean days with afternoon WBGT ≥ {28, 30, 32} °C (30-yr climatology) | Bias-corrected CMIP6 (NEX-GDDP), 5-model ensemble mean |
| **E — Exposure** | Gridded population counts (1 km), aggregated to 0.25° by sum | WorldPop SSP projections (Bondarenko et al. 2025) |
| **V — Vulnerability** | National Global Vulnerability Index (`pgvi`, 0–100) | GVI Projections Database (Huisman et al. 2025) |

Two internally coherent futures (SSP1-2.6 and SSP2-4.5, each paired with its own
socio-economic pathway) are assessed at mid- (2055) and late-century (2085) horizons.

---

## What the workflow computes

The single grid-cell operation is `b_i = P_i × A_i` (**person-days**). Population is
summed; the hazard is never summed; the GVI is joined at **country level only**
(never projected onto the grid). From there:

- **Population-weighted exposure** per country (`a_bar_c`, `charge_c`, headcounts `N(τ)`).
- **Concentration index (CI)** of the burden ranked by GVI — pop-weighted,
  CI > 0 ⇔ burden concentrated on the more vulnerable.
- **Severity sweep**: CI computed at WBGT ≥ 28, 30, 32 °C — tests whether a more
  severe heat threshold concentrates the burden further on the vulnerable
  (a signal *independent* of the GVI's built-in pathway assumptions).
- **Double dividend**: per-country, the geographic co-location of climatic relief
  (ΔA) and vulnerability relief (ΔGVI) under the sustainable pathway.

---

## Repository structure

```
.
├── Snakefile             # workflow DAG
├── config.yaml           # all paths, scenarios, horizons, thresholds
├── axev_lib.py           # processing functions (hazard, pop aggregation, zonal, CI)
├── figures_axev.R        # F1–F5 figures (reads the result CSVs)
└── results/              # generated
    ├── country_metrics.csv          # 1 row per (iso3, ssp, horizon, level)
    ├── concentration_index.csv      # CI per (ssp, horizon, level)   ← key table
    ├── concentration_curve.csv      # concentration-curve points
    ├── threshold_sensitivity.csv    # global headcount per (ssp, horizon, level, τ)
    └── double_dividend.csv          # per-country ΔA vs ΔGVI (headline level, 2085)
```

---

## Requirements

Python (processing) and R (figures), kept separate by design.

**Python** (conda-forge recommended for the geospatial stack):
```bash
conda create -n axev -c conda-forge \
  snakemake-minimal numpy pandas xarray rioxarray rasterio geopandas netcdf4
conda activate axev
```

**R** (figures): `tidyverse`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `patchwork`.

---

## Usage

1. Edit **`config.yaml`** — point the paths to your local data and verify the
   hazard variable names (`cdo showname <file>`). Hazard thresholds are defined as
   **lists of day-count variables to sum** (binned classes → cumulative levels).

2. Run the workflow:
   ```bash
   snakemake -n                              # dry-run: inspect the DAG
   snakemake --cores 4 --resources mem_mb=9000
   ```

3. Render the figures (R):
   ```r
   source("figures_axev.R")     # writes F1–F5 to results/figures/
   ```

Visualise the DAG: `snakemake --dag | dot -Tpng > dag.png`.

---

## Figures

- **F1** — Bivariate map: exposure × vulnerability (least-favourable future).
- **F2** — Concentration curves (burden vs population ranked by GVI).
- **F3** — Double dividend: co-location of climatic and socio-economic relief.
- **F4** — Robustness to the persistence threshold τ.
- **F5** — CI vs WBGT severity (the GVI-independent result).

---

## Methodological notes

- **Hazard = 30-yr climatological mean** centred on the year at which population
  and vulnerability are sampled (2055, 2085). Population and GVI are point snapshots.
- **WBGT** is the shaded afternoon index built on daily Tmax + reconstructed RH,
  i.e. the thermodynamic daily peak — longitude-agnostic, no solar-time reconstruction.
- **Resolution discipline**: population (1 km) is aggregated *up* to the 0.25° hazard
  grid by sum; the GVI (national) is joined at country level. Nothing finer-than-source
  is ever invented.
- **No SSP5**: the GVI Projections Database does not provide SSP5; the analysis is
  restricted to plausible low/intermediate pathways (SSP1-2.6, SSP2-4.5).
- **Screening-level estimates.** Absolute exceedance counts depend on a simplified
  WBGT formulation and gridded inputs; relative and distributional results
  (CI, concentration) are robust to systematic biases that cancel under ranking.
- **CI interpretation**: positive ⇒ burden on the vulnerable. When a concentration
  curve crosses the diagonal, the single CI value nets opposing segments and should
  be read alongside the curve (F2), not alone.

---

## Data sources

- **CMIP6 / NEX-GDDP-CMIP6** — NASA Earth Exchange (Thrasher et al. 2022).
- **WorldPop SSP population** — Bondarenko et al. 2025, doi:10.5258/SOTON/WP00849.
- **GVI Projections Database** — Huisman et al. 2025 (Global Data Lab), *Scientific Data*.
- **Country borders** — Natural Earth (10 m admin-0).

## Citation

If you use this workflow, please cite the associated manuscript (in preparation,
CosmosClimae) and the upstream data sources above.

## License

Code released under the MIT License. Input datasets retain their original licenses.
