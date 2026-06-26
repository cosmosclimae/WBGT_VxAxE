# ============================================================================
# Snakefile — A x E x V workflow with WBGT-level sweep (28/30/32)
#   Run:  snakemake --cores 4 --resources mem_mb=9000
#   (mem_mb caps parallel pop jobs; with windowed pop you can also use --cores 4)
# ============================================================================
import os, sys, json
import numpy as np
import pandas as pd

configfile: "config.yaml"
sys.path.insert(0, workflow.basedir)
import axev_lib as lib

OUT      = config["outdir"]
RES      = os.path.join(OUT, "results")
HSSPS    = list(config["scenarios"].keys())
HORIZONS = [str(h) for h in config["horizons"]]
LEVELS   = [str(l) for l in config["hazard_levels"].keys()]     # 28 / 30 / 32
HEADLVL  = str(config["headline_level"])
os.makedirs(RES, exist_ok=True)

REF_GRID = f"{RES}/A_{HSSPS[0]}_{HORIZONS[0]}_{HEADLVL}.nc"      # shared 0.25 grid

wildcard_constraints:
    hssp = r"ssp\d+",
    horizon = r"\d{4}",
    lvl = r"\d+",


rule all:
    input:
        f"{RES}/country_metrics.csv",
        f"{RES}/concentration_curve.csv",
        f"{RES}/concentration_index.csv",       # <- CI per ssp x horizon x level (key new table)
        f"{RES}/threshold_sensitivity.csv",
        f"{RES}/double_dividend.csv",


# ---------------------------------------------- hazard -> A grid, per level
rule hazard:
    input:  lambda wc: config["wbgt_paths"][f"{wc.hssp}__{wc.horizon}"]
    output: f"{RES}/A_{{hssp}}_{{horizon}}_{{lvl}}.nc"
    run:
        varnames = config["hazard_levels"][int(wildcards.lvl)]
        A = lib.load_hazard(input[0], varnames, config["n_years"], config["wbgt_lon_0_360"])
        lib.save_grid(A, output[0])


# ------------------------------- country id raster (built once, shared grid)
rule country_mask:
    input:  REF_GRID
    output: cid=f"{RES}/country_id.npy", ids=f"{RES}/id2iso.json"
    run:
        A = lib.load_grid(input[0])
        cid, id2iso = lib.build_country_id(A, config["borders_shp"], config["borders_iso_field"])
        np.save(output.cid, cid)
        json.dump({int(k): v for k, v in id2iso.items()}, open(output.ids, "w"))


# ----------------------------- population on grid (per ssp x horizon, NOT level)
rule pop_on_grid:
    input:  grid=REF_GRID                       # any level shares the same grid
    output: f"{RES}/P_{{hssp}}_{{horizon}}.npy"
    resources:
        mem_mb=9000
    run:
        pop_ssp = config["scenarios"][wildcards.hssp][0]
        pop_path = config["pop_template"].format(SSP=pop_ssp, YEAR=wildcards.horizon)
        A = lib.load_grid(input.grid)
        P = lib.load_pop_on_grid(pop_path, A)
        np.save(output[0], P.values)
        print(f"[{wildcards.hssp} {wildcards.horizon}] global pop = {P.values.sum()/1e9:.2f} bn")


# ----------------------------------------- zonal metrics, per ssp x horizon x level
rule zonal:
    input:
        A=f"{RES}/A_{{hssp}}_{{horizon}}_{{lvl}}.nc",
        P=f"{RES}/P_{{hssp}}_{{horizon}}.npy",
        cid=f"{RES}/country_id.npy",
        ids=f"{RES}/id2iso.json",
    output:
        metrics=f"{RES}/metrics_{{hssp}}_{{horizon}}_{{lvl}}.csv",
        sens=f"{RES}/sens_{{hssp}}_{{horizon}}_{{lvl}}.csv",
    run:
        A = lib.load_grid(input.A).values
        P = np.load(input.P)
        cid = np.load(input.cid)
        id2iso = {int(k): v for k, v in json.load(open(input.ids)).items()}
        assert A.shape == P.shape == cid.shape, f"grid mismatch {A.shape} {P.shape} {cid.shape}"

        df, sens = lib.zonal(A, P, cid, id2iso, config["taus"])

        gvi = lib.load_gvi(config["gvi_csv"], config["gvi_iso"], config["gvi_year"],
                           config["gvi_ssp"], config["gvi_val"])
        gssp = config["scenarios"][wildcards.hssp][1]
        gsub = gvi[(gvi.ssp_gvi == gssp) & (gvi.year == int(wildcards.horizon))][["iso3", "pgvi"]]
        df = df.merge(gsub, on="iso3", how="left")

        df["hazard_ssp"] = wildcards.hssp
        df["pop_ssp"] = config["scenarios"][wildcards.hssp][0]
        df["horizon"] = int(wildcards.horizon)
        df["level"] = int(wildcards.lvl)
        df.to_csv(output.metrics, index=False)

        s = pd.DataFrame(sens)
        s["hazard_ssp"] = wildcards.hssp; s["horizon"] = int(wildcards.horizon)
        s["level"] = int(wildcards.lvl)
        s.to_csv(output.sens, index=False)


# --------------------------------------------------- merge metrics + terciles
rule merge_metrics:
    input:
        expand(f"{RES}/metrics_{{hssp}}_{{horizon}}_{{lvl}}.csv",
               hssp=HSSPS, horizon=HORIZONS, lvl=LEVELS)
    output: f"{RES}/country_metrics.csv"
    run:
        m = pd.concat([pd.read_csv(f) for f in input], ignore_index=True)
        grp = ["hazard_ssp", "horizon", "level"]
        m["A_tercile"] = m.groupby(grp)["a_bar_c"].transform(lib.tercile)
        m["V_tercile"] = m.groupby(grp)["pgvi"].transform(lib.tercile)
        cols = ["iso3", "hazard_ssp", "pop_ssp", "horizon", "level", "a_bar_c", "charge_c",
                "pop_c", "pgvi", "A_tercile", "V_tercile"] + [f"N_{t}" for t in config["taus"]]
        m[cols].to_csv(output[0], index=False)


# ------------------------------------ concentration curve + index (per level)
rule concentration:
    input:  f"{RES}/country_metrics.csv"
    output:
        curve=f"{RES}/concentration_curve.csv",
        index=f"{RES}/concentration_index.csv",
    run:
        m = pd.read_csv(input[0])
        curves, ci = [], []
        for (hssp, hz, lvl), g in m.groupby(["hazard_ssp", "horizon", "level"]):
            CI, curve = lib.concentration(g)
            curve["hazard_ssp"] = hssp; curve["horizon"] = hz; curve["level"] = lvl
            curves.append(curve)
            ci.append({"hazard_ssp": hssp, "horizon": hz, "level": lvl, "CI": CI})
        pd.concat(curves, ignore_index=True).to_csv(output.curve, index=False)
        cidf = pd.DataFrame(ci).sort_values(["hazard_ssp", "horizon", "level"])
        cidf.to_csv(output.index, index=False)
        print("\n=== Concentration index (CI) — does severity raise the burden on the vulnerable? ===")
        print(cidf.to_string(index=False))


# --------------------------------------------------- threshold sensitivity
rule sensitivity:
    input:
        expand(f"{RES}/sens_{{hssp}}_{{horizon}}_{{lvl}}.csv",
               hssp=HSSPS, horizon=HORIZONS, lvl=LEVELS)
    output: f"{RES}/threshold_sensitivity.csv"
    run:
        s = pd.concat([pd.read_csv(f) for f in input], ignore_index=True)
        s.to_csv(output[0], index=False)


# ------------------------------------------- double dividend (headline level)
rule double_dividend:
    input:  f"{RES}/country_metrics.csv"
    output: f"{RES}/double_dividend.csv"
    run:
        m = pd.read_csv(input[0])
        hz = config["dd_horizon"]; lvl = config["headline_level"]
        sub = m[(m.horizon == hz) & (m.level == lvl)]
        s1 = sub[sub.hazard_ssp == "ssp126"]; s2 = sub[sub.hazard_ssp == "ssp245"]
        dd = (s1[["iso3", "pop_c", "a_bar_c", "pgvi"]]
              .merge(s2[["iso3", "a_bar_c", "pgvi"]], on="iso3", suffixes=("_ssp1", "_ssp2")))
        dd["dA_bar"] = dd.a_bar_c_ssp2 - dd.a_bar_c_ssp1
        dd["dpgvi"] = dd.pgvi_ssp2 - dd.pgvi_ssp1
        dd["pop"] = dd.pop_c
        dd.to_csv(output[0], index=False)
