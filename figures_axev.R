# ============================================================================
# figures_axev.R — reads the CSVs from the Snakemake workflow and renders F1-F5.
# Run on Windows. Edit OUTDIR. Deps: tidyverse sf rnaturalearth patchwork
# ============================================================================
library(tidyverse)
library(sf)
library(rnaturalearth)
library(patchwork)

OUTDIR <- "F:/climae/WBGT-VAE/results"      # <- adjust
FIGDIR <- file.path(OUTDIR, "figures"); dir.create(FIGDIR, showWarnings = FALSE)

metrics <- read_csv(file.path(OUTDIR, "country_metrics.csv"))
curve   <- read_csv(file.path(OUTDIR, "concentration_curve.csv"))
ci      <- read_csv(file.path(OUTDIR, "concentration_index.csv"))
dd      <- read_csv(file.path(OUTDIR, "double_dividend.csv"))
sens    <- read_csv(file.path(OUTDIR, "threshold_sensitivity.csv"))

ssp_lab  <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5")
HEADLVL  <- 32       # headline WBGT level for F1-F4

# ---------------------------------------------------------------------------
# F1 — bivariate choropleth panel: 2 SSP x 2 horizons (headline level)
#      shows the map DARKENING 2055->2085 and LIGHTENING SSP2->SSP1
# ---------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") |>
  select(iso3 = adm0_a3, geometry)

bivar_pal <- c(
  "low-low"="#e8e8e8","low-mid"="#b8d6be","low-high"="#73ae80",
  "mid-low"="#b5c0da","mid-mid"="#90b2b3","mid-high"="#5a9178",
  "high-low"="#6c83b5","high-mid"="#567994","high-high"="#2a5a5b")

mk_bivar <- function(hssp, hz){
  df <- metrics |> filter(hazard_ssp==hssp, horizon==hz, level==HEADLVL) |>
    mutate(key = paste(A_tercile, V_tercile, sep="-"))
  world |> left_join(df, by="iso3") |>
    ggplot() +
    geom_sf(aes(fill=key), color="grey70", linewidth=0.03) +
    scale_fill_manual(values=bivar_pal, na.value="white", guide="none") +
    coord_sf(crs="+proj=robin") + theme_void(base_size=15) +
    labs(title=sprintf("%s — %d", ssp_lab[hssp], hz))
}

panel <- (mk_bivar("ssp126",2055) | mk_bivar("ssp245",2055)) /
         (mk_bivar("ssp126",2085) | mk_bivar("ssp245",2085))

leg <- expand_grid(A=c("low","mid","high"), V=c("low","mid","high")) |>
  mutate(key=paste(A,V,sep="-"), A=factor(A,c("low","mid","high")), V=factor(V,c("low","mid","high"))) |>
  ggplot(aes(A,V,fill=key)) + geom_tile() +
  scale_fill_manual(values=bivar_pal, guide="none") +
  labs(x="exposure ->", y="vulnerability ->") +
  theme_minimal(base_size=15) +
  theme(axis.text=element_blank(), panel.grid=element_blank(),
        axis.title=element_text(size=12))

f1 <- panel +
  plot_annotation(
    title = sprintf("Heat-stress exposure x socio-economic vulnerability (WBGT-aft>=%d d/yr)", HEADLVL),
    subtitle = "Darkest = highly exposed AND highly vulnerable. Compare across rows (horizon) and columns (scenario).",
    theme = theme(plot.title = element_text(size=22, face="bold"),
                  plot.subtitle = element_text(size=15))) +
  inset_element(leg, left = 0, bottom = 0.005, right = 0.16, top = 0.27,
                align_to = "full")
# bigger canvas so the four maps breathe and the legend clears South America
ggsave(file.path(FIGDIR,"F1_bivariate_panel.png"), f1, width=15, height=9, dpi=300)

# ---------------------------------------------------------------------------
# F2 — concentration curves, FACETED BY WBGT LEVEL (28/30/32)
#      the curve peels away from the 45° line as the threshold gets more severe
# ---------------------------------------------------------------------------
curve2 <- curve |> mutate(grp = paste(ssp_lab[hazard_ssp], horizon),
                          lvl = factor(paste0("WBGT >= ", level, " C"),
                                       levels=paste0("WBGT >= ", sort(unique(level)), " C")))
f2 <- ggplot(curve2, aes(cum_pop, cum_charge, color=grp)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey50") +
  geom_line(linewidth=0.8) +
  facet_wrap(~lvl) +
  scale_x_continuous("Cumulative population share (ranked by GVI, least->most vulnerable)",
                     labels=scales::percent, breaks=c(0,.5,1)) +
  scale_y_continuous("Cumulative share of exposure burden (person-days)",
                     labels=scales::percent) +
  scale_color_brewer("Scenario / horizon", palette="Set1") +
  labs(caption="Below the dashed 45° line: burden concentrated on the more vulnerable. The gap widens with severity.") +
  theme_minimal(base_size=14) +
   theme(legend.position="bottom",
        plot.caption=element_text(hjust=0, color="grey30"),
        axis.title.y=element_text(size=12))
ggsave(file.path(FIGDIR,"F2_concentration_by_level.png"), f2, width=10, height=4.6, dpi=300)

# ---------------------------------------------------------------------------
# F3 — double dividend (headline level, 2085) : co-LOCATION of the two reliefs
# ---------------------------------------------------------------------------
f3 <- ggplot(dd, aes(dA_bar, dpgvi, size=pop)) +
  geom_hline(yintercept=0, color="grey70") + geom_vline(xintercept=0, color="grey70") +
  geom_point(alpha=0.55, color="#2a5a5b") +
  scale_size_area("Population", max_size=12, labels=scales::comma) +
  annotate("label", x=Inf, y=Inf, hjust=1.02, vjust=1.3, label.size=0,
           fill="white", alpha=0.7, size=3.2,
           label="upper-right: SSP1 relieves BOTH\nhazard and vulnerability in the same country") +
  labs(x="Hazard relief under SSP1  (A_bar: SSP2 - SSP1, days/yr)",
       y="Vulnerability relief under SSP1  (GVI: SSP2 - SSP1)",
       title="Co-location of climatic & socio-economic relief (2085)") +
  theme_minimal(base_size=15) +
  theme(plot.title=element_text(size=16))
ggsave(file.path(FIGDIR,"F3_double_dividend.png"), f3, width=10, height=6.5, dpi=300)

# ---------------------------------------------------------------------------
# F4 — tau-persistence sensitivity at the headline level
# ---------------------------------------------------------------------------
f4 <- sens |> filter(level==HEADLVL) |>
  mutate(scen=ssp_lab[hazard_ssp], tau=factor(tau)) |>
  ggplot(aes(tau, global_N/1e6, fill=scen)) +
  geom_col(position="dodge") + facet_wrap(~horizon) +
  scale_fill_brewer("Scenario", palette="Set1") +
  labs(x=sprintf("Persistence threshold tau (days/yr, WBGT-aft>=%d)", HEADLVL),
       y="Global population in unsafe regime (millions)",
       title="Robustness to the persistence threshold") +
  theme_minimal(base_size=15)
ggsave(file.path(FIGDIR,"F4_tau_sensitivity.png"), f4, width=9, height=4.8, dpi=300)

# ---------------------------------------------------------------------------
# F5 — *** the GVI-independent result *** : does CI rise with WBGT severity?
# ---------------------------------------------------------------------------
ci5 <- ci |> mutate(scen = paste(ssp_lab[hazard_ssp], horizon))
f5 <- ggplot(ci5, aes(level, CI, color=scen, group=scen)) +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  scale_x_continuous("WBGT severity level (°C)", breaks=sort(unique(ci5$level))) +
  scale_y_continuous("Concentration index (CI) — burden on the vulnerable") +
  scale_color_brewer("Scenario / horizon", palette="Set1") +
  labs(title="Severity concentrates the burden on the vulnerable",
       subtitle=str_wrap("Rising CI with severity = injustice driven by the hazard itself, independent of GVI assumptions (note SSP1-2.6 2085, the exception).", 70)) +
  theme_minimal(base_size=15) +
  theme(plot.title=element_text(size=17), plot.subtitle=element_text(size=12))
ggsave(file.path(FIGDIR,"F5_CI_by_severity.png"), f5, width=10, height=6, dpi=300)

message("Figures (F1-F5) written to ", FIGDIR)
