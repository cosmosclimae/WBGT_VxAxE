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
# F1 — bivariate choropleth (headline level), least-favourable future
# ---------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") |>
  select(iso3 = adm0_a3, geometry)

bivar_pal <- c(
  "low-low"="#e8e8e8","low-mid"="#b8d6be","low-high"="#73ae80",
  "mid-low"="#b5c0da","mid-mid"="#90b2b3","mid-high"="#5a9178",
  "high-low"="#6c83b5","high-mid"="#567994","high-high"="#2a5a5b")

f1df <- metrics |> filter(hazard_ssp=="ssp245", horizon==2085, level==HEADLVL) |>
  mutate(key = paste(A_tercile, V_tercile, sep="-"))
f1 <- world |> left_join(f1df, by="iso3") |>
  ggplot() +
  geom_sf(aes(fill=key), color="grey60", linewidth=0.05) +
  scale_fill_manual(values=bivar_pal, na.value="white", guide="none") +
  coord_sf(crs="+proj=robin") + theme_void() +
  labs(title="Heat-stress exposure x socio-economic vulnerability (2085, SSP2-4.5)",
       subtitle=sprintf("Darkest = both highly exposed (WBGT-aft>=%d d/yr) and highly vulnerable (GVI)", HEADLVL))
leg <- expand_grid(A=c("low","mid","high"), V=c("low","mid","high")) |>
  mutate(key=paste(A,V,sep="-"), A=factor(A,c("low","mid","high")), V=factor(V,c("low","mid","high"))) |>
  ggplot(aes(A,V,fill=key)) + geom_tile() +
  scale_fill_manual(values=bivar_pal, guide="none") +
  labs(x="exposure ->", y="vulnerability ->") +
  theme_minimal(base_size=8) + theme(axis.text=element_blank(), panel.grid=element_blank())
ggsave(file.path(FIGDIR,"F1_bivariate_2085_ssp245.png"),
       f1 + inset_element(leg, 0.02,0.02,0.24,0.34), width=10, height=5.2, dpi=300)

# ---------------------------------------------------------------------------
# F2 — concentration curves at the headline level (annotation moved, legend right)
# ---------------------------------------------------------------------------
curve2 <- curve |> filter(level==HEADLVL) |>
  mutate(grp = paste(ssp_lab[hazard_ssp], horizon))
f2 <- ggplot(curve2, aes(cum_pop, cum_charge, color=grp)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey50") +
  geom_line(linewidth=0.9) +
  scale_x_continuous("Cumulative population share (ranked by GVI, least->most vulnerable)",
                     labels=scales::percent) +
  scale_y_continuous("Cumulative share of exposure burden (person-days)",
                     labels=scales::percent) +
  scale_color_brewer("Scenario / horizon", palette="Set1") +
  labs(caption="Below the dashed 45° line: burden concentrated on the more vulnerable.") +
  theme_minimal(base_size=12) +
  theme(legend.position="right",
        plot.caption=element_text(hjust=0, color="grey30"))
ggsave(file.path(FIGDIR,"F2_concentration.png"), f2, width=7.6, height=5.6, dpi=300)

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
       title="Geographic co-location of climatic and socio-economic relief (2085)") +
  theme_minimal(base_size=12)
ggsave(file.path(FIGDIR,"F3_double_dividend.png"), f3, width=7.8, height=6, dpi=300)

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
  theme_minimal(base_size=12)
ggsave(file.path(FIGDIR,"F4_tau_sensitivity.png"), f4, width=8, height=4.5, dpi=300)

# ---------------------------------------------------------------------------
# F5 — *** the GVI-independent result *** : does CI rise with WBGT severity?
# ---------------------------------------------------------------------------
ci5 <- ci |> mutate(scen = paste(ssp_lab[hazard_ssp], horizon))
f5 <- ggplot(ci5, aes(level, CI, color=scen, group=scen)) +
  geom_line(linewidth=0.9) + geom_point(size=2.5) +
  scale_x_continuous("WBGT severity level (°C)", breaks=sort(unique(ci5$level))) +
  scale_y_continuous("Concentration index (CI) — burden on the vulnerable") +
  scale_color_brewer("Scenario / horizon", palette="Set1") +
  labs(title="Does a more severe heat threshold concentrate the burden on the vulnerable?",
       subtitle="Rising CI with severity = an injustice driven by the hazard itself, independent of GVI assumptions") +
  theme_minimal(base_size=12)
ggsave(file.path(FIGDIR,"F5_CI_by_severity.png"), f5, width=8, height=5.2, dpi=300)

message("Figures (F1-F5) written to ", FIGDIR)
