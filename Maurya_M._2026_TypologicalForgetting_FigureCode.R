# =========================================================================
#  Maurya_2026_TypologicalForgetting_FigureCode.R
#
#  Figure-generation code for the figures in:
#
#    Maurya, M. (2026). Beyond the single forgetting curve: A typologically-
#    extended Wickelgren function, a simulation study, and a pre-registered
#    test of error-type-specific memory decay. PsyArXiv preprint.
#
#  Figure 1.  Predicted divergence of error-type forgetting curves under
#             the typologically-extended Wickelgren model at prior-mean
#             parameters, with a shaded reversion gap making the
#             Confabulation reversion visible.
#  Figure 2.  Operating-characteristics surface from the 7,200-fit
#             simulation in run_simulation.py / SimCode.R.
#
#  AUTHOR:    Manik Maurya (ORCID 0009-0005-3554-693X)
#  CONTACT:   manikmaurya.in@gmail.com
#  LICENCE:   CC BY 4.0
#  REQUIRES:  R >= 4.2,  ggplot2 >= 3.4,  dplyr, tidyr, jsonlite,
#             ggtext, scales, patchwork, viridis
# =========================================================================

suppressPackageStartupMessages({
  library(ggplot2);  library(dplyr);    library(tidyr)
  library(jsonlite); library(ggtext);   library(scales)
  library(patchwork);library(viridis)
})

# ----- 0. PRIOR-MEAN PARAMETERS (matches Table 3 of the preprint) ----------

TRUTH <- list(
  lambda = c(RF = 0.50, PK = 0.85, CF = 1.30, IN = 0.80),
  psi    = c(RF = 1.40, PK = 1.00, CF = 0.75, IN = 1.00),
  sigma  = c(RF = 0.80, PK = 0.40, CF = 0.20, IN = 0.40),
  gamma  = 1.30, tau_CF = 7.0,
  A_IN   = 0.55, kappa  = 10.0,
  omega  = 2*pi/8.0, phi = 0.0
)

inv_logit <- function(x) 1 / (1 + exp(-x))

# Compute R_e(t) on the probability scale; no relearning savings (pure
# forgetting trajectory). For Confabulation, the with-and-without-reversion
# curves are returned so the reversion gap can be shaded.

R_e <- function(t, etype, truth = TRUTH, no_reversion = FALSE) {
  eta <- truth$lambda[etype] - truth$psi[etype] * log1p(t)
  if (etype == "CF" && !no_reversion) {
    eta <- eta - truth$gamma * (1 - exp(-t / truth$tau_CF))
  } else if (etype == "IN") {
    eta <- eta + truth$A_IN * exp(-t / truth$kappa) *
                  cos(truth$omega * t + truth$phi)
  }
  inv_logit(eta)
}

# =========================================================================
#  FIGURE 1 — predicted divergence of the four error-type curves
# =========================================================================

t_grid <- seq(0, 28, length.out = 1121)
fig1_df <- bind_rows(
  lapply(c("RF","PK","CF","IN"), function(et) {
    tibble(t = t_grid, type = et, R = R_e(t_grid, et))
  }),
  tibble(t = t_grid, type = "CF_norev",
         R = R_e(t_grid, "CF", no_reversion = TRUE))
)

label_map <- c(
  PK        = "Partial Knowledge (PK) — reference (classical Wickelgren)",
  RF        = "Recall Failure (RF) — fastest decay; largest savings",
  IN        = "Interference (IN) — damped oscillation (exploratory)",
  CF        = "Confabulation (CF) — high λ<sub>CF</sub>, slow ψ<sub>CF</sub>, late reversion (γ)",
  CF_norev  = "CF without reversion term (γ = 0) — reference"
)
fig1_df$lab <- factor(fig1_df$type, levels = names(label_map),
                      labels = label_map[names(label_map)])

cols <- c("PK" = "#2c3e50", "RF" = "#c0392b",
          "IN" = "#117a65", "CF" = "#7d3c98", "CF_norev" = "#c39bd3")

# Shade the gap between CF-with-reversion and CF-without-reversion
shade_df <- tibble(
  t        = t_grid,
  ymin_val = R_e(t_grid, "CF"),
  ymax_val = R_e(t_grid, "CF", no_reversion = TRUE)
)

fig1 <- ggplot(fig1_df, aes(x = t, y = R, group = type)) +
  geom_ribbon(data = shade_df,
              aes(x = t, ymin = ymin_val, ymax = ymax_val),
              inherit.aes = FALSE,
              fill = "#7d3c98", alpha = 0.10) +
  geom_line(data = filter(fig1_df, type %in% c("PK","RF","IN")),
            aes(colour = type), linewidth = 1.0) +
  geom_line(data = filter(fig1_df, type == "CF_norev"),
            aes(colour = type), linewidth = 0.7, linetype = "dashed",
            alpha = 0.65) +
  geom_line(data = filter(fig1_df, type == "CF"),
            aes(colour = type), linewidth = 1.15) +
  geom_vline(xintercept = c(0, 2, 7, 14),
             colour = "#888", linewidth = 0.25, linetype = "dotted") +
  annotate("text", x = c(0, 2, 7, 14), y = 1.03,
           label = c("0", "2d", "7d", "14d"),
           size = 2.6, colour = "#444") +
  annotate("text", x = 7.5, y = 1.08, hjust = 0.5,
           label = "registered retention assessments",
           size = 2.6, colour = "#444", fontface = "italic") +
  annotate("richtext", x = 3.6, y = 0.88, hjust = 0,
           label = "Immediate hypercorrection: **λ<sub>CF</sub> &gt; λ<sub>PK</sub>**<br><span style='color:#5b2a86;font-size:7pt'>(Butterfield &amp; Metcalfe, 2001)</span>",
           label.colour = NA, fill = NA, size = 2.7, colour = "#5b2a86") +
  annotate("richtext", x = 15.7, y = 0.56, hjust = 0,
           label = paste0("<span style='color:#5b2a86'>Reversion gap (γ-driven):</span><br>",
                          "the one-week return of<br>",
                          "high-confidence errors<br>",
                          "<span style='font-size:7pt'>(Butler, Fazio, &amp; Marsh, 2011;<br>",
                          "blocked by intervening test —<br>",
                          "Metcalfe &amp; Miele, 2014)</span>"),
           label.colour = NA, fill = NA, size = 2.7, colour = "#5b2a86") +
  scale_colour_manual(values = cols, breaks = c("PK","RF","IN","CF_norev","CF"),
                      labels = label_map[c("PK","RF","IN","CF_norev","CF")],
                      name = NULL) +
  scale_x_continuous(breaks = c(0,2,7,14,21,28),
                     limits = c(-0.5, 28.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.10)) +
  labs(x = "Retention interval t (days)",
       y = expression(paste("Predicted retention probability ", R[e], "(t)"))) +
  theme_minimal(base_size = 9.5, base_family = "serif") +
  theme(
    legend.text       = element_markdown(size = 7.5),
    legend.position   = c(0.30, 0.10),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line         = element_line(linewidth = 0.3)
  )

ggsave("figures/Figure1_divergence.pdf", fig1,
       width = 7.6, height = 4.7, device = cairo_pdf)
ggsave("figures/Figure1_divergence.png", fig1,
       width = 7.6, height = 4.7, dpi = 300)
message("Figure 1 written.")

# =========================================================================
#  FIGURE 2 — operating-characteristics surface from the 7,200-fit run
# =========================================================================
#
#  Reads the JSON output written by run_simulation.py (or the
#  brms-equivalent grid_sensitivity() in SimCode.R written to .rds).
#  If a brms .rds is preferred, swap fromJSON() for readRDS() below.

read_power_grid <- function(path = "figures/power_grid.json") {
  if (!file.exists(path))
    stop("power_grid.json not found. Run run_simulation.py first.")
  fromJSON(path) |> as_tibble()
}
pw <- read_power_grid()

# ----- Panel A: heatmaps for each pre-registered contrast (4-point design)
pw4 <- pw |> filter(design == "4pt")
panel_meta <- tribble(
  ~metric,         ~title,                                ~status,
  "power_P1a",     "P1a:  ψ_RF > ψ_PK",                   "well-powered",
  "power_P1b",     "P1b:  ψ_PK > ψ_CF",                   "low power",
  "power_P2",      "P2:   γ > 0  (CF reversion)",         "low power",
  "power_P3",      "P3:   σ_RF > σ_PK",                   "well-powered",
  "M3_pref_rate",  "M3 over M0  (model preference)",      "well-powered"
)

heat_panel <- function(m, title, status) {
  ggplot(pw4, aes(x = factor(n), y = factor(cf_rate),
                  fill = .data[[m]])) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", .data[[m]])),
              size = 2.2,
              colour = ifelse(pw4[[m]] > 0.55, "white", "#222")) +
    scale_fill_gradientn(
      colours = c("#f4f0fa","#e6c8e8","#d18ccf","#9a4ba8","#3b1144"),
      limits  = c(0, 1), name = "power") +
    labs(x = "Sample size n", y = "CF base rate", title = title,
         subtitle = sprintf("[%s]", status)) +
    theme_minimal(base_size = 8, base_family = "serif") +
    theme(legend.position = "none",
          plot.title    = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(
                            size = 7,
                            colour = ifelse(status == "well-powered",
                                            "#27ae60", "#c0392b"),
                            face = "bold"),
          panel.grid    = element_blank())
}
panels <- mapply(heat_panel, panel_meta$metric,
                 panel_meta$title, panel_meta$status,
                 SIMPLIFY = FALSE)
top_row <- wrap_plots(panels, ncol = 5)

# ----- Panel B: design comparison at n=90, cf=0.15
pw_compare <- pw |> filter(n == 90, cf_rate == 0.15) |>
  pivot_longer(c(power_P1a, power_P1b, power_P2, power_P3, M3_pref_rate),
               names_to = "metric", values_to = "power") |>
  mutate(metric = recode(metric,
                         "power_P1a" = "P1a", "power_P1b" = "P1b",
                         "power_P2"  = "P2",  "power_P3"  = "P3",
                         "M3_pref_rate" = "M3 vs M0"),
         metric = factor(metric, levels = c("P1a","P1b","P2","P3","M3 vs M0")),
         design = factor(design, levels = c("3pt","4pt","5pt")))

bottom_panel <- ggplot(pw_compare,
                       aes(x = metric, y = power, fill = design)) +
  geom_col(position = position_dodge(0.78), width = 0.7,
           colour = "#333", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", power)),
            position = position_dodge(0.78), vjust = -0.4, size = 2.4) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             colour = "#27ae60", linewidth = 0.4) +
  annotate("text", x = 0.6, y = 0.83, label = "conventional 0.80",
           colour = "#27ae60", size = 2.6, hjust = 0) +
  scale_fill_manual(values = c("3pt"="#bdc3c7","4pt"="#7d3c98","5pt"="#27ae60"),
                    name = "Design") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Recovery rate",
       title = "Panel B. Retention-design comparison at registered n = 90, CF base rate = .15") +
  theme_minimal(base_size = 9, base_family = "serif") +
  theme(plot.title    = element_text(size = 9, face = "bold"),
        legend.position = "top",
        panel.grid.major.x = element_blank())

fig2 <- top_row / bottom_panel + plot_layout(heights = c(3.0, 2.4))

ggsave("figures/Figure2_power_surface.pdf", fig2,
       width = 11.2, height = 6.8, device = cairo_pdf)
ggsave("figures/Figure2_power_surface.png", fig2,
       width = 11.2, height = 6.8, dpi = 300)
message("Figure 2 written.")

# =========================================================================
#  End of file.   manikmaurya.in@gmail.com   |   cognivia.vercel.app
# =========================================================================
