# =============================================================================
# Figure | Genetic correlation heatmap (LDSC vs HDL)
# =============================================================================
# Split-cell heatmap of genome-wide genetic correlation (rg) between the
# three migraine phenotypes and seven cardiovascular diseases. Each cell is
# divided in two: the left half shows the LDSC estimate, the right half the
# HDL estimate. Significance is annotated with asterisks.
#
# This script is self-contained: the rg / P values below are the published
# estimates (Table 1 of the manuscript). Only public packages are required.
#   install.packages(c("ggplot2", "dplyr"))
# =============================================================================

library(ggplot2)
library(dplyr)

OUT_DIR <- "results/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Data: published rg and P values ---------------------------------------
dat <- data.frame(
  migraine = c(rep("Migraine", 7), rep("MO", 7), rep("MA", 7)),
  cvd = rep(c("AF", "CAD", "HF", "HTN", "IS", "MI", "PAD"), 3),
  ldsc_rg = c(0.2006, 0.3042, 0.1231, 0.2733, 0.1949, 0.1914, 0.2563,
              0.0600, 0.1746, 0.1634, 0.2858, 0.1156, 0.1530, 0.3174,
              0.0931, 0.2059, 0.1555, 0.2471, 0.3192, 0.2073, 0.3056),
  ldsc_p  = c(0.0085, 8.0441e-06, 0.37, 4.6604e-05, 0.2245, 0.0271, 0.1191,
              0.0725, 5.1099e-06, 0.0278, 8.9713e-14, 0.0712, 0.0023, 0.0015,
              0.0159, 1.9183e-08, 0.0496, 4.9095e-10, 3.5469e-07, 2.4206e-06, 0.0035),
  hdl_rg  = c(0.0654, 0.1159, 0.1491, 0.1122, 0.0519, 0.0947, 0.1245,
              0.0851, 0.2415, 0.2625, 0.3407, 0.1534, 0.2452, 0.3377,
              0.0851, 0.2485, 0.2207, 0.2610, 0.3671, 0.2983, 0.3033),
  hdl_p   = c(9.23e-02, 1.60e-03, 4.59e-03, 5.88e-03, 5.03e-01, 2.02e-01, 7.93e-02,
              4.89e-02, 1.33e-07, 1.70e-03, 1.23e-14, 3.72e-02, 6.55e-03, 4.57e-03,
              4.89e-02, 2.14e-08, 5.67e-04, 1.50e-09, 2.11e-05, 5.19e-04, 4.81e-03),
  stringsAsFactors = FALSE
)

# Bonferroni threshold: 0.05 / 7 cardiovascular traits
bonf <- 0.05 / 7

# significance annotation: ** Bonferroni, * nominal
add_sig <- function(p) ifelse(p < bonf, "**", ifelse(p < 0.05, "*", ""))

# ---- Reshape to long format: one row per migraine x cvd x method ------------
dat_long <- bind_rows(
  dat %>% transmute(migraine, cvd, method = "LDSC",
                    rg = ldsc_rg, p = ldsc_p,
                    label = paste0(sprintf("%.2f", ldsc_rg), add_sig(ldsc_p))),
  dat %>% transmute(migraine, cvd, method = "HDL",
                    rg = hdl_rg, p = hdl_p,
                    label = paste0(sprintf("%.2f", hdl_rg), add_sig(hdl_p)))
)

dat_long$migraine <- factor(dat_long$migraine, levels = c("Migraine", "MO", "MA"))
dat_long$cvd      <- factor(dat_long$cvd,
                            levels = c("HTN", "CAD", "MI", "AF", "IS", "HF", "PAD"))
dat_long$method   <- factor(dat_long$method, levels = c("LDSC", "HDL"))

# ---- Numeric coordinates for the left / right half-cells --------------------
dat_long <- dat_long %>%
  mutate(
    x_num = as.numeric(cvd),
    y_num = as.numeric(migraine),
    x_pos = ifelse(method == "LDSC", x_num - 0.25, x_num + 0.25),  # LDSC=left
    text_col = ifelse(rg > 0.22, "white", "grey20")               # readable text
  )

# ---- Plot -------------------------------------------------------------------
p <- ggplot(dat_long) +
  # left / right half-cell colour blocks
  geom_tile(aes(x = x_pos, y = y_num, fill = rg),
            width = 0.48, height = 0.95, color = NA) +
  # cell outline
  geom_tile(aes(x = x_num, y = y_num),
            width = 0.98, height = 0.95, fill = NA,
            color = "white", linewidth = 0.8) +
  # central divider between the two methods
  geom_segment(data = dat_long %>% filter(method == "LDSC"),
               aes(x = x_num, xend = x_num,
                   y = y_num - 0.475, yend = y_num + 0.475),
               color = "white", linewidth = 0.5) +
  # rg value + significance asterisks
  geom_text(aes(x = x_pos, y = y_num, label = label, color = text_col),
            size = 2.8, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(0, 0.40),
                       name = expression(italic(r)[g])) +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:7, labels = levels(dat_long$cvd),
                     expand = c(0.06, 0.06)) +
  scale_y_continuous(breaks = 1:3, labels = levels(dat_long$migraine),
                     expand = c(0.12, 0.12), trans = "reverse") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid  = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(30, 10, 30, 10)
  )

# ---- Method legend (LDSC = left, HDL = right) and significance note ---------
legend_df <- data.frame(
  x = c(3.0, 4.5), y = c(0.25, 0.25),
  label = c("LDSC (left)", "HDL (right)"),
  fill  = c("#7a9ec4", "#c47a7a")
)

p2 <- p +
  geom_rect(data = legend_df,
            aes(xmin = x - 0.15, xmax = x + 0.15,
                ymin = y - 0.08, ymax = y + 0.08),
            fill = legend_df$fill, inherit.aes = FALSE) +
  geom_text(data = legend_df,
            aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 3, color = "grey40",
            fontface = "bold", inherit.aes = FALSE) +
  annotate("text", x = 4, y = 3.7,
           label = "*P < 0.05    **P < 0.05/7 (Bonferroni)",
           size = 3, color = "grey50") +
  coord_cartesian(ylim = c(3.8, 0.15), clip = "off")

# ---- Save -------------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure_rg_heatmap.pdf"), p2,
       width = 8, height = 4.5, dpi = 300)
ggsave(file.path(OUT_DIR, "Figure_rg_heatmap.png"), p2,
       width = 8, height = 4.5, dpi = 300, bg = "white")

cat("Saved: Figure_rg_heatmap.pdf / .png in", OUT_DIR, "\n")
