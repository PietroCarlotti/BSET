# Load necessary libraries
library(dplyr)
library(ggplot2)
library(BSET)

# Clear console and environment
rm(list = ls())
graphics.off()
cat("\014")

# Load the data
hb <- read.csv("~/Desktop/dcct_data_more.csv")

# Filter and wrangle the data
hb_high_risk <- hb %>%
  filter(
    PRIMARY == 0,
    SMOKER == 1,
    BMI.x >= 25.0
  ) %>%
  rename(
    Y = CHANGE_HBA18,
    S = CHANGE_HBA06
  ) %>%
  filter(
    !is.na(Y),
    !is.na(S)
  ) %>%
  mutate(
    Z = if_else(GROUP == "EXPERIMENTAL", 1, 0),
    MALE = if_else(SEX == "M", 1, 0)
  )

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###################################################
# Theta vs delta distance plot - Binary covariate #
###################################################

# Parameters
Delta <- 5

# Grid
d_vals <- seq(-10, 10, length.out = 500)
grid_X_binary <- expand.grid(
  d_Y = d_vals,
  d_S = d_vals
)

# Distance
grid_X_binary$distance <- with(grid_X_binary,
  (1/4) * abs(
    pnorm(Delta + d_Y) - pnorm(Delta + d_S) +
    pnorm(Delta - d_Y) - pnorm(Delta - d_S)
  )
)

# Distance plot
distance_plot_X_binary <- ggplot(grid_X_binary, mapping = aes(x = d_Y, y = d_S, fill = distance)) +
  geom_raster() +
  scale_fill_gradientn(
    name = expression("|" * theta - delta * "|"),
    colors = c("#08306B", "#2171B5", "#238B45", "#CCCC00", "#D94801", "#A50F15"),
    limits = c(0, 0.25)
  ) +
  labs(x = expression(d[Y]), y = expression(d[S])) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(
    legend.key.height = unit(1.2, "cm"),
    panel.grid = element_blank()
  )

ggsave(
  filename = "distance_plot_X_binary.pdf",
  plot = distance_plot_X_binary,
  width = 6,
  height = 5
)

####################################################
# Discrepancy heatmap - Gaussian covariate setting #
####################################################

# Parameters
Delta <- 2
m <- 2
v <- 2

# Grid
beta_vals <- seq(-12, 9, length.out = 500)
grid_X_Gaussian <- expand.grid(
  beta_Y0 = beta_vals,
  beta_S0 = beta_vals
)

# Distance
h_Delta_m_v <- function(x) {
  pnorm(m * Delta / sqrt(1 + v^2 * ((x + Delta)^2 + x^2)))
}

grid_X_Gaussian$discrepancy <- abs(h_Delta_m_v(grid_X_Gaussian$beta_Y0) - h_Delta_m_v(grid_X_Gaussian$beta_S0))

sup_distance <- abs(pnorm(m * Delta / sqrt(1 + v^2 * Delta^2 / 2)) - 0.5)

# Distance plot
distance_plot_X_Gaussian <- ggplot(grid_X_Gaussian, aes(x = beta_Y0, y = beta_S0, fill = discrepancy)) +
  geom_raster() +
  scale_fill_gradientn(
    name = expression("|" * theta - delta * "|"),
    colors = c("#08306B", "#2171B5", "#238B45", "#CCCC00", "#D94801", "#A50F15"),
    limits = c(0, sup_distance)
  ) +
  labs(
    x = expression(beta[Y[0]]),
    y = expression(beta[S[0]])
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(
    legend.key.height = unit(1.2, "cm"),
    panel.grid = element_blank()
  )

ggsave(
  filename = "distance_plot_X_Gaussian.pdf",
  plot = distance_plot_X_Gaussian,
  width = 6,
  height = 5
)

######################################################
# DCCT primary outcome distribution by treatment arm #
######################################################

plot_data <- hb_high_risk %>%
  mutate(Treatment = if_else(Z == 1, "Treatment", "Control"))

p_Y <- ggplot(plot_data, aes(x = Y, color = Treatment, fill = Treatment)) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.2,
    binwidth = 0.5,
    boundary = 0
  ) +
  scale_color_manual(values = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  scale_fill_manual(values  = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  labs(
    x = "Change in HbA1c at 4.5 years",
    y = "Frequency",
    color = NULL,
    fill  = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "DCCT_Y.pdf",
  plot = p_Y,
  width = 6,
  height = 4.5,
  device = cairo_pdf
)

################################################
# DCCT surrogate distribution by treatment arm #
################################################

p_S <- ggplot(plot_data, aes(x = S, color = Treatment, fill = Treatment)) +
  geom_histogram(
    mapping = aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.2,
    binwidth = 0.5,
    boundary = 0
  ) +
  scale_color_manual(values = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  scale_fill_manual(values  = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  labs(
    x = "Change in HbA1c at 1.5 years",
    y = "Frequency",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "DCCT_S.pdf",
  plot = p_S,
  width = 6,
  height = 4.5,
  device = cairo_pdf
)

############################
# DCCT S vs Y scatter plot #
############################

scatter_data <- hb_high_risk %>%
  mutate(Treatment = if_else(Z == 1, "Treatment", "Control"))

p_scatter <- ggplot(scatter_data, aes(x = S, y = Y, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(
    formula = y ~ x,
    method = "loess",
    se = FALSE,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c("Treatment" = "#0072B2", "Control" = "#D55E00"),
    guide  = guide_legend(override.aes = list(linetype = 0))
  ) +
  labs(
    x = "Change in HbA1c at 1.5 years",
    y = "Change in HbA1c at 4.5 years",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "DCCT_scatterplot.pdf",
  plot = p_scatter,
  width = 7.5,
  height = 5,
  device = cairo_pdf
)

#################
# DCCT CDF plot #
#################

cdf_data <- hb_high_risk %>%
  mutate(Treatment = if_else(Z == 1, "Treatment", "Control")) %>%
  tidyr::pivot_longer(
    cols = c(Y, S),
    names_to = "Outcome",
    values_to = "value"
  ) %>%
  mutate(
    Outcome = if_else(Outcome == "Y", "Primary outcome", "Surrogate"),
    Group = paste(Treatment, Outcome, sep = " — ")
  )

p_cdf <- ggplot(cdf_data, aes(x = value, color = Group, linetype = Group)) +
  stat_ecdf(linewidth = 0.8) +
  scale_color_manual(values = c(
    "Treatment — Primary outcome" = "#0072B2",
    "Control — Primary outcome" = "#0072B2",
    "Treatment — Surrogate" = "#D55E00",
    "Control — Surrogate" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "Treatment — Primary outcome" = "solid",
    "Control — Primary outcome" = "dashed",
    "Treatment — Surrogate" = "solid",
    "Control — Surrogate" = "dashed"
  )) +
  labs(
    x = "Change in HbA1c",
    y = "Empirical CDF",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.91, "cm")
  )

ggsave(
  filename = "DCCT_CDF.pdf",
  plot = p_cdf,
  width = 7.5,
  height = 5,
  device = cairo_pdf
)

#####################################
# DCCT posterior without covariates #
#####################################

set.seed(1)

# BSET_no_X <- BSET(
#   data = hb_high_risk,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   beta = 0.3,
#   parallel = FALSE,
#   plot = TRUE
# )
# 
# ggsave(
#   filename = "DCCT_no_X_theta_posterior_plot.pdf",
#   plot = BSET_no_X$theta_posterior_plot,
#   width = 10,
#   height = 8,
#   device = cairo_pdf
# )

##################################
# DCCT posterior with covariates #
##################################

set.seed(1)

# BSET_X <- BSET(
#   data = hb_high_risk,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   X = c("AGE","MALE"),
#   beta = 0.3,
#   plot = TRUE
# )
# 
# ggsave(
#   filename = "DCCT_X_theta_posterior_plot.pdf",
#   plot = BSET_X$theta_posterior_plot,
#   width = 10,
#   height = 8,
#   device = cairo_pdf
# )
