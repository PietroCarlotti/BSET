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

#######################
# Discrepancy heatmap #
#######################

Delta <- 5
d_vals <- seq(-10, 10, length.out = 500)
grid <- expand.grid(
  d_Y = d_vals,
  d_S = d_vals
)

grid$discrepancy <- with(grid,
  (1/4) * abs(
    pnorm(Delta + d_Y) - pnorm(Delta + d_S) +
    pnorm(Delta - d_Y) - pnorm(Delta - d_S)
  )
)

discrepancy_heatmap <- ggplot(grid, aes(x = d_Y, y = d_S, fill = discrepancy)) +
  geom_raster() +
  scale_fill_gradientn(
    name   = expression("|" * theta - delta * "|"),
    colors = c("#08306B", "#2171B5", "#238B45", "#CCCC00", "#D94801", "#A50F15"),
    limits = c(0, 0.25)
  ) +
  labs(x = expression(d[Y]), y = expression(d[S])) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(legend.key.height = unit(1.2, "cm"))

ggsave(
  filename = "discrepancy_heatmap.pdf",
  plot = discrepancy_heatmap,
  width = 6,
  height = 5
)


######################################################
# DCCT primary outcome distribution by treatment arm #
######################################################

plot_data <- hb_high_risk %>%
  mutate(Treatment = if_else(Z == 1, "Treatment", "Control"))

p_Y <- ggplot(plot_data, aes(x = Y, color = Treatment, fill = Treatment)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 alpha = 0.2, binwidth = 0.5, boundary = 0) +
  scale_color_manual(values = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  scale_fill_manual(values  = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  labs(
    x     = "Change in HbA1c from baseline to 4.5 years",
    y     = "Density",
    color = NULL,
    fill  = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "DCCT_Y.pdf",
  plot     = p_Y,
  width    = 6,
  height   = 4.5,
  device   = cairo_pdf
)

################################################
# DCCT surrogate distribution by treatment arm #
################################################

p_S <- ggplot(plot_data, aes(x = S, color = Treatment, fill = Treatment)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 alpha = 0.2, binwidth = 0.5, boundary = 0) +
  scale_color_manual(values = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  scale_fill_manual(values  = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  labs(
    x     = "Change in HbA1c from baseline to 1.5 years",
    y     = "Density",
    color = NULL,
    fill  = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "DCCT_S.pdf",
  plot     = p_S,
  width    = 6,
  height   = 4.5,
  device   = cairo_pdf
)

#############################
# DCCT S vs Y scatter plot  #
#############################

scatter_data <- hb_high_risk %>%
  mutate(Treatment = if_else(Z == 1, "Treatment", "Control"))

p_scatter <- ggplot(scatter_data, aes(x = S, y = Y, color = Treatment, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(formula = y ~ x, method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_color_manual(values = c("Treatment" = "#0072B2", "Control" = "#D55E00")) +
  scale_shape_manual(values = c("Treatment" = 16, "Control" = 17)) +
  labs(
    x = "Surrogate (S): Change in HbA1c at 1.5 years",
    y = "Primary outcome (Y): Change in HbA1c at 4.5 years",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave(
  filename = "DCCT_scatter.pdf",
  plot = p_scatter,
  width = 6,
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
