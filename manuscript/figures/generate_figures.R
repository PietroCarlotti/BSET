# Load necessary libraries
library(ggplot2)
library(BSET)

# Clear console and environment
rm(list = ls())
graphics.off()
cat("\014")

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#######################
# Discrepancy heatmap #
#######################

Delta  <- 5
d_vals <- seq(-10, 10, length.out = 500)
grid   <- expand.grid(d_Y = d_vals, d_S = d_vals)

grid$discrepancy <- with(grid,
  (1/4) * abs(
    pnorm(Delta + d_Y) - pnorm(Delta + d_S) +
    pnorm(Delta - d_Y) - pnorm(Delta - d_S)
  )
)

p_heatmap <- ggplot(grid, aes(x = d_Y, y = d_S, fill = discrepancy)) +
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

ggsave("discrepancy.pdf", plot = p_heatmap, width = 6, height = 5)


#####################################
# DCCT primary outcome distribution #
#####################################



###############################
# DCCT surrogate distribution #
###############################



#####################################
# DCCT posterior without covariates #
#####################################



##################################
# DCCT posterior with covariates #
##################################


