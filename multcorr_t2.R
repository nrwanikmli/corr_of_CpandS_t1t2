library(tidyverse)
library(GGally)  # For ggpairs function
library(ggplot2) # For ggplot2 functionality
theme_set(theme_bw())
# Data for Entropy and Heat Capacity of 30 benzenoids hydrocarbons
x1 <- c(269.722, 334.155, 389.475, 395.882, 444.724, 447.437, 457.958, 455.839, 
        450.418, 399.491, 499.831, 513.857, 508.537, 507.395, 506.076, 512.523, 
        500.734, 520.307, 509.210, 513.879, 511.770, 509.611, 462.545, 463.738, 
        468.712, 555.409, 472.295, 554.784, 468.796, 551.708) #Entropy
x2 <- c(83.019, 133.325, 184.194, 183.654, 235.165, 233.497, 234.568, 234.638, 
        233.558, 200.815, 286.182, 285.056, 284.037, 284.088, 285.148, 284.595, 
        284.870, 284.503, 284.785, 284.740, 284.233, 284.552, 251.175, 250.568, 
        251.973, 336.098, 267.543, 337.204, 285.041, 368.518) #Heat Capacity
# 30x3 array for number of edges in each dudv partition: (2,2), (2,3), (3,3)
coef_y <- t(matrix(c(
  # (2,2)
  6, 6, 6, 7, 6, 8, 7, 8, 9, 6, 6, 7, 8, 8, 7, 10, 9, 9, 9, 8, 9, 8, 8, 8, 7, 9, 7, 6, 6, 6,
  # (2,3)
  0, 4, 8, 6, 12, 8, 10, 8, 6, 8, 16, 14, 12, 12, 14, 8, 10, 10, 10, 12, 10, 12, 8, 8, 10, 14, 10, 20, 12, 16,
  # (3,3)
  0, 1, 2, 3, 3, 5, 4, 5, 6, 5, 4, 5, 6, 6, 5, 8, 7, 7, 7, 6, 7, 6, 8, 8, 7, 8, 10, 5, 12, 19
), ncol = 30, byrow = TRUE))
# Number of vertices (v_exp)
v_exp <- c(6, 10, 14, 14, 18, 18, 18, 18, 18, 16, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 20, 20, 20, 26, 22, 26, 24, 32)
# List of index-computing functions
t_2 <- t(cbind(
  4 / (v_exp - 2)^2,
  6 / ((v_exp - 2) * (v_exp - 3)),
  9 / (v_exp - 3)^2
))
# The 'y' function will now perform the matrix multiplication
y <- function(a = 1) {
  alphaone <- t_2 ^ a
  as.numeric(coef_y %*% rowSums(alphaone))
}
# Create the multiple correlation function (depends only on a)
mult_rho <- function(.a = 1) {
  dat <- tibble(y = y(a = .a), x1 = x1, x2 = x2)
  fit <- lm(y ~ x1 + x2, data = dat)
  sqrt(summary(fit)$r.squared)
}
# Objective function for maximization
myfun <- function(a) {
  rho_value <- mult_rho(a)
  if (is.nan(rho_value) || rho_value == 0) {
    return(Inf)  # Return a large value to avoid issues with log(0)
  }
  return(-log(1 + rho_value))  # Negate for maximization
}
# Optimizing using nlminb
res <- nlminb(-abs(rnorm(1)), myfun)
print(res)
# The optimal value of 'a' and multiple correlation
alpha_hat <- res$par
rho_hat <- mult_rho(alpha_hat)
# The optimal values of y in a table
dat <- tibble(y = y(a = res$par), x1 = x1, x2 = x2)
# Plotting with GGpairs
GGally::ggpairs(dat,
                lower = list(continuous = wrap("smooth", method = "lm", se = FALSE)),
                diag = list(continuous = wrap("densityDiag"))
)
# Plot the relationship of 'alpha' and multiple correlation
plot_df <- tibble(a = seq(-20, 20, length = 500)) %>%
  mutate(mcor = unlist(purrr::map(a, myfun)))
ggplot(plot_df, aes(a, exp(-mcor) - 1)) +
  geom_line() +
  geom_point(x = alpha_hat, y = mult_rho(alpha_hat)) +
  annotate("text", x = alpha_hat + 2.2, y = mult_rho(alpha_hat) + 0.015, 
           label = paste0("(", round(alpha_hat, 3), ", ", round(rho_hat, 3), ")")) +
  geom_vline(xintercept = alpha_hat, linetype = "dashed", col ="red3") +
  labs(x = expression(alpha), y = expression("Multiple correlation "~R(alpha))) +
  theme_bw()