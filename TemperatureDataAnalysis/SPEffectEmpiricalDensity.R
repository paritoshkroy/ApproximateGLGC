# module load gcc/11.3.0
# module load r/4.2.1

rm(list=ls())
graphics.off()
library(tidyverse)

load("NNNN_GLGC_Temps.RData")
kde_df_nnnn <- kde_df %>% mutate(Method = 1) %>% select(-Key)
rm(kde_df)

load("NNHS1_GLGC_Temps.RData")
kde_df_nnhs1 <- kde_df %>% mutate(Method = 2) %>% select(-Key)
rm(kde_df)

load("NNHS2_GLGC_Temps.RData")
kde_df_nnhs2 <- kde_df %>% mutate(Method = 3) %>% select(-Key)
rm(kde_df)

load("HSHS1_GLGC_Temps.RData")
kde_df_hshs1 <- kde_df %>% mutate(Method = 4) %>% select(-Key)
rm(kde_df)

load("HSHS2_GLGC_Temps.RData")
kde_df_hshs2 <- kde_df %>% mutate(Method = 5) %>% select(-Key)

##################################################################
## Empirical kernel density and Confidence Intervals
##################################################################

# R(K) for a normal
Rk <- 1 / (2 * sqrt(pi))

# Compute the kde (NR bandwidth)
kdePm <- density(z_summary$post.mean, from = min(z_summary$post.mean), to = max(z_summary$post.mean), n = nrow(z_summary), bw = "nrd")

# Selected bandwidth
hPm <- kdePm$bw

# Estimate the variance
var_kdePm_hat <- kdePm$y * Rk / (nrow(z_summary) * hPm)

# CI with estimated variance
alpha <- 0.05
z_alpha2 <- qnorm(1 - alpha / 2)

lciPm <- kdePm$y - z_alpha2 * sqrt(var_kdePm_hat)
uciPm <- kdePm$y + z_alpha2 * sqrt(var_kdePm_hat)

# Plot estimate, CIs and expectation
kde_df_hshs2 <- tibble(x = kdePm$x, d = kdePm$y, lci = lciPm, uci = uciPm) %>% mutate(Method = 5)

kde_df <- rbind(kde_df_nnnn, kde_df_nnhs1, kde_df_nnhs2, kde_df_hshs1, kde_df_hshs2)
kde_df <- kde_df %>% mutate(MethodFactor = factor(Method, labels = c("NNNN (5)", "NNHS (5, 47)", "NNHS (5, 54)", "HSHS (47)", "HSHS (54)")))

ggplot(data = kde_df, group = MethodFactor) +
  geom_ribbon(aes(x = x, ymin = lci, ymax = uci, fill = MethodFactor), alpha = 0.35) +
  geom_line(aes(x = x, y = d, col = MethodFactor), linetype = "dashed", linewidth = 0.35) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.7),
        legend.title = element_blank())
ggsave(filename = "SPEffect_TemperatureDataAnalysis.png", height = 4, width = 8)

