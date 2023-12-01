rm(list=ls())
library(tidyverse)


load("Full_Setup1.RData")
fixed_summary <- fixed_summary %>% mutate(Method = "Full") %>% mutate(lscale = 0.1) %>% select(Method, lscale, everything())

load("Short_Results.rda")
full_fixed_summary_df <- full_fixed_summary_df %>% filter(Setup == 4)
full_fixed_summary_df <- full_fixed_summary_df %>% select(-Setup,-y) %>% select(Method,lscale, everything())

full_fixed_summary_df <- rbind(fixed_summary, full_fixed_summary_df)


fnames <- c("NNNN_Setup1.RData", "NNNN_Setup2.RData", "NNNN_Setup3.RData", "NNNN_Setup4.RData", "NNNN_Setup5.RData", "NNNN_Setup6.RData")
nnnn_fixed_summary_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  fixed_summary <- fixed_summary %>% 
    mutate(lscale = true[variable=="ell1"]) %>% 
    mutate(Method = "NNNN") %>%
    mutate(m = m) %>%
    select(Method, m, lscale, everything())
})
nnnn_fixed_summary_df <- do.call(rbind, nnnn_fixed_summary_list)


fnames <- c("NNHS_Setup1.RData","NNHS_Setup2.RData","NNHS_Setup3.RData", "NNHS_Setup4.RData", "NNHS_Setup5.RData", "NNHS_Setup6.RData", "NNHS_Setup7.RData", "NNHS_Setup8.RData")
nnhs_fixed_summary_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  fixed_summary <- fixed_summary %>% 
    mutate(lscale = true[variable=="ell1"]) %>% 
    mutate(Method = "NNHS") %>%
    mutate(m = m) %>%
    mutate(m1 = m1) %>%
    select(Method, m, m1, lscale, everything())
})
nnhs_fixed_summary_df <- do.call(rbind, nnhs_fixed_summary_list)


fnames <- c("HSHS_Setup1.RData", "HSHS_Setup2.RData", "HSHS_Setup3.RData", "HSHS_Setup4.RData", "HSHS_Setup5.RData", "HSHS_Setup6.RData"  )
hshs_fixed_summary_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  fixed_summary <- fixed_summary %>% 
    mutate(lscale = true[variable=="ell1"]) %>% 
    mutate(Method = "HSHS") %>%
    mutate(m1 = m1) %>%
    select(Method, m1, lscale, everything())
})
hshs_fixed_summary_df <- do.call(rbind, hshs_fixed_summary_list)

save(hshs_fixed_summary_df, nnhs_fixed_summary_df, nnnn_fixed_summary_df, full_fixed_summary_df, file = "collected_fixed_summary.RData")


### Prediction score
rm(list=ls())
load("Full_Setup1.RData")
scores_df <- scores_df %>% mutate(Method = "Full") %>% mutate(lscale = 0.1) %>% select(Method, lscale, everything())

load("Short_Results.rda")
full_scores_df <- full_scores_df %>% filter(Setup == 4)
full_scores_df <- full_scores_df %>% select(-Setup,-y) %>% select(Method,lscale, everything())

full_scores_df <- rbind(scores_df, full_scores_df)
full_scores_df

fnames <- c("NNNN_Setup1.RData", "NNNN_Setup2.RData", "NNNN_Setup3.RData", "NNNN_Setup4.RData", "NNNN_Setup5.RData", "NNNN_Setup6.RData")
nnnn_scores_df_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  scores_df <- scores_df %>% 
    mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"]) %>% 
    mutate(Method = "NNNN") %>%
    mutate(m = m) %>%
    select(Method, m, lscale, everything())
})
nnnn_scores_df <- do.call(rbind, nnnn_scores_df_list)


fnames <- c("NNHS_Setup1.RData","NNHS_Setup2.RData","NNHS_Setup3.RData", "NNHS_Setup4.RData", "NNHS_Setup5.RData", "NNHS_Setup6.RData", "NNHS_Setup7.RData", "NNHS_Setup8.RData")
nnhs_scores_df_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  scores_df <- scores_df %>% 
    mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"]) %>% 
    mutate(Method = "NNHS") %>%
    mutate(m = m) %>%
    mutate(m1 = m1) %>%
    select(Method, m, m1, lscale, everything())
})
nnhs_scores_df <- do.call(rbind, nnhs_scores_df_list)


fnames <- c("HSHS_Setup1.RData", "HSHS_Setup2.RData", "HSHS_Setup3.RData", "HSHS_Setup4.RData", "HSHS_Setup5.RData", "HSHS_Setup6.RData"  )
hshs_scores_df_list <- lapply(1:length(fnames), function(l) {
  load(fnames[l])
  scores_df <- scores_df %>% 
    mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"]) %>% 
    mutate(Method = "HSHS") %>%
    mutate(m1 = m1) %>%
    select(Method, m1, lscale, everything())
})
hshs_scores_df <- do.call(rbind, hshs_scores_df_list)

save(full_scores_df, nnnn_scores_df, nnhs_scores_df, hshs_scores_df, file = "collected_scores.RData")

### Graphs and tables
rm(list=ls())
graphics.off()
library(tidyverse)
load("./TuningParameterSensitivity/collected_fixed_summary.RData")
out_df <- rbind(hshs_fixed_summary_df %>% mutate(Method = paste0("HSHS",m1)) %>% select(-m1),
      nnhs_fixed_summary_df %>% mutate(Method = paste0("NN",m,"HS",m1)) %>% select(-m,-m1),
      nnnn_fixed_summary_df %>% mutate(Method = paste0("NNNN",m)) %>% select(-m),
      full_fixed_summary_df %>% mutate(Method = "Exact"))
table(out_df$Method)
out_df <- out_df %>% mutate(Method = recode(Method, `Exact` = 1, `NNNN5` = 2, `NNNN10` = 3, `NNNN15` = 4, `NN5HS22` = 5, `NN5HS32` = 6, `NN10HS22` = 7, `NN10HS32` = 8, `NN5HS42` = 9, `NN5HS52` = 10, `NN10HS42` = 11, `NN10HS52` = 12, `HSHS22` = 13, `HSHS27` = 14, `HSHS32` = 15, `HSHS42` = 16, `HSHS47` = 17, `HSHS52` = 18))
out_df <- out_df %>% mutate(Method = factor(Method, labels = c("Exact", "NNNN5", "NNNN10", "NNNN15", "NN5HS22", "NN5HS32", "NN10HS22", "NN10HS32", "NN5HS42", "NN5HS52", "NN10HS42", "NN10HS52", "HSHS22", "HSHS27", "HSHS32", "HSHS42", "HSHS47", "HSHS52")))
table(out_df$Method)

out_df <- out_df %>% 
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau")))

ggplot(out_df %>% filter(lscale == 0.1), aes(x= Method)) + 
  geom_point(aes(y = mean), shape = 1, size = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25, linewidth = 0.25) +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed) + 
  xlab("") +
  ylab("Posterior mean (95% CI)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, angle = 30, vjust = 0.5, hjust=0.5))
ggsave(filename = "./TuningParameterSensitivity/collected_fixed_summary_lscale0.1.png", height = 5, width = 11)



ggplot(out_df %>% filter(lscale == 0.4), aes(x= Method)) + 
  geom_point(aes(y = mean), shape = 1, size = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25, linewidth = 0.25) +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed) + 
  xlab("") +
  ylab("Posterior mean (95% CI)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, angle = 30, vjust = 0.5, hjust=0.5))
ggsave(filename = "./TuningParameterSensitivity/collected_fixed_summary_lscale0.4.png", height = 5, width = 11)

