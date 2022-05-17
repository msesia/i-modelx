options(width=300)

library(tidyverse)
library(latex2exp)

idir <- "results/crt_robust_1/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

###################
## Sharp methods ##
###################

results %>%
    filter(ci.type=="sharp", signal_std==1) %>%
    pivot_wider(names_from="bound", values_from="bound.value") %>%
    mutate(Signal=signal_mean*(null==FALSE)*prop.treat)


df <- results %>%
    filter(ci.type=="sharp", signal_std==1) %>%
    pivot_wider(names_from="bound", values_from="bound.value") %>%
    mutate(Signal=tau_k, Length=Upper-Lower, Coverage=(Lower<=Signal)*(Upper>=Signal)) %>%
    gather(Length, Coverage, key="Key", value="Value") %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, num_inter, delta_inter, rho, k.rel, alpha, treat.threshold, prop.treat, model.class, ci.type, j, multivariate, null, Key) %>%
    summarise(Size=n(), Value.se = 2*sd(Value)/sqrt(n()), Value=mean(Value))
df

df.nominal <- tibble(Key="Coverage", Value=0.9, n=1000, multivariate=c(TRUE,FALSE))
df.dummy <- tibble(Key="Coverage", Value=0, n=1000, multivariate=c(TRUE,FALSE), Statistics="Lasso")

pp <- df %>%
    filter(num_inter==0) %>%
    mutate(Null=ifelse(null, "Null", "Non-Null")) %>%
    mutate(Statistics=factor(model.class, c("glmnet", "glmnet.fast", "glmnet.inter"), c("Lasso", "Lasso (fast)", "Lasso (interactions)"))) %>%
    ggplot(aes(x=n, y=Value, color=Statistics, shape=Statistics)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
    geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Null, scale="free_y") +
    scale_x_continuous(trans='log10') +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_robust_0.pdf"), width=6, height=3, units="in")

pp <- df %>%
    filter(num_inter==1) %>%
    mutate(Null=ifelse(null, "Null", "Non-Null")) %>%
    mutate(Statistics=factor(model.class, c("glmnet", "glmnet.fast", "glmnet.inter"), c("Lasso", "Lasso (fast)", "Lasso (interactions)"))) %>%
    ggplot(aes(x=n, y=Value, color=Statistics, shape=Statistics)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
    geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Null, scale="free_y") +
    scale_x_continuous(trans='log10') +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_robust_1.pdf"), width=6, height=3, units="in")

pp <- df %>%
    filter(num_inter==2) %>%
    mutate(Null=ifelse(null, "Null", "Non-Null")) %>%
    mutate(Statistics=factor(model.class, c("glmnet", "glmnet.fast", "glmnet.inter"), c("Lasso", "Lasso (fast)", "Lasso (interactions)"))) %>%
    ggplot(aes(x=n, y=Value, color=Statistics, shape=Statistics)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
    geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Null, scale="free_y") +
    scale_x_continuous(trans='log10') +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_robust_2.pdf"), width=6, height=3, units="in")


######################
## Quantile methods ##
######################

df <- results %>%
    filter(ci.type=="quantile", signal_std==1) %>%
    mutate(Distance = abs(bound.value-tau_k), Miss=(bound=="lower")*(bound.value>tau_k)+(bound=="upper")*(bound.value<tau_k)>0, Coverage=!Miss) %>%
    gather(Coverage, Distance, value="Value", key="Key") %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, num_inter, delta_inter, rho, k.rel, alpha, treat.threshold, prop.treat, model.class, ci.type,
             j, multivariate, null, bound, Key) %>%
    summarize(Num.exp=n(), Value.se=2*sd(Value)/sqrt(n()), Value=mean(Value))

df.nominal <- tibble(Key="Coverage", Value=0.9, n=1000, multivariate=c(TRUE,FALSE))
df.dummy <- tibble(Key="Coverage", Value=0.1, n=1000, multivariate=c(TRUE,FALSE))

pp <- df %>%
    mutate(Value = pmax(-20, pmin(Value, 20))) %>%
    mutate(Interactions = factor(num_inter, c(0,1,2), c("No interactions", "1 interaction", "2 interactions"))) %>%
    filter(multivariate==TRUE, j==2, null==FALSE, k.rel==90, bound=="lower") %>%    
    ggplot(aes(x=n, y=Value)) +
    geom_point() +
    geom_line() +
    geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Interactions, scales="free") +
    scale_x_continuous(trans='log10') +
    ylim(0,NA) +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_robust_q_nonnull.pdf"), width=6, height=3, units="in")

pp <- df %>%
    mutate(Value = pmax(-20, pmin(Value, 20))) %>%
    mutate(Interactions = factor(num_inter, c(0,1,2), c("No interactions", "1 interaction", "2 interactions"))) %>%
    filter(multivariate==TRUE, j==1, null==TRUE, k.rel==10, bound=="upper") %>%
    ggplot(aes(x=n, y=Value)) +
    geom_point() +
    geom_line() +
    geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Interactions, scales="free") +
    scale_x_continuous(trans='log10') +
    ylim(0,NA) +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_robust_q_null.pdf"), width=6, height=3, units="in")
