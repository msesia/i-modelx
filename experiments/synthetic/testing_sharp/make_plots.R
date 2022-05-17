options(width=160)

library(tidyverse)

idir <- "results/estimation_sharp/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

df <- results %>%
    filter(Status=="converged") %>%
    mutate(method = as.factor(method),
           fast = as.factor(fast),
           Coverage=(Lower<=ate)*(Upper>=ate), Length=Upper-Lower) %>%
    gather(Coverage, Length, key="Key", value="Value") %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, interaction_strength, rho, method, fast, Key) %>%
    summarize(Value.se=2*sd(Value)/sqrt(n()), Value=mean(Value))

df.nominal <- tibble(Key="Coverage", Value=0.9)

df %>%
    filter(signal_mean==2, signal_std==0) %>%
    ggplot(aes(x=n, y=Value, color=fast)) +
    geom_point(alpha=0.5) +
    geom_line(alpha=0.5) +
    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se)) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~interaction_strength, labeller=labeller(.default = label_both), scales="free") +
    scale_x_continuous(trans='log10') +
    theme_bw()


df %>%
    filter(signal_mean==2, interaction_strength==0) %>%
    ggplot(aes(x=signal_std, y=Value, color=fast)) +
    geom_point(alpha=0.5) +
    geom_line(alpha=0.5) +
    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se)) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~n, labeller=labeller(.default = label_both), scales="free") +
    scale_x_continuous(trans='log10') +
    theme_bw()
