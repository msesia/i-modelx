options(width=160)

library(tidyverse)

idir <- "results/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

## Store combined results
if(FALSE) {
    out.file <- "results_release/estimation_linear.txt"
    results %>% write_delim(out.file, delim=" ")
}

if(FALSE) {
    results %>%
        ggplot(aes(x=Method, y=Estimate)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=Beta)) +
        facet_grid(n~j) +
        theme_bw()
}

df <- results %>%
    mutate(method = as.factor(method)) %>%
    group_by(n, p, p_causal, signal_mean, signal_std, rho, j, causal, k, tail, tau_k,
             null, method, multivar, seed, seed.model) %>%
    filter(p==100, signal_mean==2, signal_std==1, k==750, delta==1)


df %>%
    ggplot(aes(x=pval)) +
    geom_histogram() +
    facet_wrap(multivar~causal, labeller=labeller(.default = label_both)) +
    theme_bw()


# Measure FDR, Power
df2 <- results %>%
    mutate(method = as.factor(method)) %>%
    group_by(n, p, p_causal, signal_mean, signal_std, rho, k, delta, tail,
             method, multivar, seed, seed.model) %>%
    mutate(pval.a = p.adjust(pval, method="BH"),
           rejected = pval.a<=0.1) %>%   
    summarise(FDP = ifelse(sum(rejected)>0,1-mean(causal[rejected]),0),
              Power = ifelse(sum(rejected)>0, sum(causal[rejected])/sum(causal), 0)) %>%
    group_by(n, p, p_causal, signal_mean, signal_std, rho, k, delta, tail,
             method, multivar) %>%
    summarise(FDR = mean(FDP), Power=mean(Power))
    

df2 %>%
    ggplot(aes(x=delta, y=Power, color=multivar, linetype=multivar)) +
    geom_point() +
    geom_line() +
    facet_grid(k~signal_mean) +
    xlim(0,7) +
    theme_bw()


df2 %>%
    ggplot(aes(x=delta, y=FDR, color=multivar, linetype=multivar)) +
    geom_point() +
    geom_line() +
    facet_grid(k~signal_mean) +
    geom_hline(yintercept=0.1) +
    xlim(0,7) +
    theme_bw()

