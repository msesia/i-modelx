options(width=160)

library(tidyverse)

idir <- "results/testing_transfer/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

df <- results %>%
    mutate(Method = as.factor(Method), FDR=FDP) %>%
    gather(FDR, Power, True, Causal.prop, Beta.sd, key="Key", value="Value") %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, Method, Key, Test, significance, r, num_cluster) %>%
    summarize(Value.se=2*sd(Value,na.rm=T)/sqrt(n()), Value=mean(Value,na.rm=T)) %>%
    filter(!((Method=="Vanilla")&(r!=1)))

r.plot <- 3
significance.plot <- "accumulation"
#significance.plot <- "seqstep-0.4"
signal_mean.plot <- 4

method.levels <- c("Exact", "Naive", "Split", "Vanilla")
method.labels <- c("RSKF", "Naive", "Split", "Vanilla")
color.scale <- c("blue2", "#009E73", "#CC79A7", "orange")
shape.scale <- c(15, 4, 8, 1)

df.nominal <- tibble(Key="FDR", Value=0.1, Method="Vanilla")
df.ghost.1 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=0, n=1000, Method="Vanilla")
df.ghost.2 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=1, n=1000, Method="Vanilla")
df.ghost <- rbind(df.ghost.1, df.ghost.2)

keys.levels <- c("FDR", "Power", "Causal.prop", "Beta.sd")
keys.labels <- c("FDR", "Power", "Homeg. (non-null prop.)", "Heterog. (effect s.d.)")
df.ghost <- df.ghost %>%
    mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
           Method=factor(Method, levels=method.levels, labels=method.labels))
df.nominal <- df.nominal %>%
    mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
           Method=factor(Method, levels=method.levels, labels=method.labels))

pp <- df %>%
    mutate(Test = factor(Test, c(FALSE, TRUE), c("Within population", "Covariate shift"))) %>%
    filter((r==r.plot)|((Method=="Vanilla")&(r==1)), significance %in% c("KF", significance.plot)) %>%
    filter(delta_inter==0.05, signal_mean==signal_mean.plot, Key %in% keys.levels) %>%
    mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
           Method=factor(Method, levels=method.levels, labels=method.labels)) %>%
    ggplot(aes(x=n, y=Value, color=Method, shape=Method)) +
    geom_point(alpha=1) +
    geom_line(alpha=1) +
    geom_point(alpha=0, data=df.ghost, aes(x=n,y=Value)) +
#    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    facet_wrap(Test~Key, scales="free", nrow=2) +
    scale_x_continuous(trans='log10') +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))

pp %>% ggsave(file="../figures/experiment_transfer_1.pdf", width=8, height=4.5, units="in")
