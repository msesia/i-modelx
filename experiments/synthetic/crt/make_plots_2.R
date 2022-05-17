options(width=160)

library(tidyverse)

idir <- "results/crt_2/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

df <- results %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, prop.treat, rho, multivariate, reference, null, tail.sign, k.rel, k, delta) %>%
    summarize(Num.exp=n(), `Prob. of rejection` = mean(pval<0.05), Power.se = 2*sd(pval<0.05)/sqrt(n()))

df.nominal <- tibble(Key="`Prob. of rejection`", null=TRUE, `Prob. of rejection`=0.05)

pp <- df %>%
    filter(n==500, p==10, k.rel==90) %>%
    ggplot(aes(x=prop.treat, y=`Prob. of rejection`, color=multivariate)) +
    geom_point() +
    geom_line() +
#    geom_errorbar(aes(ymin=Power-Power.se, ymax=Power+Power.se), width=0.1) +
    facet_grid(null~reference) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
#    scale_x_continuous(trans='log10') +
#    scale_y_continuous(lim=c(0,1)) +
    theme_bw()


## Plot: multivariate vs. univariate
df.nominal <- tibble(Key="`Prob. of rejection`", null=TRUE, `Prob. of rejection`=0.05) %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null")))
pp <- df %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null"))) %>%
    mutate(Reference = factor(reference, c("treated","control","auto"), c("Treatment", "Control", "Adaptive"))) %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
    filter(n==1000, p==10, k.rel==90) %>%
    ggplot(aes(x=prop.treat, y=`Prob. of rejection`, color=Reference, linetype=Reference, shape=Reference)) +
    geom_point(alpha=0.75) +
    geom_line(alpha=0.75) +
#    geom_errorbar(aes(ymin=Power-Power.se, ymax=Power+Power.se), width=0.1) +
    facet_grid(Null~Statistics) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
#    scale_x_continuous(trans='log10') +
#    scale_y_continuous(lim=c(0,1)) +
    xlab("Proportion of treated samples") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
          )
pp %>% ggsave(file="../figures/experiment_crt_2.pdf", width=6, height=3, units="in")
