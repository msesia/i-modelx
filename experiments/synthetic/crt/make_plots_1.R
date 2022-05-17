options(width=160)

library(tidyverse)

idir <- "results/crt_1/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

df <- results %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, multivariate, reference, null, tail.sign, k.rel, k, delta) %>%
    summarize(Num.exp=n(), `Prob. of rejection` = mean(pval<0.05), Power.se = 2*sd(pval<0.05)/sqrt(n()))

df.nominal <- tibble(Key="`Prob. of rejection`", null=TRUE, `Prob. of rejection`=0.05)

color.scale <- c("tan2", "coral4", "indianred")
shape.scale <- c(1, 0, 8)

df %>%
    filter(p==10, delta==2, k.rel==90) %>%
    ggplot(aes(x=n, y=`Prob. of rejection`, color=multivariate)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=`Prob. of rejection`-Power.se, ymax=`Prob. of rejection`+Power.se), width=0.1) +
    facet_grid(null~reference) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
    scale_x_continuous(trans='log10') +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
#    scale_y_continuous(lim=c(0,1)) +
    theme_bw()


## Plot 1: multivariate vs. univariate
df.nominal <- tibble(Key="Prob. of rejection", null=TRUE, `Prob. of rejection`=0.05) %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null")))
pp <- df %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null"))) %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
    filter(p==10, k.rel==90, delta==2, reference=="auto") %>%
    ggplot(aes(x=n, y=`Prob. of rejection`, color=Statistics)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=`Prob. of rejection`-Power.se, ymax=`Prob. of rejection`+Power.se), width=0.1) +
    facet_grid(.~Null) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
    scale_x_continuous(trans='log10') +
#    scale_y_continuous(lim=c(0,1)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("Sample size") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp

pp %>% ggsave(file="../figures/experiment_crt_1.pdf", width=5, height=2, units="in")

## Plot 3:

df <- results %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, multivariate, reference, null, tail.sign, k.rel, k, delta, j) %>%
    summarize(Num.exp=n(), `Prob. of rejection` = mean(pval<0.05), Power.se = 2*sd(pval<0.05)/sqrt(n()))

delta.plot <- c(0, 1, 2)
pp <- df %>%
    filter(p==10, n==1000, delta %in% delta.plot, j%in%c(1,2), reference=="auto") %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null"))) %>%
    mutate(delta = sprintf("c: %d", delta)) %>%
    mutate(delta = factor(delta, sprintf("c: %d", delta.plot))) %>%
    ggplot(aes(x=k.rel, y=`Prob. of rejection`, color=Statistics)) +
    geom_point() +
    geom_line() +
#    geom_errorbar(aes(ymin=Power-Power.se, ymax=Power+Power.se), width=0.1) +
    facet_grid(Null~delta) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
#    scale_x_continuous(trans='log10') +
#    scale_y_continuous(lim=c(0,1)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("100*q/n") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp

pp %>% ggsave(file="../figures/experiment_crt_3.pdf", width=7, height=3, units="in")


## Plot 4:

df <- results %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, multivariate, reference, null, tail.sign, k.rel, k, delta, j) %>%
    summarize(Num.exp=n(), `Prob. of rejection` = mean(pval<0.05), Power.se = 2*sd(pval<0.05)/sqrt(n()))

df.nominal <- tibble(Key="`Prob. of rejection`", null=TRUE, `Prob. of rejection`=0.05) %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null")))

n.plot <- c(100,200,500,1000)
null.threshold <- 5.3
pp <- df %>%
    filter(p==10, n %in% n.plot, k.rel==90, j==2, reference=="auto") %>%
    mutate(Null = factor(null, c(TRUE,FALSE), c("Null", "Non-Null"))) %>%
    mutate(n = sprintf("Sample size: %d", n)) %>%
    mutate(n = factor(n, sprintf("Sample size: %d", n.plot))) %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
    ggplot(aes(x=delta, y=`Prob. of rejection`, color=Statistics)) +
    geom_point() +
    geom_line() +
#    geom_errorbar(aes(ymin=Power-Power.se, ymax=Power+Power.se), width=0.1) +
    facet_grid(.~n) +
    geom_hline(data=df.nominal, aes(yintercept=`Prob. of rejection`), linetype=2) +
    geom_vline(xintercept=null.threshold, linetype=2) +
#    scale_x_continuous(trans='log10') +
#    scale_y_continuous(lim=c(0,1)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("c") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10))
pp

pp %>% ggsave(file="../figures/experiment_crt_4.pdf", width=7, height=2, units="in")


## method.levels <- c("Exact", "Naive", "Split", "Vanilla")
## method.labels <- c("HAKF", "Naive", "Split", "Vanilla")
## color.scale <- c("blue2", "#009E73", "#CC79A7", "orange")
## shape.scale <- c(15, 4, 8, 1)

## df.nominal <- tibble(Key="FDR", Value=0.1, Method="Vanilla")
## df.ghost.1 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=0, n=1000, Method="Vanilla")
## df.ghost.2 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=1, n=1000, Method="Vanilla")
## df.ghost <- rbind(df.ghost.1, df.ghost.2)

## keys.levels <- c("FDR", "Power", "Causal.prop", "Beta.sd")
## keys.labels <- c("FDR", "Power", "Homeg. (non-null prop.)", "Heterog. (effect s.d.)")
## df.ghost <- df.ghost %>%
##     mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
##            Method=factor(Method, levels=method.levels, labels=method.labels))
## df.nominal <- df.nominal %>%
##     mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
##            Method=factor(Method, levels=method.levels, labels=method.labels))
    
## pp <- df %>%
##     filter(delta_inter==0.05, signal_mean==4, Key %in% keys.levels) %>%
##     mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
##            Method=factor(Method, levels=method.levels, labels=method.labels)) %>%
##     ggplot(aes(x=n, y=Value, color=Method, shape=Method)) +
##     geom_point(alpha=0.5) +
##     geom_line(alpha=0.5) +
##     geom_point(alpha=0, data=df.ghost, aes(x=n,y=Value)) +
## #    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
##     geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
##     facet_wrap(~Key, scales="free", nrow=1) +
##     scale_x_continuous(trans='log10') +
##     scale_color_manual(values=color.scale) +
##     scale_shape_manual(values=shape.scale) +
##     xlab("Sample size") +
##     ylab("") +
##     theme_bw() +
##     theme(legend.position = "bottom",
##           plot.title = element_text(size = 10),
##           axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
##           strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
##           legend.title=element_text(size = 10))
## pp

## pp %>% ggsave(file="../figures/experiment_heterogeneous_1.pdf", width=8, height=2.5, units="in")
