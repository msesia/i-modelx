options(width=200)

library(tidyverse)
library(latex2exp)

idir <- "results/crtci_1/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
}))

df <- results %>%
    filter(bound=="lower") %>%
    mutate(alpha=0.1) %>%
    mutate(Distance = abs(bound.value-tau_k), Miss=(bound=="lower")*(bound.value>tau_k)+(bound=="upper")*(bound.value<tau_k)>0, Coverage=!Miss) %>%
    mutate(Distance = pmax(-20, pmin(Distance, 20))) %>%
    gather(Coverage, Distance, value="Value", key="Key") %>%
    group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, prop.treat, multivariate, reference, j, null, bound, k.rel, alpha, Key) %>%
    summarize(Num.exp=n(), Value.se=2*sd(Value)/sqrt(n()), Value=mean(Value))

color.scale <- c("tan2", "coral4", "indianred")
shape.scale <- c(1, 0, 8)

for(prop.treat.plot in c(0.24,0.76)) {

    ## Plot 1: length and coverage vs sample size
    df.dummy <- tibble(Key="Coverage", Value=0.5, n=1000, multivariate=c(FALSE,TRUE)) %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    
    df.nominal <- tibble(Key="Coverage", Value=0.9, n=1000, multivariate=c(FALSE,TRUE)) %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    

    pp <- df %>%
        filter(j %in% c(1,3), signal_mean==2) %>%
        filter(p==20, alpha==0.1, k.rel==90, prop.treat==prop.treat.plot, reference=="auto") %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate","Multivariate"))) %>%
        mutate(Null = factor(null, c(FALSE,TRUE), c("Null", "Non-Null"))) %>%    
        ggplot(aes(x=n, y=Value, color=Statistics, shape=Statistics)) +
        geom_point() +    
        geom_line() +
        geom_point(data=df.dummy, aes(x=n, y=Value), alpha=0) +
        geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
        facet_wrap(Key~Null, scale="free") +
        geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
        scale_x_continuous(trans='log10') +
        scale_y_continuous(trans='log10') +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        xlab("Sample size") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "right",
              plot.title = element_text(size = 10),
              axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
              strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
              legend.title=element_text(size = 10))
    pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_1_%s.pdf", prop.treat.plot), width=7, height=4.5, units="in")


    ## Plot 2: length and coverage vs desired quantile (for upper and lower bounds)

    df.dummy <- tibble(Key="Coverage", Value=0.5, k.rel=90, multivariate=c(FALSE,TRUE)) %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    
    df.nominal <- tibble(Key="Coverage", Value=0.9, k.rel=90, multivariate=c(FALSE,TRUE)) %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    

    df <- results %>%
        filter(j==3) %>%
        mutate(alpha=0.1) %>%
        mutate(Distance = abs(bound.value-tau_k), Miss=(bound=="lower")*(bound.value>tau_k)+(bound=="upper")*(bound.value<tau_k)>0, Coverage=!Miss) %>%
                                        #    mutate(Distance = pmax(-20, pmin(Distance, 20))) %>%
        gather(Coverage, Distance, value="Value", key="Key") %>%
        group_by(n, p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, prop.treat, multivariate, reference, j, null, bound, k.rel, alpha, tau_k, Key) %>%
        summarize(Num.exp=n(), Value.se=2*sd(Value)/sqrt(n()), Value=mean(Value))

    pp <- df %>%
        mutate(Value = pmax(-20, pmin(Value, 20))) %>%
        filter(j %in% c(3), prop.treat==prop.treat.plot, Key=="Distance") %>%
        filter(p==20, n==2000, alpha==0.1, reference=="auto") %>%
        mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
        mutate(Bound = factor(bound, c("lower","upper"), c("Lower", "Upper"))) %>%    
        mutate(Effect = factor(signal_mean, c("-2","2"), c("ATE: -2", "ATE: +2"))) %>%    
        mutate(prop.treat = sprintf("Proportion treated: %d%%", prop.treat*100)) %>%
        ggplot(aes(x=k.rel, y=Value, color=Statistics, shape=Statistics)) +
        geom_point() +    
        geom_line() +
        geom_point(data=df.dummy, aes(x=k.rel, y=Value), alpha=0) +
                                        #    geom_errorbar(aes(ymin=Value-Value.se, ymax=Value+Value.se), width=0.1) +
        facet_wrap(Bound~Effect, scales="free") +
                                        #    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
                                        #    scale_y_continuous(trans="log") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        xlab("100*q/n") +
        theme_bw() +
        theme(legend.position = "right",
              plot.title = element_text(size = 10),
              axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
              strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
              legend.title=element_text(size = 10))
    pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_2_%s.pdf", prop.treat.plot), width=7, height=4.5, units="in")
}


## Plot 3: confidence bands
library(scales)
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

df.dummy <- tibble(Key="Coverage", Value=0.5, k.rel=90, multivariate=c(FALSE,TRUE)) %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    
df.nominal <- tibble(Key="Coverage", Value=0.9, k.rel=90, multivariate=c(FALSE,TRUE)) %>%
    mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate")))    

for(p.plot in c(20)) {
    for(prop.treat.plot in c(0.24, 0.76)) {
        df <- results %>%
            filter(p==p.plot, j %in% c(2,3), prop.treat==prop.treat.plot, signal_mean==2, reference=="auto") %>%
            mutate(alpha=0.1) %>%
            mutate(bound.value = pmax(-20, pmin(bound.value, 20))) %>%
            group_by(p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, prop.treat, multivariate, reference, j, bound, k.rel, alpha) %>%
            mutate(tau_k=mean(tau_k)) %>%
            group_by(n,p, p_causal, p_causal_covar, signal_mean, signal_std, delta_inter, rho, prop.treat, multivariate, reference, j, bound, k.rel, alpha) %>%
            summarize(Num.exp=n(), Bound.se=2*sd(bound.value)/sqrt(n()), Bound.mean=mean(bound.value), tau_k=mean(tau_k)) %>%
            filter((bound=="lower")&(k.rel==90) | (bound=="upper")&(k.rel==10))

        pp <- df %>%
            mutate(Statistics = factor(multivariate, c(FALSE,TRUE), c("Univariate", "Multivariate"))) %>%
            filter(! ((j == 3) && (bound=="upper"))) %>%
            mutate(Bound = factor(bound, c("lower","upper"), c(TeX("Lower bound on top 10%$"), TeX("Upper bound on bottom 10%$")))) %>%    
            mutate(j = factor(j, c(2,3), c("ILE ~ N(0,0)", "ILE ~ N(2,1)"))) %>%
            ggplot(aes(x=n, y=Bound.mean, color=Statistics, shape=Statistics)) +
            geom_point() +
            geom_line() +
            geom_hline(aes(yintercept=tau_k), linetype=2) +
            facet_wrap(j~Bound, nrow=1) +
            scale_x_continuous(trans='log10') +
            scale_y_continuous(trans="S_sqrt", breaks=c(-20,-10,-2,0,2,10,20)) +
            ylab("Individual linear effects") +
            xlab("Sample size") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            theme_bw() +
            theme(legend.position = "right",
                  plot.title = element_text(size = 10),
                  axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
                  strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
                  legend.title=element_text(size = 10),
                  panel.spacing.y = unit(1, "lines"))
        pp

        pp %>% ggsave(file=sprintf("../figures/experiment_crt_ci_3_bands_p%d_treat%s.pdf", p.plot, prop.treat.plot), width=7, height=2.5, units="in")
    }
}
