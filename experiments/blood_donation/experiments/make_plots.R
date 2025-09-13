options(width=160)

library(tidyverse)

idir <- "results/"
ifile.list <- list.files(idir)

summary <- do.call("rbind", lapply(ifile.list, function(ifile) {
    if(startsWith(ifile, "summary_")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
    } else {
        tibble()
    }
}))

method.levels <- c("skf.cp", "naive", "split", "vanilla")
method.labels <- c("SKF", "Naive", "Split", "Vanilla")

keys.levels <- c("FDR", "Power", "Causal.prop", "Beta.sd")
keys.labels <- c("FDR", "Power", "Homeg. (non-null prop.)", "Heterog. (effect s.d.)")

color.scale <- c("blue2", "#009E73", "#CC79A7", "orange")
shape.scale <- c(15, 4, 8, 1)

df <- summary %>%
    mutate(FDR=FDP) %>% select(-FDP, -True) %>%
    gather(Power, FDR, Causal.prop, Beta.sd, key="Key", value="Value") %>%
    group_by(Method, a, n, num.int, Offset, Key) %>%
    summarise(Value=mean(Value,na.rm=T))

df.ghost.1 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=0, n=1000, Method="vanilla")
df.ghost.2 <- tibble(Key=c("Power", "FDR","Causal.prop"), Value=0.35, n=1000, Method="vanilla")
df.ghost <- rbind(df.ghost.1, df.ghost.2)
df.ghost <- df.ghost %>%
    mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
           Method=factor(Method, levels=method.levels, labels=method.labels))

df.nominal <- tibble(Key="FDR", Value=0.1, Method="Vanilla")
df.nominal <- df.nominal %>%
     mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
            Method=factor(Method, levels=method.levels, labels=method.labels))

pp <- df %>%
    filter(a==0.4, Offset==0, Key!="True", Method!="skf") %>%
    mutate(Key=factor(Key, levels=keys.levels, labels=keys.labels),
           Method=factor(Method, levels=method.levels, labels=method.labels)) %>%
    ggplot(aes(x=n, y=Value, color=Method, shape=Method)) +
    geom_point() +
    geom_line() +
    geom_point(alpha=0, data=df.ghost, aes(x=n,y=Value)) +
    facet_wrap(.~Key, scales="free_y", nrow=1) +
    scale_x_continuous(breaks=c(100, 1000, 10000, 80000), trans='log10') +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          legend.margin=margin(0,0,0,0),
          plot.margin = margin(0,0.3,0,0, "cm"))
pp %>% ggsave(file="figures/experiment_blood.pdf", width=8, height=2, units="in")


##########################
## Confidence intervals ##
##########################

options(width=160)

library(tidyverse)

idir <- "results2/"
ifile.list <- list.files(idir)

results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    if(startsWith(ifile, "n")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
    } else {
        tibble()
    }
}))

df <- results %>%
    filter(reference=="auto") %>%
    mutate(Method = sprintf("%s_multi%s_int%s_subset%s", model.class, multivariate, interactions, subset)) %>%
    select(-model.class, -multivariate, -interactions, -subset) %>%
    mutate(bound.value = pmax(-5, bound.value), bound.value=pmin(5, bound.value),
           Cover = ifelse(bound=="lower", bound.value<=tau_k, bound.value>=tau_k), Distance=abs(bound.value-tau_k))

df.1 <- df %>%
    filter(ci.type=="quantile") %>%
    group_by(a, n, treatment, bound, k.rel, tau_k, ci.type, alpha, Method) %>%
    summarise(Num=n(), Value=mean(bound.value), Coverage = mean(Cover), Distance=mean(Distance)) %>%    
    arrange(ci.type, n, bound, k.rel, Method)

df.1 %>%
    print(n=100)


#################
## Single plot ##
#################

color.scale <- c("tan2", "coral4", "indianred")
shape.scale <- c(1, 0, 8)

method.values <- c("glmnet_multiFALSE_intFALSE_subsetFALSE", "glmnet_multiTRUE_intFALSE_subsetFALSE", "glmnet_multiTRUE_intTRUE_subsetFALSE")
method.labels <- c("Univariate", "Multivariate", "Multivariate w. interactions")

plot.k <- 90
plot.j <- 1
plot.bound <- "lower"
df.nominal <- tibble(Key="Coverage", Value=0.9) %>% mutate(Key = factor(Key, c("Coverage", "Distance", "Bound")))
df.nominal.2 <- tibble(Key="Bound", Value=ifelse(plot.j==5, 1, 2)) %>% mutate(Key = factor(Key, c("Coverage", "Distance", "Bound")))
df.ghost <- tibble(Key="Coverage", Value=0.5) %>% mutate(Key = factor(Key, c("Coverage", "Distance", "Bound")))
df.ghost.2 <- tibble(Key="Bound", Value=-2) %>% mutate(Key = factor(Key, c("Coverage", "Distance", "Bound")))
pp <- df.1 %>%
    filter(treatment==plot.j, bound==plot.bound, k.rel==plot.k) %>%
    mutate(Bound=Value) %>%
    mutate(Statistics = factor(Method, method.values, method.labels)) %>%
    gather(Coverage, Distance, Bound, key="Key", value="Value") %>%
    mutate(Key = factor(Key, c("Coverage", "Distance", "Bound"))) %>%
    ggplot(aes(x=n, y=Value, color=Statistics, shape=Statistics)) +
    geom_point() +
    geom_line() +
    facet_wrap(.~Key, scale="free") +
    geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
    geom_hline(data=df.nominal.2, aes(yintercept=Value), linetype=3) +
    geom_hline(data=df.ghost, aes(yintercept=Value), alpha=0) +
    geom_hline(data=df.ghost.2, aes(yintercept=Value), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_x_continuous(trans='log10', breaks=c(100,1000,10000,80000)) +
    xlab("Sample size") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right", 
          plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          legend.margin=margin(0,0,0,0),
          plot.margin = margin(0,0.3,0,0, "cm"))
pp %>% ggsave(file=sprintf("figures/experiment_blood_continuous_%d_%s_k%d.pdf", plot.j, plot.bound, plot.k), width=8, height=2, units="in")


##########################
## Table of discoveries ##
##########################

options(width=160)

library(tidyverse)
library(kableExtra)

ifile <- "results/discoveries_synthetic_a0.4_n80000_seed100_skf.cp.csv"
discoveries <- read_delim(ifile, delim=" ")
df <- discoveries %>% filter(Offset==0) %>%
    arrange(Treatment, desc(W)) %>%
    mutate(W = round(100*W,2), Causal.prop=round(Causal.prop,2), Beta.sd=round(Beta.sd,2),
           n.sub = sprintf("%s",formatC(n.sub, format="d", big.mark=","))) %>%
    mutate(Label = str_replace_all(Label, "student", "Student"),
           Label = str_replace_all(Label, "married", "Married"),
           Label = str_replace_all(Label, "resident", "Resident"),
           Label = str_replace_all(Label, "male", "Male"),
           Label = str_replace_all(Label, "rh_neg", "Rh-negative"),          
           Label = str_replace_all(Label, ";", " and"),          
           ) %>%
    mutate(`Sub-population`=Label,
           Truth = ifelse(Null, "Null", "Non-null"),
           `Homeg.` = Causal.prop,
           `Heterog.` = Beta.sd,
           Samples = n.sub,
           Treatment = factor(Treatment, 1:5, c("Reminder", "Individual reward", "Friends request", "Group reward", "Small group gift"))) %>%
    select(Treatment, `Sub-population`, Truth, W, `Homeg.`, `Heterog.`, Samples)
df    

df %>%
    kbl("latex", booktabs = T) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
