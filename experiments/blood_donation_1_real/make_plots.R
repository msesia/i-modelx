options(width=160)

library(tidyverse)

idir <- "results/"
ifile.list <- list.files(idir)

discoveries <- do.call("rbind", lapply(ifile.list, function(ifile) {
    if(startsWith(ifile, "discoveries_")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
    } else {
        tibble()
    }
}))

df <- discoveries %>% arrange(num.int, seed, Method)

df %>%
    filter(Method=="skf.cp", Offset==0, num.int==1) %>%
    group_by(Treatment, Label, Variables) %>% summarise(Num=n()) %>%
    arrange(Treatment, desc(Num)) %>%
    filter(Num>=10) %>%
    print(n=100)


df %>%
    filter(Offset==0) %>%
    group_by(num.int, Method, seed) %>%
    summarise(Discoveries=n()) %>%
    ggplot(aes(x=Discoveries)) +
    geom_histogram() +
    facet_grid(num.int~Method) +
    theme_bw()


#########
## CRT ##
#########

idir <- "results_crt/"
ifile.list <- list.files(idir)

results.crt <- do.call("rbind", lapply(ifile.list, function(ifile) {
    if(startsWith(ifile, "pvalues_")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())
    } else {
        tibble()
    }
}))

results.crt$pval.adj <- p.adjust(results.crt$pval, method="BH")

results.crt %>% arrange(pval)

results.crt
