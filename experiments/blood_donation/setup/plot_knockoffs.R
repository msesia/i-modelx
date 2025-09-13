#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

seed = 1        ## Random seed (for synthetic data)

library(kableExtra)
suppressMessages(library(tidyverse))
options("width"=200)

ifile <- sprintf("../blood_donation_data/knockoffs/dt_donor_mar042018_simple_N5_seed1.csv")
data <- read_delim(ifile, delim="\t")

X <- data.matrix(data[,1:5])
Xk <- data.matrix(data[,6:10])
X.all <- cbind(X, Xk)


Sigma <- round(cov(scale(X.all)),3)

s.names <- c(sprintf("$X_%d$", 1:5), sprintf("$\\tilde{X}_%d$", 1:5))
rownames(Sigma) <- s.names
colnames(Sigma) <- s.names

for(j in 1:5) {
    Sigma[j,j+5] <- cell_spec(Sigma[j,j+5], "latex", color = "red")
}
for(j in 1:5) {
    for(k in 1:5) {
        Sigma[j+5,k] <- cell_spec(Sigma[j+5,k], "latex", color = "white")
    }
}

tb <- kable(Sigma, "latex", align="r", booktabs=TRUE, escape=FALSE) %>%
    add_header_above(c(" "=1, "Treatments" = 5, "Knockoffs" = 5))
 #   row_spec(1:5, color = "black", background = "#D7261E")
tb
tb %>% save_kable("tables/knockoffs_table.tex")



means <- round(apply(X.all, 2, mean),3)
sds <- round(apply(X.all, 2, sd),3)

df <- as_tibble(t(tibble(Mean=means, `Standard deviation`=sds))) %>%
    mutate(` `=c("Mean", "Standard deviation")) %>%
    select(` `, everything())
colnames(df) <- c(" ", s.names)
tb2 <- kable(df, "latex", align="r", booktabs=TRUE, escape=FALSE) %>%
    add_header_above(c(" "=1, "Treatments" = 5, "Knockoffs" = 5))
tb2
tb2 %>% save_kable("tables/knockoffs_table_2.tex")
