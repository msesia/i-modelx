suppressMessages(library(cowplot))
suppressMessages(library(latex2exp))

rm(list = ls())

# get statistics
Knockoffs = read.csv("lasso_statistics_aLKF.csv")

font.size <- 20
dot.size <- 2.5

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

pre_process <- function(gwas) {
  gwas.don <- gwas %>%
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwas, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  
  return(gwas.don)
}


manhattan_knock <- function(gwas.don, axisdf, limit.left, limit.right, yintercept=0, use.Mb=FALSE) {
  
  # Rescale things to Mb
  limit.left <- limit.left/1e6
  limit.right <- limit.right/1e6
  axisdf$center <- axisdf$center/1e6
  gwas.don$BP <- gwas.don$BP/1e6
  gwas.don$BPcum <- gwas.don$BPcum/1e6
  
  # Make LMM Manhattan plot
  y.max <- 1.1 * max(gwas.don$W)
  
  color_values <- c(
    setNames(rep(c("lightgrey", "darkgrey"), 22), as.character(1:22)),
    "female" = "#FF007F", 
    "male" = "#008CFF"  
  )
  
  
  p.manhattan <- gwas.don %>%     
    mutate(color_group = factor(ifelse(Interaction == "sex" & environment == 1 & W >= W.thresh, "female",
                                       ifelse(Interaction == "sex" & environment == 2 & W >= W.thresh, "male",
                                              as.factor(CHR))))) %>%
    
    ggplot(aes(x=BPcum, y=W)) +
    
    # Show all points
    geom_point( aes(color=color_group), alpha=0.8, size=dot.size) +
    scale_color_manual(values = color_values) +
    
    # Show significance threshold
    geom_hline(yintercept=yintercept, linetype="dashed", color = "red") +
    
    # Custom axes:
    scale_x_continuous(label=axisdf$CHR, breaks=axisdf$center, limits=c(limit.left,limit.right),
                       expand=c(0.01,0.01)) +
    scale_y_continuous(name="Test statistics", labels=function(x) round(x,3)) +
    xlab("Chromosome") + ylab("Importance") +
    
    # Custom the theme:
    theme_bw() +
    theme(legend.position="none", panel.border = element_blank(), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
          panel.grid.major = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor = element_line(size = 0.1, colour = "lightgray"),
          text = element_text(size=font.size)
    )
  
  return(p.manhattan)
}

# Transform knockoffs results into pseudo p-values
W.thresh <- knockoff.threshold(Knockoffs$W)
Knockoffs <- Knockoffs %>% 
  filter(W>0) %>% 
  mutate(SNP=SNP.lead, BP=BP.lead) %>%
  select(CHR, SNP, BP, W, Interaction, environment)


Knockoffs.don <- pre_process(Knockoffs)
don <- Knockoffs.don
# Compute plot limits
limit.left <- min(don$BPcum)
limit.right <- max(don$BPcum)
# Find centers of each chromosome
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make Knockoffs Manhattan plot
ytrans="log10"
p.knockoffs <- manhattan_knock(Knockoffs.don, axisdf, limit.left, limit.right, 
                               yintercept=W.thresh)


p.knockoffs
ggsave(paste0("manhattan_plot" ,".pdf"), width = 14, height = 6)


