
# This file creates the plots for the simulations on the real UK Biobank genotypes #

#### vary signal amplitude ########## 


rm(list = ls())

for(num_cov in c(1, 2)) {
  
  for(sparsity in c(0.04)) {
    
    
    if(num_cov == 1) {
      color.scale <- c("blue2","cyan2", "#009E73", "#CC79A7", "orange", 
                       "#BA55D3", "grey60", "grey80")
      shape.scale <- c(15, 15, 4, 8, 17, 16, 20, 20)
      linetype.scale <- c(1,1,1,1, 1, 1, 1, 1)
      
      color.scale.no.uwalkf <- c("blue2", "#009E73", "#CC79A7", "orange", 
                                 "#BA55D3", "grey60", "grey80")
      shape.scale.no.uwalkf <- c(15, 4, 8, 17, 16, 20, 20)
      
      color.scale.main <-  c("blue2", "#009E73", "#CC79A7", "orange", "#BA55D3")
      shape.scale.main <- c(15, 4, 8, 17, 16)
    }
    
    if(num_cov == 2) {
      color.scale <- c("blue2","cyan2", "#009E73", "#CC79A7", "orange", 
                       "#BA55D3", "grey60", "grey60", "grey80", "grey80")
      shape.scale <- c(15, 15, 4, 8, 17, 16, 20, 20, 20, 20)
      linetype.scale <- c(1,1,1,1, 1, 1, 1, 1,1, 1)
      
      color.scale.no.uwalkf <- c("blue2", "#009E73", "#CC79A7", "orange", 
                                 "#BA55D3", "grey60", "grey60", "grey80", "grey80")
      shape.scale.no.uwalkf <- c(15, 4, 8, 17, 16, 20, 20, 20, 20)
      
      color.scale.main <-  c("blue2", "#009E73", "#CC79A7", "orange", "#BA55D3")
      shape.scale.main <- c(15, 4, 8, 17, 16)
    }
    
    
    indir <- ""
    
    outdir <- ""
    
    my_population = "whitenonbritish"
    resolution = 1
    
    amp = c(5,10, 15,20, 25,30, 35,40)
    
    fdrs = c(0.1)
    B = 100
    mainprop = c(0.5)
    
    prop_z0 = 1
    prop_z1 = 0
    
    random_sample_prop = 0.3
    
    if(num_cov == 2) {
      method.levels = c( "sskf_lasso_",
                         "sskf_lasso_weight0.25_",
                         "naive_lasso_", 
                         "split_lasso_", 
                         "vanilla_lasso_",
                         "new_combined_separate_lasso_", 
                         "new_z00_separate_lasso_", 
                         "new_z10_separate_lasso_", 
                         "new_z01_separate_lasso_", 
                         "new_z11_separate_lasso_") 
      method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                           paste0("sskf_lasso_weight0.25_mp", mainprop),
                           paste0("naive_lasso_mp", mainprop),  
                           paste0("split_lasso_mp", mainprop), 
                           paste0("vanilla_lasso_mp", mainprop),
                           paste0("new_combined_separate_lasso_mp", mainprop), 
                           paste0("new_z00_separate_lasso_mp", mainprop),
                           paste0("new_z01_separate_lasso_mp", mainprop), 
                           paste0("new_z10_separate_lasso_mp", mainprop),
                           paste0("new_z11_separate_lasso_mp", mainprop))  #,
      
      method.labels <- c("uw-aLKF",
                         "aLKF",
                         "LKF-Naive", 
                         "LKF-Split", 
                         "KF-Global", 
                         "Fixed-LKF", 
                         "Group 1",
                         "Group 2",
                         "Group 3",
                         "Group 4")
      
    }
    
    
    if(num_cov == 1) {
      
      method.levels = c( "sskf_lasso_",
                         "sskf_lasso_weight0.25_"
                         ,"naive_lasso_", 
                         "split_lasso_", 
                         "vanilla_lasso_", 
                         "new_combined_separate_lasso_", 
                         "new_z0_separate_lasso_", 
                         "new_z1_separate_lasso_") 
      method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                           paste0("sskf_lasso_weight0.25_mp", mainprop),
                           paste0("naive_lasso_mp", mainprop),  
                           paste0("split_lasso_mp", mainprop), 
                           paste0("vanilla_lasso_mp", mainprop),
                           paste0("new_combined_separate_lasso_mp", mainprop), 
                           paste0("new_z0_separate_lasso_mp", mainprop), 
                           paste0("new_z1_separate_lasso_mp", mainprop))  
      
      method.labels <- c("uw-aLKF",
                         "aLKF",
                         "LKF-Naive", 
                         "LKF-Split", 
                         "KF-Global", 
                         "Fixed-LKF",
                         "Group 1",
                         "Group 2")
    }
    
    all_res = c()
    
    file_existance = c()
    
    for(p in mainprop) {
      for(m in method.levels) {
        for(a in amp) {
          for(s in sparsity) {
            for(f in fdrs) {
              
              fln = paste0(indir,m ,my_population, "_res", resolution,
                           "_sim_B",B, "_num_cov", num_cov, "_","s",s, "_a", a, "_f", f,"_mp", p,
                           "_z1", prop_z1,"_z0", prop_z0, "_r", random_sample_prop, ".txt")
              
              # check for file existance
              file_existance = rbind(file_existance, data.frame(method = m, amp = a, sparsity = s, fdr = f, mp = p, file_exists = file.exists(fln)))
              
              if(file.exists(fln)) {
                res = read.delim(fln, sep = " ")
                res$method = paste0(m, "mp", p)
                
                if(m == "vanilla_lasso_") {
                  res = res %>% mutate(homogeneity_complement = NA)
                }
                
                all_res = rbind(all_res, res)
              }
              
            }
          }
        }
      }
    }
    
    
    print(file_existance %>% filter(file_exists == FALSE))
    
    all_res_means_orig = all_res %>%
      dplyr::select(-c("propmain")) %>%
      mutate(fdp = ifelse(is.na(fdp), 0, fdp)) %>%
      group_by(sparsity, amp, fdr, method) %>%
      summarise_all(mean, na.rm = TRUE)
    
    all_res_means = all_res_means_orig %>% 
      dplyr::select(c(sparsity, amp, fdr, method, 
                      fdp, power_cond_association, power_cond_association_deconstruct_main,
                      power_cond_association_deconstruct_int, homogeneity, heritability))
    
    measure.labels = c("FDP", "Power", "Global power", "Local power", "Homogeneity", "Heritability")
    measure.levels = c("fdp", "power_cond_association","power_cond_association_deconstruct_main",
                       "power_cond_association_deconstruct_int",
                       "homogeneity", "heritability")

    all_res_means = all_res_means %>%
      pivot_longer(!c(sparsity, amp, fdr, method), names_to = "measure", values_to = "value") %>%
      mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels),
             method = factor(method, levels=method.levels.mp, labels=method.labels))
    
    df.nominal <- tibble(measure=c("fdp"), value=c(0.1), method=c("vanilla_lasso_"))
    df.nominal <- df.nominal %>%
      mutate(measure=factor(measure, levels=measure.levels, labels=measure.labels),
             method=factor(method, levels=method.levels.mp, labels=method.labels))
    
    all_res_means %>%
      filter(!(method %in% c("Group 1", "Group 2", "Group 3", "Group 4", "uw-aLKF"))) %>%
      filter(!(measure %in% c("Heritability"))) %>%
      mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
      ggplot(aes(x = amp, y = value, col = method, shape = method, linetype = method)) +
      geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
      facet_wrap(~ measure, nrow = 1) + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) + # 
      labs(x = "Signal amplitude", y = "", color = "Method", shape = "Method", linetype = "Method") +
      theme(legend.position = "bottom") +
      scale_color_manual(values=color.scale.main) +
      scale_shape_manual(values=shape.scale.main)  +
      scale_linetype_manual(values = linetype.scale) +
      scale_x_continuous(breaks = seq(0, 50, 10)) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
            strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
            legend.title=element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.margin=ggplot2::margin(0,0,0,0),
            plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
      theme(legend.key.size =  unit(0.4, "in")) +
      guides(color = guide_legend(ncol = 5)) 
    
    ggsave(paste0(outdir, my_population, "_res", resolution,
                  "_varyamp_z1", prop_z1,  "_z0", 
                  prop_z0, "_r",
                  random_sample_prop,
                  "_M", mainprop, "_with_prop_", 
                  num_cov,"cov_fdr0.1_rep", B, "_s", sparsity, ".pdf"), width=6.5, height=2.5)
    
    # for appendix
    all_res_means %>%
      mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
      ggplot(aes(x = amp, y = value, col = method, shape = method, linetype = method)) +
      geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
      facet_wrap(~ measure, nrow = 1) + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) + # 
      labs(x = "Signal amplitude", y = "", color = "Method", shape = "Method", linetype = "Method") +
      theme(legend.position = "bottom") +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale)  +
      scale_linetype_manual(values = linetype.scale) +
      scale_x_continuous(breaks = seq(0, 50, 10)) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
            strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
            legend.title=element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.margin=ggplot2::margin(0,0,0,0),
            plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
      theme(legend.key.size =  unit(0.4, "in")) +
      guides(color = guide_legend(ncol = 5)) 
    
    
    ggsave(paste0(outdir, my_population, "_res", resolution,
                  "_varyamp_z1", prop_z1,  "_z0", 
                  prop_z0, "_r",
                  random_sample_prop,
                  "_M", mainprop, "_with_prop_", 
                  num_cov,"cov_fdr0.1_rep", B, "_s", sparsity,"_appendix.pdf"), width=8, height=3)
    
    
  }
  
  
}



######### ___other sparsity levels for appendix #########

rm(list = ls())

for(num_cov in c(1, 2)) {
  
  sparsity = c(0.02, 0.08)
  
  if(num_cov == 1) {
    color.scale <- c("blue2","cyan2", "#009E73", "#CC79A7", "orange", 
                     "#BA55D3", "grey60", "grey80")
    shape.scale <- c(15, 15, 4, 8, 17, 16, 20, 20)
    linetype.scale <- c(1,1,1,1, 1, 1, 1, 1)
    
    color.scale.no.uwalkf <- c("blue2", "#009E73", "#CC79A7", "orange", 
                               "#BA55D3", "grey60", "grey80")
    shape.scale.no.uwalkf <- c(15, 4, 8, 17, 16, 20, 20)
    
    color.scale.main <-  c("blue2", "#009E73", "#CC79A7", "orange", "#BA55D3")
    shape.scale.main <- c(15, 4, 8, 17, 16)
  }
  
  if(num_cov == 2) {
    color.scale <- c("blue2","cyan2", "#009E73", "#CC79A7", "orange", 
                     "#BA55D3", "grey60", "grey60", "grey80", "grey80")
    shape.scale <- c(15, 15, 4, 8, 17, 16, 20, 20, 20, 20)
    linetype.scale <- c(1,1,1,1, 1, 1, 1, 1,1, 1)
    
    color.scale.no.uwalkf <- c("blue2", "#009E73", "#CC79A7", "orange", 
                               "#BA55D3", "grey60", "grey60", "grey80", "grey80")
    shape.scale.no.uwalkf <- c(15, 4, 8, 17, 16, 20, 20, 20, 20)
    
    color.scale.main <-  c("blue2", "#009E73", "#CC79A7", "orange", "#BA55D3")
    shape.scale.main <- c(15, 4, 8, 17, 16)
  }
  
  
  indir <- ""
  
  outdir <- ""
  
  my_population = "whitenonbritish"
  resolution = 1
  
  amp = c(5,10, 15,20, 25,30, 35,40)
  
  fdrs = c(0.1)
  B = 100
  mainprop = c(0.5)
  
  prop_z0 = 1
  prop_z1 = 0
  
  random_sample_prop = 0.3
  
  if(num_cov == 2) {
    method.levels = c( "sskf_lasso_",
                       "sskf_lasso_weight0.25_"
                       ,"naive_lasso_", 
                       "split_lasso_", 
                       "vanilla_lasso_",
                       "new_combined_separate_lasso_", 
                       "new_z00_separate_lasso_", 
                       "new_z10_separate_lasso_", 
                       "new_z01_separate_lasso_", 
                       "new_z11_separate_lasso_") #, "vanillaint_lasso_"
    method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                         paste0("sskf_lasso_weight0.25_mp", mainprop),
                         paste0("naive_lasso_mp", mainprop),  
                         paste0("split_lasso_mp", mainprop), 
                         paste0("vanilla_lasso_mp", mainprop),
                         paste0("new_combined_separate_lasso_mp", mainprop), 
                         paste0("new_z00_separate_lasso_mp", mainprop),
                         paste0("new_z01_separate_lasso_mp", mainprop), 
                         paste0("new_z10_separate_lasso_mp", mainprop),
                         paste0("new_z11_separate_lasso_mp", mainprop))  #,
    
    method.labels <- c("uw-aLKF",
                       "aLKF",
                       "LKF-Naive", 
                       "LKF-Split", 
                       "KF-Global", 
                       "Fixed-LKF", 
                       "Group 1",
                       "Group 2",
                       "Group 3",
                       "Group 4")
    
  }
  
  
  if(num_cov == 1) {
    
    method.levels = c( "sskf_lasso_",
                       "sskf_lasso_weight0.25_"
                       ,"naive_lasso_", 
                       "split_lasso_", 
                       "vanilla_lasso_", 
                       "new_combined_separate_lasso_", 
                       "new_z0_separate_lasso_", 
                       "new_z1_separate_lasso_") #, "vanillaint_lasso_"
    method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                         paste0("sskf_lasso_weight0.25_mp", mainprop),
                         paste0("naive_lasso_mp", mainprop),  
                         paste0("split_lasso_mp", mainprop), 
                         paste0("vanilla_lasso_mp", mainprop),
                         paste0("new_combined_separate_lasso_mp", mainprop), 
                         paste0("new_z0_separate_lasso_mp", mainprop), 
                         paste0("new_z1_separate_lasso_mp", mainprop))  #,
    
    method.labels <- c("uw-aLKF",
                       "aLKF",
                       "LKF-Naive", 
                       "LKF-Split", 
                       "KF-Global", 
                       "Fixed-LKF",
                       "Group 1",
                       "Group 2")
  }
  
  all_res = c()
  
  file_existance = c()
  
  for(p in mainprop) {
    for(m in method.levels) {
      for(a in amp) {
        for(s in sparsity) {
          for(f in fdrs) {
            
            fln = paste0(indir,m ,my_population, "_res", resolution,
                         "_sim_B",B, "_num_cov", num_cov, "_","s",s, "_a", a, "_f", f,"_mp", p,
                         "_z1", prop_z1,"_z0", prop_z0, "_r", random_sample_prop, ".txt")
            file_existance = rbind(file_existance, data.frame(method = m, amp = a, sparsity = s, fdr = f, mp = p, file_exists = file.exists(fln)))
            
            if(file.exists(fln)) {
              res = read.delim(fln, sep = " ")
              res$method = paste0(m, "mp", p)
              if(m == "vanilla_lasso_") {
                res = res %>% mutate(homogeneity_complement = NA)
              }
              
              all_res = rbind(all_res, res)
            }
            
          }
        }
      }
    }
  }
  
  
  print(file_existance %>% filter(file_exists == FALSE))
  
  all_res_means_orig = all_res %>%
    dplyr::select(-c("propmain")) %>%
    mutate(fdp = ifelse(is.na(fdp), 0, fdp)) %>%
    group_by(sparsity, amp, fdr, method) %>%
    summarise_all(mean, na.rm = TRUE)
  
  all_res_means = all_res_means_orig %>% 
    dplyr::select(c(sparsity, amp, fdr, method, 
                    fdp, power_cond_association, power_cond_association_deconstruct_main,
                    power_cond_association_deconstruct_int, homogeneity, heritability))
  
  measure.labels = c("FDP", "Power", "Global power", "Local power", "Homogeneity", "Heritability")
  measure.levels = c("fdp", "power_cond_association","power_cond_association_deconstruct_main",
                     "power_cond_association_deconstruct_int",
                     "homogeneity", "heritability")

  all_res_means = all_res_means %>%
    pivot_longer(!c(sparsity, amp, fdr, method), names_to = "measure", values_to = "value") %>%
    mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels),
           method = factor(method, levels=method.levels.mp, labels=method.labels))
  
  df.nominal <- tibble(measure=c("fdp"), value=c(0.1), method=c("vanilla_lasso_"))
  df.nominal <- df.nominal %>%
    mutate(measure=factor(measure, levels=measure.levels, labels=measure.labels),
           method=factor(method, levels=method.levels.mp, labels=method.labels))
  
  all_res_means %>%
    filter(!(method %in% c("Group 1", "Group 2", "Group 3", "Group 4", "uw-aLKF"))) %>%
    filter(!(measure %in% c("Heritability"))) %>%
    mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
    ggplot(aes(x = amp, y = value, col = method, shape = method, linetype = method)) +
    geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
    facet_grid(sparsity_lab ~measure, switch = "y") + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) + # 
    labs(x = "Signal amplitude", y = "", color = "Method", shape = "Method", linetype = "Method") +
    theme(legend.position = "bottom") +
    scale_color_manual(values=color.scale.main) +
    scale_shape_manual(values=shape.scale.main)  +
    scale_linetype_manual(values = linetype.scale) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.margin=ggplot2::margin(0,0,0,0),
          plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
    theme(legend.key.size =  unit(0.4, "in")) +
    guides(color = guide_legend(ncol = 5)) 
  
  
  ggsave(paste0(outdir, my_population, "_res", resolution,
                "_varyamp_z1", prop_z1,  "_z0", 
                prop_z0, "_r",
                random_sample_prop,
                "_M", mainprop, "_with_prop_", 
                num_cov,"cov_fdr0.1_rep", B, "_s", paste0(sparsity, collapse = "_"), ".pdf"), width=8, height=5)
  
  # for appendix
  all_res_means %>%
    mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
    ggplot(aes(x = amp, y = value, col = method, shape = method, linetype = method)) +
    geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
    facet_grid(sparsity_lab ~measure, switch = "y") + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) + # 
    labs(x = "Signal amplitude", y = "", color = "Method", shape = "Method", linetype = "Method") +
    theme(legend.position = "bottom") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale)  +
    scale_linetype_manual(values = linetype.scale) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.margin=ggplot2::margin(0,0,0,0),
          plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
    theme(legend.key.size =  unit(0.4, "in")) +
    guides(color = guide_legend(ncol = 5)) 
  
  
  ggsave(paste0(outdir, my_population, "_res", resolution,
                "_varyamp_z1", prop_z1,  "_z0", 
                prop_z0, "_r",
                random_sample_prop,
                "_M", mainprop, "_with_prop_", 
                num_cov,"cov_fdr0.1_rep", B, "_s", paste0(sparsity, collapse = "_"), "_appendix.pdf"), width=8, height=5)
  
  
  
}

##########  vary proportion global effects #######

rm(list = ls())

for(num_cov in c(1, 2)) {
  
  for(sparsity in c(0.04)) {
    
    
    if(num_cov == 1) {
      color.scale <- c("blue2","cyan2",  "orange", "#BA55D3", "grey60",  "grey80")
      shape.scale <- c(15, 15,  17, 16, 20, 20)
      linetype.scale <- c(1,1,1,1, 1, 1)
    }
    
    if(num_cov == 2) {
      color.scale <- c("blue2","cyan2",  "orange", "#BA55D3", "grey60", "grey60", "grey80", "grey80")
      shape.scale <- c(15, 15,  17, 16, 20, 20, 20, 20)
      linetype.scale <- c(1,1,1,1, 1, 1,1, 1)
    }
    
    
    indir <- ""
    
    outdir <- ""
    
    my_population = "whitenonbritish"
    resolution = 1
    
    amp = c(15) 
    
    mainprop = c(0,0.1,0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1) # 
    
    
    if(num_cov == 1) {
      
      method.levels = c( "sskf_lasso_", 
                         "sskf_lasso_weight0.25_",
                         "vanilla_lasso_",
                         "new_combined_separate_lasso_", 
                         "new_z0_separate_lasso_", 
                         "new_z1_separate_lasso_") 
      
      method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                           paste0("sskf_lasso_weight0.25_mp", mainprop), 
                           paste0("vanilla_lasso_mp", mainprop), 
                           paste0("new_combined_separate_lasso_mp", mainprop), 
                           paste0("new_z0_separate_lasso_mp", mainprop), 
                           paste0("new_z1_separate_lasso_mp", mainprop))
      
    }
    
    
    
    if(num_cov == 2) {
      
      method.levels = c( "sskf_lasso_", 
                         "sskf_lasso_weight0.25_",
                         "vanilla_lasso_",
                         "new_combined_separate_lasso_", 
                         "new_z00_separate_lasso_", 
                         "new_z10_separate_lasso_", 
                         "new_z01_separate_lasso_", 
                         "new_z11_separate_lasso_") 
      
      method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                           paste0("sskf_lasso_weight0.25_mp", mainprop), 
                           paste0("vanilla_lasso_mp", mainprop), 
                           paste0("new_combined_separate_lasso_mp", mainprop), 
                           paste0("new_z00_separate_lasso_mp", mainprop),
                           paste0("new_z01_separate_lasso_mp", mainprop), 
                           paste0("new_z10_separate_lasso_mp", mainprop),
                           paste0("new_z11_separate_lasso_mp", mainprop))
      
    }
    
    fdrs = c(0.1)
    B = 100
    
    prop_z0 = 1
    prop_z1 = 0
    
    random_sample_prop = 0.3
    
    all_res = c()
    
    file_existance = c()
    
    for(p in mainprop) {
      for(m in method.levels) {
        for(a in amp) {
          for(s in sparsity) {
            for(f in fdrs) {
              
              fln = paste0(indir,m ,my_population, "_res", resolution, "_sim_B",B, 
                           "_num_cov",num_cov,"_s",s, "_a", a, "_f", f,"_mp", p,
                           "_z1", prop_z1,"_z0", prop_z0,  "_r",random_sample_prop,".txt")
              file_existance = rbind(file_existance, data.frame(method = m, amp = a, sparsity = s, fdr = f, mp = p, file_exists = file.exists(fln)))
              
              file_existance = rbind(file_existance, data.frame(method = m, amp = a, sparsity = s, fdr = f, mp = p, file_exists = file.exists(fln)))
              if(file.exists(fln)) {
                res = read.delim(fln, sep = " ")
                res$method = paste0(m, "mp", p)
                
                
                
                if(m == "vanilla_lasso_") {
                  res = res %>% mutate(homogeneity_complement = NA)
                }
                
                
                
                all_res = rbind(all_res, res)
              }
              
            }
          }
        }
      }
    }
    
    
    
    print(file_existance %>% filter(file_exists == FALSE))
    
    
    
    all_res_means_orig = all_res %>% 
      dplyr::select(-c(amp)) %>%
      mutate(fdp = ifelse(is.na(fdp), 0, fdp), 
             fdp_interactions = ifelse(is.na(fdp_interactions), 0, fdp_interactions), 
             fdp_main = ifelse(is.na(fdp_main), 0, fdp_main)) %>%
      group_by(propmain, sparsity, fdr, method) %>% 
      summarise_all(mean, na.rm = TRUE)
    
    
    # drop not needed columns 
    all_res_means = all_res_means_orig %>% 
      dplyr::select(-c("homogeneity_complement"))
    
    # FDP MEASURES 
    all_res_means = all_res_means %>% 
      dplyr::select(c(sparsity, propmain, fdr, method, 
                      fdp, power_cond_association, 
                      power_cond_association_deconstruct_main, 
                      power_cond_association_deconstruct_int, 
                      homogeneity, heritability))
    
    
    measure.labels = c("FDP", "Power", "Global power", "Local power", "Homogeneity", "Heritability")
    measure.levels = c("fdp", "power_cond_association","power_cond_association_deconstruct_main", 
                       "power_cond_association_deconstruct_int", 
                       "homogeneity", "heritability")
    
    
    if(num_cov == 2) {
      all_res_means = all_res_means %>% 
        pivot_longer(!c(sparsity, propmain, fdr, method), names_to = "measure", values_to = "value") %>% 
        mutate(method_orig = method, 
               method_mod = ifelse(grepl("sskf_lasso_weight0.25", method), "aLKF", 
                                   ifelse(grepl("vanilla_lasso", method), "KF-Global", 
                                          ifelse(grepl("sskf_lasso_", method), "uw-aLKF", 
                                                 ifelse(grepl("combined_separate_", method), "Fixed-LKF", 
                                                        ifelse(grepl("z00_separate", method), "Group 1", 
                                                               ifelse(grepl("z01_separate", method), "Group 2", 
                                                                      ifelse(grepl("z10_separate", method), "Group 3", 
                                                                             ifelse(grepl("z11_separate", method), "Group 4", NA))))))))) %>%
        mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels), 
               method = factor(method_mod, levels = c("uw-aLKF", "aLKF", "KF-Global", "Fixed-LKF", "Group 1", "Group 2", "Group 3", "Group 4")))
      
    }
    
    
    if(num_cov == 1) {
      all_res_means = all_res_means %>% 
        pivot_longer(!c(sparsity, propmain, fdr, method), names_to = "measure", values_to = "value") %>% 
        mutate(method_orig = method, 
               method_mod = ifelse(grepl("sskf_lasso_weight0.25", method), "aLKF", 
                                   ifelse(grepl("vanilla_lasso", method), "KF-Global", 
                                          ifelse(grepl("sskf_lasso_", method), "uw-aLKF", 
                                                 ifelse(grepl("combined_separate_", method), "Fixed-LKF", 
                                                        ifelse(grepl("z0_separate", method), "Group 1", 
                                                               ifelse(grepl("z1_separate", method), "Group 2", NA))))))) %>%
        mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels), 
               method = factor(method_mod, levels = c("uw-aLKF", "aLKF", "KF-Global", "Fixed-LKF", "Group 1", "Group 2")))
      
    }
    
    df.nominal <- tibble(measure=c("fdp"), value=c(0.1), method=c("KF-Global"))
    df.nominal <- df.nominal %>%
      mutate(measure=factor(measure, levels=measure.levels, labels=measure.labels),
             method=factor(method))
    
    all_res_means %>% 
      mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
      ggplot(aes(x = propmain, y = value, col = method, shape = method, linetype = method)) +
      geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
      facet_wrap(~ measure, nrow = 1) + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) +  
      labs(x = "Proportion of global effects", y = "", color = "Method", shape = "Method", linetype = "Method") + 
      theme(legend.position = "bottom") +  
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale)  +
      scale_linetype_manual(values = linetype.scale) + 
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
            strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
            legend.title=element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.margin=ggplot2::margin(0,0,0,0),
            plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
      scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
      theme(legend.key.size =  unit(0.4, "in")) +
      guides(color = guide_legend(ncol = 4)) 
    
    ggsave(paste0(outdir, my_population, "_res", resolution, "_varyprop_z1", prop_z1,
                  "_z0", prop_z0, "_r", random_sample_prop, 
                  "_amp", amp, "_fdr0.1_rep", B,"_num_cov",num_cov, "_s", sparsity ,".pdf"), width = 8.5, height = 2.7)
    
  }
  
  
}


######### ___other sparsity levels for appendix ############ 


rm(list = ls())

for(num_cov in c(1, 2)) {
  
  if(num_cov == 1) {
    color.scale <- c("blue2","cyan2",  "orange", "#BA55D3", "grey60",  "grey80")
    shape.scale <- c(15, 15,  17, 16, 20, 20)
    linetype.scale <- c(1,1,1,1, 1, 1)
  }
  
  if(num_cov == 2) {
    color.scale <- c("blue2","cyan2",  "orange", "#BA55D3", "grey60", "grey60", "grey80", "grey80")
    shape.scale <- c(15, 15,  17, 16, 20, 20, 20, 20)
    linetype.scale <- c(1,1,1,1, 1, 1,1, 1)
  }
  
  indir <- ""
  
  outdir <- ""
  
  my_population = "whitenonbritish"
  resolution = 1
  
  amp = c(15) 
  
  sparsity = c(0.02, 0.08)
  
  mainprop = c(0,0.1,0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1) # 
  
  
  if(num_cov == 1) {
    
    method.levels = c( "sskf_lasso_", 
                       "sskf_lasso_weight0.25_",
                       "vanilla_lasso_",
                       "new_combined_separate_lasso_", 
                       "new_z0_separate_lasso_", 
                       "new_z1_separate_lasso_") 
    
    method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                         paste0("sskf_lasso_weight0.25_mp", mainprop), 
                         paste0("vanilla_lasso_mp", mainprop), 
                         paste0("new_combined_separate_lasso_mp", mainprop), 
                         paste0("new_z0_separate_lasso_mp", mainprop), 
                         paste0("new_z1_separate_lasso_mp", mainprop))
    
  }
  
  
  
  if(num_cov == 2) {
    
    method.levels = c( "sskf_lasso_", 
                       "sskf_lasso_weight0.25_",
                       "vanilla_lasso_",
                       "new_combined_separate_lasso_", 
                       "new_z00_separate_lasso_", 
                       "new_z10_separate_lasso_", 
                       "new_z01_separate_lasso_", 
                       "new_z11_separate_lasso_") 
    
    
    
    method.levels.mp = c(paste0("sskf_lasso_mp", mainprop),
                         paste0("sskf_lasso_weight0.25_mp", mainprop), 
                         paste0("vanilla_lasso_mp", mainprop), 
                         paste0("new_combined_separate_lasso_mp", mainprop), 
                         paste0("new_z00_separate_lasso_mp", mainprop),
                         paste0("new_z01_separate_lasso_mp", mainprop), 
                         paste0("new_z10_separate_lasso_mp", mainprop),
                         paste0("new_z11_separate_lasso_mp", mainprop))
    
  }
  
  fdrs = c(0.1)
  B = 100
  
  prop_z0 = 1
  prop_z1 = 0
  
  random_sample_prop = 0.3
  
  all_res = c()
  
  file_existance = c()
  
  for(p in mainprop) {
    for(m in method.levels) {
      for(a in amp) {
        for(s in sparsity) {
          for(f in fdrs) {
            
            fln = paste0(indir,m ,my_population, "_res", resolution, "_sim_B",B, 
                         "_num_cov",num_cov,"_s",s, "_a", a, "_f", f,"_mp", p,
                         "_z1", prop_z1,"_z0", prop_z0,  "_r",random_sample_prop,".txt")
            file_existance = rbind(file_existance, data.frame(method = m, amp = a, sparsity = s, num_cov = num_cov,
                                                              fdr = f, mp = p, file_exists = file.exists(fln)))
            
            if(file.exists(fln)) {
              res = read.delim(fln, sep = " ")
              res$method = paste0(m, "mp", p)
              
              
              
              if(m == "vanilla_lasso_") {
                res = res %>% mutate(homogeneity_complement = NA)
              }
              
              
              
              all_res = rbind(all_res, res)
            }
            
          }
        }
      }
    }
  }
  
  
  
  print(file_existance %>% filter(file_exists == FALSE))
  
  
  
  all_res_means_orig = all_res %>% 
    dplyr::select(-c(amp)) %>%
    mutate(fdp = ifelse(is.na(fdp), 0, fdp), 
           fdp_interactions = ifelse(is.na(fdp_interactions), 0, fdp_interactions), 
           fdp_main = ifelse(is.na(fdp_main), 0, fdp_main)) %>%
    group_by(propmain, sparsity, fdr, method) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  
  # drop not needed columns 
  all_res_means = all_res_means_orig %>% 
    dplyr::select(-c("homogeneity_complement"))
  
  # FDP MEASURES 
  all_res_means = all_res_means %>% 
    dplyr::select(c(sparsity, propmain, fdr, method, 
                    fdp, power_cond_association, 
                    power_cond_association_deconstruct_main, 
                    power_cond_association_deconstruct_int, 
                    homogeneity, heritability))
  
  
  measure.labels = c("FDP", "Power", "Global power", "Local power", "Homogeneity", "Heritability")
  measure.levels = c("fdp", "power_cond_association","power_cond_association_deconstruct_main", 
                     "power_cond_association_deconstruct_int", 
                     "homogeneity", "heritability")
  
  
  if(num_cov == 2) {
    all_res_means = all_res_means %>% 
      pivot_longer(!c(sparsity, propmain, fdr, method), names_to = "measure", values_to = "value") %>% 
      mutate(method_orig = method, 
             method_mod = ifelse(grepl("sskf_lasso_weight0.25", method), "aLKF", 
                                 ifelse(grepl("vanilla_lasso", method), "KF-Global", 
                                        ifelse(grepl("sskf_lasso_", method), "uw-aLKF", 
                                               ifelse(grepl("combined_separate_", method), "Fixed-LKF", 
                                                      ifelse(grepl("z00_separate", method), "Group 1", 
                                                             ifelse(grepl("z01_separate", method), "Group 2", 
                                                                    ifelse(grepl("z10_separate", method), "Group 3", 
                                                                           ifelse(grepl("z11_separate", method), "Group 4", NA))))))))) %>%
      mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels), 
             method = factor(method_mod, levels = c("uw-aLKF", "aLKF", "KF-Global", "Fixed-LKF", "Group 1", "Group 2", "Group 3", "Group 4")))
    
  }
  
  
  if(num_cov == 1) {
    all_res_means = all_res_means %>% 
      pivot_longer(!c(sparsity, propmain, fdr, method), names_to = "measure", values_to = "value") %>% 
      mutate(method_orig = method, 
             method_mod = ifelse(grepl("sskf_lasso_weight0.25", method), "aLKF", 
                                 ifelse(grepl("vanilla_lasso", method), "KF-Global", 
                                        ifelse(grepl("sskf_lasso_", method), "uw-aLKF", 
                                               ifelse(grepl("combined_separate_", method), "Fixed-LKF", 
                                                      ifelse(grepl("z0_separate", method), "Group 1", 
                                                             ifelse(grepl("z1_separate", method), "Group 2", NA))))))) %>%
      mutate(measure = factor(measure, levels = measure.levels, labels = measure.labels), 
             method = factor(method_mod, levels = c("uw-aLKF", "aLKF", "KF-Global", "Fixed-LKF", "Group 1", "Group 2")))
    
  }
  
  df.nominal <- tibble(measure=c("fdp"), value=c(0.1), method=c("KF-Global"))
  df.nominal <- df.nominal %>%
    mutate(measure=factor(measure, levels=measure.levels, labels=measure.labels),
           method=factor(method))
  
  all_res_means %>% 
    mutate(sparsity_lab = paste0("sparsity snp = ", round(sparsity/4, 3))) %>%
    ggplot(aes(x = propmain, y = value, col = method, shape = method, linetype = method)) +
    geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
    #facet_wrap(~ measure, nrow = 1) + geom_line() + geom_point() +  
    facet_grid(sparsity_lab ~ measure, switch = "y") + geom_line(linewidth=0.7, alpha = 0.6) + geom_point(size=2, alpha = 0.6) +
    labs(x = "Proportion of global effects", y = "", color = "Method", shape = "Method", linetype = "Method") + 
    theme(legend.position = "bottom") +  
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale)  +
    scale_linetype_manual(values = linetype.scale) + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.margin=ggplot2::margin(0,0,0,0),
          plot.margin = ggplot2::margin(0,0.3,0,0, "cm")) +
    theme(legend.key.size =  unit(0.4, "in")) + #, panel.spacing = unit(0.13, "in"), 
    guides(color = guide_legend(ncol = 4)) + 
    scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))
  
  ggsave(paste0(outdir, my_population, "_res", resolution, "_varyprop_z1", prop_z1,
                "_z0", prop_z0, "_r", random_sample_prop, 
                "_amp", amp, "_fdr0.1_rep", B,"_num_cov",num_cov, "_s", paste0(sparsity, collapse = "_") ,".pdf"), width = 8.5, height = 4)
  
}
