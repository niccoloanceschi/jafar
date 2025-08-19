
rm(list = ls())

library(fields)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multiview)
library(knitr)
library(kableExtra)

if(T){
  mth_names = c('BIP',   'CoopL',  'IntegL', 'BSFP',   'JAFAR', 'JAFAR_T',  'JFR',    'JFR_T')
  v_colors = c('#8D2FE5','#DC267F','#FE6100','#FFB000','#648FFF','#2DC4DC','#2B5B43','#57B38C')
  v_pch <- c(8, 22, 23, 21, 24, 25, 15, 17)
  v_pch <- c(8, 24, 25, 21, 22, 27, 23, 9)
}

v_meth = c('JAFAR','JFR','BSFP','BIP','IntegL','CoopL')

mth_ci = which(v_meth %in% c('JAFAR','JFR','BSFP','IntegL'))

idx_match = match(v_meth,mth_names)
mth_names = mth_names[idx_match]
v_colors = v_colors[idx_match]
v_pch = v_pch[idx_match]

ris_sim <- readRDS('~/Desktop/jafar_sec3_sim_vs_n/Sec3_Simulations_results.rds')
list2env(ris_sim, envir = .GlobalEnv)  # puts a, b, c into the global environment
# rm(ris_sim)

n_values = unique(time_run[,'n'])
s_values = unique(time_run[,'r'])

# Response predictions -----

which_data <- list(train=mse_train,test=mse_test)

col_targets <- c('mse_y','R_squared','coverage90s')
my_labs <- c('Response MSE', 'R-Squared', '90% CI Empirical Coverage')
file_name <- c('mse', 'Rsquared','90ci_cov')

for(cc in 1:length(col_targets)){
  for(dd in 1:length(which_data)){
    result <- as.data.frame(as.matrix(as.data.frame(which_data[[dd]]) %>%
                                        group_by(n, method) %>%
                                        summarize(
                                          Median = median(!!sym(col_targets[cc]),na.rm=T),
                                          Quantile1 = quantile(!!sym(col_targets[cc]), 0.25,na.rm=T),
                                          Quantile2 = quantile(!!sym(col_targets[cc]), 0.75,na.rm=T)
                                        )))
    
    sel_cols = v_colors
    sel_meth = v_meth
    sel_pch = v_pch
    if(col_targets[cc]=='coverage90s'){
      result <- result[result[,'method'] %in% mth_ci,]
      sel_cols = v_colors[mth_ci]
      sel_meth = v_meth[mth_ci]
      sel_pch = v_pch[mth_ci]
    }
    
    combined_plot <- ggplot() + theme_bw() +
      geom_ribbon(data = result, aes(x=n, y = Median, ymin = Quantile1, ymax = Quantile2, fill=factor(method), color=factor(method)),
                  show.legend = FALSE,alpha=0.05, linewidth=0.5, linetype = "dotted") +
      geom_point(data = result, aes(x = n, y = Median, color=factor(method),
                                    shape=factor(method), fill=factor(method)),size=1.6) +
      theme(aspect.ratio=0.5, legend.position="bottom") +
      labs(x = "n (train)", y = my_labs[cc]) +
      scale_x_continuous(breaks=n_values) +
      scale_colour_manual(name="Methods:",values=sel_cols, labels=sel_meth) +
      scale_shape_manual(name="Methods:",values=sel_pch, labels=sel_meth) +
      scale_fill_manual(name="Methods:",values=sel_cols, labels=sel_meth) +
      guides(colour = guide_legend(override.aes = list(size=1))) + 
      guides(colour = guide_legend(override.aes = list(size=2)))
    if(col_targets[cc]=='mse_y'){combined_plot <- combined_plot + scale_y_log10()}
    if(col_targets[cc]=='coverage90s'){
      combined_plot <- combined_plot + scale_y_continuous(limits = c(0.5,1), breaks=seq(0,1,by=0.1)) +
        # ylim(0, 1) + scale_y_continuous(breaks=seq(0,1,by=0.1)) +
        geom_hline(yintercept=0.90, color = "black",linewidth =0.3,lty=1)
    }
    
    ggsave(file.path('~/Desktop/',paste0('sim_',names(which_data)[dd],'_',file_name[cc],'.pdf')), combined_plot, height=3, width=4.7)
    
    show(combined_plot)
  }
}


# cor frob norm -----

mth_cor = which(v_meth %in% c('JAFAR','JFR','BSFP','BIP'))

sel_cols = v_colors[mth_cor]
sel_meth = v_meth[mth_cor]
sel_pch = v_pch[mth_cor]

result <- cor_fn2[cor_fn2[,'exact']==0,!colnames(cor_fn2) %in% "exact"]

result <- as.data.frame(as.matrix(as.data.frame(result) %>%
                                    group_by(n, method, m) %>%
                                    summarize(
                                      Median = median(cor_fn2,na.rm=T),
                                      Quantile1 = quantile(cor_fn2, 0.25,na.rm=T),
                                      Quantile2 = quantile(cor_fn2, 0.75,na.rm=T)
                                    )))

custom_labeller <- function(value){
  m = as.numeric(strsplit(as.character(value), "")[[1]])
  if(length(m)==1){
    lab_cor = paste0("$cor(X_{",m[1],"})$")
  } else {
    lab_cor = paste0("$cor(X_{",m[1],"}, X_{",m[2],"})$")
  }
  return(lab_cor)
}

m_values = unique(result$m)
m_labs = sapply(m_values, function(vv) custom_labeller(vv))

result$m = factor(result$m, levels=m_values)
result$m = plyr::mapvalues(result$m,
                           from = m_values, 
                           to = latex2exp::TeX(m_labs))

combined_plot <- ggplot() + theme_bw() +
  geom_ribbon(data = result, aes(x=n, y = Median, ymin = Quantile1, ymax = Quantile2, fill=factor(method), color=factor(method)),
              show.legend = FALSE,alpha=0.1, linewidth=0.4, linetype = "dotted") +
  geom_point(data = result, aes(x = n, y = Median, color=factor(method),
                                shape=factor(method), fill=factor(method)),size=1.4) +
  facet_wrap(~ m, ncol= 3, strip.position="top", labeller = label_parsed, scales = "free") +
  theme(aspect.ratio=0.4, legend.position="bottom") +
  labs(x = "n (train)", y = "Frobenius Norm") +
  scale_x_continuous(breaks=n_values) +
  scale_colour_manual(name="Methods:",values=sel_cols, labels=sel_meth) +
  scale_shape_manual(name="Methods:",values=sel_pch, labels=sel_meth) +
  scale_fill_manual(name="Methods:",values=sel_cols, labels=sel_meth) +
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + scale_y_log10()

ggsave(file.path('~/Desktop/',paste0('sim_corRec.pdf')), combined_plot, height=2.5, width=9)

show(combined_plot)


# n fact ----

if(F){
  n_fact_g <- as.data.frame(n_fact) %>% group_by(n, method) %>% 
    summarize(K=median(K,na.rm=T), K1=median(K1,na.rm=T),
              K2=median(K2,na.rm=T), K3=median(K3,na.rm=T)) %>%
    mutate(
      K  = if_else(method == 2, K - 1, K),
      K  = if_else(method == 1, K - 1, K),
      K1 = if_else(method == 1, K1 - 1, K1),
      K2 = if_else(method == 1, K2 - 1, K2),
      K3 = if_else(method == 1, K3 - 1, K3)) %>%
    mutate(Ktot = K+K1+K2+K3)
  
  time_run_g <- as.data.frame(time_run) %>% group_by(n, method) %>% 
    summarize(timeRun=median(timeRun,na.rm=T)) %>%
    mutate(timeRun=round(timeRun,1))
  
  ess_mcmc_g <- as.data.frame(ess_mcmc) %>% group_by(n, method) %>% 
    summarize(nTot=median(nTot,na.rm=T), nBurn=median(nBurn,na.rm=T)) 
  
  merged_nfact <- n_fact_g %>%
    left_join(time_run_g, by = c("n","method")) %>%
    left_join(ess_mcmc_g,     by = c("n","method")) %>% 
    mutate(nTot = if_else(method == 4, 6000, nTot),
           nBurn = if_else(method == 4, 3000, nBurn)) %>%
    mutate(method = recode(method, `1` = "JAFAR", `2` = "JFR", `3` = "BSFP", `4` = "BIP")) %>%
    select(n, method, nTot, nBurn, Ktot, timeRun, everything(), -K, -K1, -K2, -K3)
  View(merged_nfact)
  
  # Basic LaTeX table
  merged_nfact %>%
    kable(format = "latex", booktabs = TRUE, 
          digits = 1,  # rounds numeric columns
          caption = "MCMC runtime vs samples and latent dimensions by method and n") %>%
    kable_styling(latex_options = c("hold_position", "striped"))
}


# ess ----

if(F){

  ## (0) drop unused columns
  time_run2 <- as.data.frame(time_run) %>% group_by(n, method) %>% 
    summarize(timeRun=median(timeRun,na.rm=T)) 
  
  ess_mcmc2 <- as.data.frame(ess_mcmc) %>% group_by(n, method) %>% 
    summarize(ess_s2=median(ess_s2,na.rm=T),
              nTot=median(nTot,na.rm=T),
              nBurn=median(nBurn,na.rm=T),
              nThin=median(nThin,na.rm=T)) %>%
    mutate(nEff = (nTot - nBurn) / nThin)
  
  ess_views2 <- as.data.frame(ess_views) %>% group_by(n, method,m) %>%
    summarize(median_ess=median(median_ess,na.rm=T))
  
  ## (0.1) wide format for ess_views: one col per m
  ess_views_wide <- ess_views2 %>%
    pivot_wider(id_cols = c(n, method),
                names_from = m,
                values_from = median_ess,
                names_prefix = "ess_m")
  
  merged <- ess_mcmc2 %>%
    left_join(ess_views_wide, by = c("n","method")) %>%
    left_join(time_run2,     by = c("n","method"))
  
  merged_cut <- merged %>% mutate(timeEff = timeRun * (nTot - nBurn) / nTot) %>%
    select(-nThin)         
  View(merged_cut)
  
  merged_cut <- merged %>% mutate(timeEff = timeRun * (nTot - nBurn) / nTot) %>%
    mutate(ess_s2=round(100*ess_s2/nEff,1), ess_m1=round(100*ess_m1/nEff,1),
           ess_m2=round(100*ess_m2/nEff,1), ess_m3=round(100*ess_m3/nEff,1)) %>%
    select(-nThin, -nTot, -nBurn, -timeRun) %>%
    select(n, method, nEff, timeEff, everything(), -ess_s2, ess_s2) %>% 
    mutate(method = recode(method, `1` = "JAFAR", `2` = "JFR", `3` = "BSFP")) %>%
    mutate(timeEff = round(timeEff,1))
  
  View(merged_cut)
  
  # Basic LaTeX table
  merged_cut %>%
    kable(format = "latex", booktabs = TRUE, 
          digits = 1,  # rounds numeric columns
          caption = "MCMC ESS and time efficiency by method and n") %>%
    kable_styling(latex_options = c("hold_position", "striped"))
}
