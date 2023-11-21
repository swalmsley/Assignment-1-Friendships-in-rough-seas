
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Options
tar_option_set()

# Seed
tar_option_set(seed = 1234)

# Variables
set_cmdstan_path(path='C:/Users/sjfwa/AppData/Local/R/cmdstan-2.33.1') # cmdstan path

# Targets
list(
  
  ### ### ### ### ### 
  # load in data ### ### ### ### Q: file system here possible with operation? relevant for stan models in cmdstanr approach too...
  tar_target(raw_photoData, read_data('./input/LV_SS_Master_1988-2021.csv')),
  tar_target(raw_envData, read_data('./input/Scotian_Shelf_1988-2019_edited_HW_environment.csv')),
  
  
  # process photo data
  tar_target(side_raw, process_photoID(raw_photoData, 'Left', 'Gully')), # only left-sided photos from the Gully
  tar_target(side, side_raw[year==2016,,]), # trimming to 2019 for now
  
  
  # prepare binary association data
  tar_target(dyad_list, build_binary_df(side)),
  tar_target(dyad_data, dyad_list[bothPresent==1,together:=binary_association(side[day==day,,], A, B, 10), by='index']), # association within 10 minutes
  tar_target(assoc_data, dyad_data[is.na(together),together:=0,]), # convert remaining NAs to '0's
  tar_target(assoc_alive, both_catalogue_alive_assocData(assoc_data, side)), # only when both animals are known to be alive
  tar_target(assoc_labelled, assoc_alive[,dyad_label:=.GRP,by=c('A','B')]), # add integer labels for each dyad (necessary for Stan model)
  
  # merge in sea state at level of observation
  tar_target(assoc_final, merge(assoc_labelled, raw_envData[Time==0.5,c('Date','Beaufort'),], by.x='day',by.y='Date')),
  
  # plot DAG
  tar_target(Figure1, save_figure('./Figures/Figure1_DAG.png',w=4,h=4,create_dag())),
  

  tar_target(data_list, list(N = nrow(assoc_final),
                             N_dyads = length(unique(assoc_final$dyad_label)),
                             Together = assoc_final$together,
                             SeaState = standardize(assoc_final$Beaufort),
                             dyad = assoc_final$dyad_label)),
  tar_stan_mcmc(ewm, './Stan-Models/m1.stan',
                data=data_list,
                parallel_chains=4),
  
  
  # plot diagnostics
  tar_target(trankplot, mcmc_rank_overlay(as_draws(ewm_mcmc_m1), pars=c('beta', 'a_bar', 'sigma', 'a[1]','a[10]','a[100]'))),
  tar_target(Figure3, save_figure('./Figures/Figure3_trankplot.png',w=8,h=4,trankplot)),
  
  
  # pull out edges
  tar_target(edge_params, data.table(ewm_summary_m1[2085:4164,,])),
  tar_target(edges, edge_params[,weight:=inv_logit(mean),]),
  
  
  # next, merge on true proportions for comparison, if I decide I need to plot those
  tar_target(Figure2, save_figure('./Figures/Figure2_priorPredictive.png',w=8,h=4,prior_predictive())),
  
  
  # visualize model results
  tar_target(Figure4, save_figure('./Figures/Figure4_ShrinkagePanels.png',w=6,h=8,
                                  (shrinkPlot(assoc_final, edges, TRUE) / shrinkPlot(assoc_final, edges, FALSE) + plot_annotation(tag_levels = 'A')))),
  tar_target(Figure5, save_figure('./Figures/Figure5_SeaStateEffect.png',w=6,h=4,posterior_density(ewm_draws_m1$beta, 'Effect of Sea State'))),
  

  # write output
  tar_quarto(
    render,
    file.path('Submission', 'outputs-Assignment1.qmd'))

)
