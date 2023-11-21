
# Custom functions to be called in '_targets.R'



################################################
# DAGs and graphics
################################################



# DAG ---------------------------------------------------------------------
create_dag <- function() {
  dag <- dagitty("dag {
  SocialTie_dyadic -> Together_unobserved
  Together_unobserved -> Together_observed
  SeaState -> Together_observed
                 }")
  
  coordinates(dag) <- list( x=c(SocialTie_dyadic=0, SeaState=0, Together_unobserved=1, Together_observed=1),
                            y=c(SocialTie_dyadic=1, SeaState=2, Together_unobserved=1, Together_observed=2))
  
  return(drawdag(dag))
}



################################################
# General functions
################################################


# Read data ---------------------------------------------------------------
read_data <- function(path) {
  d <- data.table(read.csv(path))
  return(d)
}
# testing



# Save figure -------------------------------------------------------------
save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=1000)
  print(call)
  dev.off()
}


# Posterior density -------------------------------------------------------
posterior_density <- function(samples,xlab) {
  
  dt <- data.frame(samples)
  ggplot(dt, aes(x=samples)) +
    geom_density(color='white',linewidth=1, fill=rangi2, alpha=0.5) +
    geom_vline(xintercept=0, linetype='dashed')+
    labs(x=xlab, y='Density')+
    theme_classic()
  
}


################################################
# Project-specific functions
################################################



# Process photo-ID observations -------------------------------------------
process_photoID <- function(lv_data, chosen_side, chosen_canyon){
  
  # check that Gully is keyworded appropriately in 2021 data
  
  # extract side information
  lv_data[,side:=ifelse(grepl("Left Side", Keyword.export),"Left","Right"),] # extract side from keyword list
  lv_data[,numSides:=length(unique(side)),by=Title] # identify IDs that are 2-sided
  
  # fix up date information
  lv_data[,dateTime:=as.POSIXct(Date.Original),] # reformat dateTimes
  lv_data[,day:=format(dateTime, format = "%Y-%m-%d")]
  lv_data[,year:=year(day),] # create a column of unique dates for use later
  lv_data[,numDaysByYear:=length(unique(day)),by=year]
  
  # remove observations with multiple or unknown whales
  lv_data <- lv_data[Title!="see crops",,] 
  lv_data <- lv_data[Title!="unk",,] 
  
  # add extra columns to indicate number of observations
  lv_data[,obsByYear:=.N,by=year] # count number of observations by year
  lv_data[,numDaysByTitle:=length(unique(day)),by=Title] # count number of days each individual was detected
  
  # add columns 
  lv_data[, canyon:=NA,]
  lv_data[, canyon:=ifelse(grepl('Gully', Keyword.export),'Gully',NA),]
  lv_data[, canyon:=ifelse(grepl('Shortland', Keyword.export),'Shortland',canyon),]
  lv_data[, canyon:=ifelse(grepl('Haldimand', Keyword.export),'Haldimand',canyon),]
  ## Note there are 1609 photos without canyon information
  # many of these are multiple whales photos, others are not... in any case this should be remedied in the next iteration of the master catalogue
  # also some UNK - need to decide what to do with these
  
  
  # subset by location and side  
  output <- lv_data[canyon==chosen_canyon,,]
  output <- output[side==chosen_side,,]
  
  return(output)
  
}


# Build binary dataframe --------------------------------------------------
build_binary_df <- function(photoData) {
  
  # Note: currently hard-coding days as sampling periods
  dyads <- data.table(t(combn(photoData[,unique(Title),], 2))) # generate all dyads
  colnames(dyads) <- c('A', 'B') # rename ID columns
  df <- data.table(A = as.character(), B = as.character(), day = as.character(), A_observed = as.integer(), B_observed = as.integer())   # set up empty data frame for filling 
  
  for (d in unique(photoData$day)) { # create a new row for each observation day and each possible dyad in the dataset
    temp <- dyads[, day:=d, ]
    temp[, A_observed:=ifelse(A %in% photoData[day==d,unique(Title),],1,0)]
    temp[, B_observed:=ifelse(B %in% photoData[day==d,unique(Title),],1,0)]
    temp <- temp[A_observed==1 | B_observed==1,,] # prune out rows (day-dyad combinations) where neither member of dyad was observed
    df <- rbindlist(list(df, temp))
  }
  
  df[,index:=.I,] # add index for row number
  df[,bothPresent:=ifelse((A_observed==1 & B_observed==1),1,0),]
  
  return(df)
  
}


# Compute binary associations ---------------------------------------------
binary_association <- function(photos_sampling_period, ID1, ID2, time_limit) {
  
  # Pre-filter the data
  ID1_obs <- photos_sampling_period[(Title %in% ID1),,] # observations of ID 1 within sampling period
  ID2_obs <- photos_sampling_period[(Title %in% ID2),,] # observations of ID 2 within sampling period
  
  # Proceed only if there are any observations for ID1 or ID2 in the given sampling period (otherwise return 0 right away)
  if (nrow(ID1_obs) == 0 || nrow(ID2_obs) == 0) return(0)
  
  # Calculate the time difference between all combinations of observations
  time_diff <- outer(ID1_obs$dateTime, ID2_obs$dateTime, FUN = function(x, y) abs(as.numeric(difftime(x, y, units = 'mins'))))
  
  # Check if any time difference is within the time_limit
  if (any(time_diff <= time_limit)) return(1)
  
  return(0)
  
}



# Shrinkage plot ----------------------------------------------------------
shrinkPlot <- function(assoc_final, edges, zero_only) {
  
  assoc_final[,propTogether:=(sum(together)/.N),by=dyad_label]
  edges[,dyad_label:=as.numeric(str_extract(variable, "\\d+")),]
  edges <- merge(edges, unique(assoc_final[,c('dyad_label','propTogether', 'A','B')]))
  
  if (zero_only) (edges2 <- edges[propTogether==0,,])
  if (!zero_only) (edges2 <- edges[propTogether>0,,])
  
  sub <- edges2[sample(1:nrow(edges2), 10, replace=FALSE),,]
  sub[,dyad_index:=.I,]
  
  ggplot(sub, aes(x=as.factor(dyad_index), y=propTogether)) +
    
    geom_point(aes(y=weight),pch=1,color='black', size=3) +
    
    geom_errorbar(aes(y=weight,ymin=inv_logit(q5), ymax=inv_logit(q95)),width=.2,linewidth=0.75,alpha=0.75)+
    
    #geom_hline(yintercept=mean(edges$propTogether), color=rangi2) +
    
    geom_hline(yintercept=mean(edges$weight), color='black', linetype='dashed') +
    
    geom_hline(yintercept=mean(edges$propTogether), color=rangi2, linetype='dashed') +
    
    geom_point(pch=18,color=rangi2, size=4) +
    
    labs(x='Dyad index (random selection)', y='Proportion of time spent together') +
    
    ylim(0,1) + theme_classic()
  
}



# Coefficients table (cmdstan) --------------------------------------------
coef_table <- function(mcmc_object) {
  
  # create table
  precis_output <- precis(mcmc_object)
  stats.table <- data.table(data.frame(precis_output))
  stats.table <- cbind(rownames(precis_output), stats.table)
  #stats.table[,c('Rhat','Bulk_ESS','Tail_ESS'):=NULL,]
  names(stats.table) <- c("Parameter", "Mean", "SD", "CI-5.5", "CI-94.5", "R-hat", "ESS")
  
  return(nice_table(stats.table))
  
}


  
# Create prior predictive checks ------------------------------------------
prior_predictive <- function() {
  
  ## number of samples to simulate
  N <- 1000
  
  ## prior predictive for dyad effects
  a_bar_prior <- rnorm(N,-1,1.5)
  sigma_prior <- dexp(1)
  z_prior <- rnorm(N,0,1)
  dyad_prior <- a_bar_prior + z_prior*sigma_prior


  ## prior predictive for slope term
  
  beta <- rnorm(N,0,0.5)
  
  ## put it all together
  
  par(mfrow=c(1,3))
  
  plot(density(dyad_prior), xlim=c(-7,7), ylim=c(0,1), 
       xlab='Log odds scale',
       ylab='Density',
       main='Dyad effect (log odds scale)')
  
  plot(density(inv_logit(dyad_prior)), xlim=c(0,1), ylim=c(0,4), 
       xlab='Probability scale',
       ylab='Density',
       main='Dyad effect (prob. scale)')
  
  plot(density(beta), xlim=c(-3,3), ylim=c(0,1), 
       xlab='Log odds scale',
       ylab='Density',
       main='Beta (log odds scale)')
  
}




# Prune data to when both animals are known -------------------------------
both_catalogue_alive_assocData <- function(assocData, photoData) {
  
  assocData[,A_minYear:=photoData[Title==A,min(year),],by=A] # min years
  assocData[,B_minYear:=photoData[Title==B,min(year),],by=B]
  
  assocData[,A_maxYear:=photoData[Title==A,max(year),],by=A] # max years
  assocData[,B_maxYear:=photoData[Title==B,max(year),],by=B]
  
  assocData[,year:=year(day),] # add year column
  
  assocData[,bothAlive:=ifelse(all(year >= unique(A_minYear),
                                   year >= unique(B_minYear), 
                                   year <= unique(A_maxYear), 
                                   year <= unique(B_maxYear)),
                               TRUE,FALSE),
            by=c('A','B','year')]
  
  return(assocData)
  
}




print('Functions successfully loaded')