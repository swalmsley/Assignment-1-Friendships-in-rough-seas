
// DATA
data {
  int<lower=0> N; //number of observations
  int<lower=0> N_dyads; //number of dyads
  array[N] int Together; //whether association was detected
  vector[N] SeaState; //sea state
  array[N] int dyad; //dyad in question for given observation
}

// PARAMETERS
parameters {
  real a_bar; // mean of population of dyad-specific intercepts
  real<lower=0> sigma; // sd of population of dyad-specific intercepts
  vector[N_dyads] z; // z for each dyad
  real beta; // slope term sea state
} 

// MODEL
model {
  
  // internal parameters
  vector[N] p;
  
  // priors
  a_bar ~ normal(-1, 1.5);
  z ~ normal(0, 1);
  sigma ~ exponential(1);
  beta ~ normal(0, 0.5);

  // likelihood
  p = inv_logit(a_bar + z[dyad]*sigma  + beta*SeaState);
  
  Together ~ binomial(1, p); 
  
}

// GENERATED QUANTITIES
generated quantities {
  
  vector[N_dyads] a;
  a = a_bar + z*sigma;
  
}
