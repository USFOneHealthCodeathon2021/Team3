// estimate ancestral node locations on a sphere assuming presence of large outlier movements
// donut prior based on https://discourse.mc-stan.org/t/divergence-treedepth-issues-with-unit-vector/8059/3
// Fisher distribution inspired by same forum
// spherical cauchy in: https://arxiv.org/pdf/1510.07679.pdf
functions {
  real spherical_cauchy_lpdf(matrix y, matrix mu, vector scale) {
    real ss = square(scale);
    return(3 * sum(log1m(ss) - log1p(ss - 2 * scale .* columns_dot_product(y,mu)')));
  } // dropped constants
  real fisher_lpdf(matrix y, matrix mu, vector kappa) {
    return(sum(log(kappa)
               - log(2 * pi())
               - log_diff_exp(kappa, -kappa)
               + kappa .* columns_dot_product(y,mu)'));
  }
  vector link_scale(vector scale) {
    return(scale ./ (scale + 1));
  }
}
data {
  int<lower=2> N; // number of tips
  int<lower=1> NI; // number of internal nodes
  matrix[3,N] loc_mu; // tip location geocoded estimates
  vector<lower=0>[N] kappa; // precision of geocoding estimate (1/variance)
  vector<lower=0>[N+NI-1] time; // time between adjacent tips/nodes
  int self[N+NI-1]; // index for each tip/node in a single vector
  int ancestor[N+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_prior; // prior dispersal rate estimate
}
transformed data {
  vector[N+NI] kdonut = rep_vector(fmax(mean(kappa),100), N+NI);
  for(n in 1:N) kdonut[n] = fmax(kappa[n],100);
}
parameters {
  real<lower=0> sigma_raw; // normalized dispersal rate
  vector[3] locvec[N+NI]; // tip/node location, unnormalized vector
}
transformed parameters {
  vector[N+NI] normalizing;
  matrix[3,N+NI] loc; // tip/node location
  for(n in 1:(N+NI)) {
    normalizing[n] = sqrt(dot_self(locvec[n]));
    loc[,n] = locvec[n] / normalizing[n];
  }
}
model {
  sigma_raw ~ std_normal();
  normalizing ~ gamma(kdonut, kdonut);
  target += -sum(2 * log(normalizing)); // unit vector jacobian
  loc[,1:N] ~ fisher(loc_mu, kappa);
  loc[,self] ~ spherical_cauchy(loc[,ancestor], link_scale(sqrt(time) * sigma_prior * sigma_raw));
}
