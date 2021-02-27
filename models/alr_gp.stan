// estimate ancestral node locations on a sphere assuming presence of large outlier movements
functions {
  real wrapped_cauchy_lpdf(real y, real mu, real gamma) {
    return(log(sinh(gamma))
           - log(cosh(gamma) - cos(y-mu))
           - log(2 * pi()));
  }
  real dist_sphere(real phi1, real phi2, real theta1, real theta2) {
    real deltaPhi = phi1 - phi2;
    real d;
    d = atan2(sqrt(square(cos(theta2) * sin(deltaPhi))
                   + square(cos(theta1) * sin(theta2)
                            - sin(theta1) * cos(theta2) * cos(deltaPhi))),
              sin(theta1) * sin(theta2)
              + cos(theta1) * cos(theta2) * cos(deltaPhi));
    return(d);
  }
}
data {
  int<lower=2> N; // number of tips
  int<lower=1> NI; //number of internal nodes
  int<lower=1> NY; //number of covariate observations
  vector<lower=-pi(), upper=pi()>[N] phi_mu; // tip longitude estimate in radians
  vector<lower=-pi()/2, upper=pi()/2>[N] theta_mu; // tip latitude estimate in radians
  vector<lower=0>[N] kappa; // precision of geocoding estimate (1/variance)
  vector<lower=0>[N+NI-1] time; //temporal distances between adjacent tips/nodes
  int self[N+NI-1]; // index for each tip/node in a single vector
  int ancestor[N+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_prior;
  vector[NY] y_obs; // observed covariates
  int mut[N+NI-1]; // number of mutations along a branch
}
parameters {
  real<lower=0> sigma_raw; //dispersal rate
  unit_vector[3] loc[N+NI]; //tip and node locations
  vector[N+NI+NY] s; // log likelihood of sampling a virus at a given location
  vector[N+NI] y_pred; // interpolated covariate values
}
transformed parameters {
  vector[N+NI] phi; //longitude in radians
  vector[N+NI] theta; //latitude in radians
  vector[N] e; //error distances between geocoded and estimated node locations
  vector[N+NI-1] d; //geographic distances between adjacent tips/nodes
  vector[NY+N+NI] y = append_row(y_obs, y_pred); // observed and estimated values of covariates
  for(n in 1:(N+NI)) {
    phi[n] = atan2(loc[n,2], loc[n,1]);
    theta[n] = atan2(sqrt(dot_self(loc[n,1:2])), loc[n,3]); 
  }
  for(n in 1:N) {
    e[n] = dist_sphere(phi[n], phi_mu[n], theta[n], theta_mu[n]); // can this be vectorized?
  }
  for(n in 1:(N+NI-1)) {
    d[n] = dist_sphere(phi[ancestor[n]], phi[self[n]], theta[ancestor[n]], theta[self[n]]);
  }
}
model {
  matrix[NY+N+NI,NY+N+NI] L_y = cholesky_decompose(gp_kernel()); // fill in rest of gp code (estimate sigma)
  matrix[N+NI,N+NI] L_s = cholesky_decompose(gp_kernel()); // fill in rest of gp code (no direct obs; no sigma)
  sigma_raw ~ std_normal();
  e ~ von_mises(0, kappa);
  for(n in 1:(N+NI-1))
    d[n] ~ wrapped_cauchy(0, sqrt(time[n]) * sigma_prior * sigma_raw); // can this be vectorized?
  1 ~ bernoulli_logit(segment(s, 1, N+NI)); // observation for every phylogenetic node
  0 ~ bernoulli_logit(segment(s, N+NI+1, NY)); // no observation at every point without a phylogenetic node
  y ~ multi_normal_cholesky(zeros_vector(N), L_y); // gaussian process fit to observed covariates and predict values at node locations (should i model an intercept?)
  s ~ multi_normal_cholesky(alpha + beta * y, L_s); // gaussian process to estimate smooth probability of sampling at a given location based on covariates (double check efficiency of logistic gp in manual)
  mut ~ poisson_log(intercept + log(time) + zeta * y); // number of mutations as a function of time and covariates (maybe use builtin glm function - would it be more honest to model an overdispersion parameter, perhaps even with gp?)
}
