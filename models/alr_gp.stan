// estimate ancestral node locations on a sphere assuming presence of large outlier movements (e.g. flights) and noisy geocoded tip estimates
// predict covariate values at estimated node locations
// estimate association between covariates and both likelihood of sampling locations (e.g. not likely to occur in low population density such as oceans) and number of mutations leading to nodes (the more interesting part)
// should generalize to multiple covariates so for instance important confounders like pop density can be accounted for 
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
  vector dist_sphere(vector phi1, vector phi2, vector theta1, vector theta2) {
    int N = rows(phi1);
    vector[N] deltaPhi = phi1 - phi2;
    vector[N] cdP = cos(deltaPhi);
    vector[N] ct2 = cos(theta2);
    vector[N] ct1 = cos(theta1);
    vector[N] st2 = sin(theta2);
    vector[N] st1 = sin(theta1);
    vector[N] y = sqrt(square(ct2 .* sin(deltaPhi)) + square(ct1 .* st2 - st1 .* ct2 .* cdP));
    vector[N] x = st1 .* st2 + ct1 .* ct2 .* cdP;
    vector[N] d;
    for(n in 1:N) d[n] = atan2(y[n], x[n]);
    return(d);
  }
  matrix fill_sym(vector lt, int N, real c) {
    matrix[N,N] s_mat;
    int iter = 1;
    for(j in 1:(N-1)) {
      s_mat[j,j] = c;
      for(i in (j+1):N) {
        s_mat[i,j] = lt[iter];
        s_mat[j,i] = lt[iter];
        iter += 1;
      }
    }
    s_mat[N,N] = c;
    return(s_mat);
  }
  matrix gp_L(vector d, real alpha, real rho, real sigma) {
    int N = cols(d);
    matrix[N,N] L 
      = cholesky_decompose(
          fill_sym(square(alpha_gp) * exp(-square(d) / (2 * square(rho))),
                   N,
                   square(alpha_gp) + square(sigma)));
  }
}
data {
  int<lower=2> NT; // number of tips
  int<lower=1> NI; // number of internal nodes
  int<lower=1> NY; // number of covariate observations
  matrix[3,NT] loc_mu; // tip location geocoded estimates
  vector<lower=0>[NT] kappa; // precision of geocoding estimate (1/variance)
  vector<lower=0>[NT+NI-1] time; // time between adjacent tips/nodes
  int self[NT+NI-1]; // index for each tip/node in a single vector
  int ancestor[NT+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_d_prior; // prior dispersal rate estimate
  vector<lower=-pi(), upper=pi()>[NY] phi_y; // covariate longitude 
  vector<lower=0, upper=pi()>[NY] theta_y; // covariate latitude 
  vector[NY] y_obs; // observed covariates
  int mut[NT+NI-1]; // number of mutations along a branch
  real<lower=0> rho_prior; // expected GP length-scale
  real<lower=0> scale_gp_prior; // expected scale of covariate
}
transformed data {
  int NN = NT + NI; // number of phylogenetic data points (tips plus nodes)
  int N = NN + NY; // number of geographic points
  int ND = choose(N,2); // number of pairwise geographic distances
  int present[N] = append_array(rep_array(1,NN), rep_array(0,NY));
  vector[NN-1] stime = sqrt(time);
  vector[NN-1] ltime = log(time);
  vector[NN] kdonut = rep_vector(fmax(mean(kappa),100), NN);
  for(n in 1:N) kdonut[n] = fmax(kappa[n],100);
}
parameters {
  real<lower=0> sigma_d_raw; // normalized dispersal rate
  vector[3] loc[NN]; // tip and node locations, unnormalized vector
  real alpha_y; // covariate mean
  vector[NN] y_pred; // interpolated covariate values
  real alpha_s; // mean likelihood of sampling a virus among all coordinates in model
  real beta_s; // logistic effect of covariate on the likelihood of sampling a virus
  real alpha_mut; // mean mutation rate
  vector[2] beta_mut; // log effect of time and covariate on mutation rate
  real<lower=0> rho_raw; // normalized length scale of GP correlation
  real<lower=0> alpha_gp_raw; // normalized covariate variance explained by GP
  real<lower=0> sigma_gp_raw; // normalized covariate residual variance
}
transformed parameters {
  vector[NN] normalizing;
  matrix[3,NN] loc; // tip and node locations
  vector[N] y = append_row(y_obs, y_pred); // observed and estimated values of covariates
  vector[N] phi; // tip/node longitudes in radians
  vector[N] theta; // tip/node latitudes in radians
  real<lower=0> rho = rho_raw * rho_prior; // gp length scale
  real<lower=0> alpha_gp = alpha_gp_raw * scale_gp_prior; // covariate variance explained by GP
  real<lower=0> sigma_gp = sigma_gp_raw * scale_gp_prior; // covariate residual variance
  for(n in 1:(NN)) {
    normalizing[n] = sqrt(dot_self(locvec[n]));
    loc[,n] = locvec[n] / normalizing[n];
    phi[n] = atan2(loc[2,n], loc[1,n]);
    theta[n] = atan2(sqrt(dot_self(loc[1:2,n])), loc[3,n]); 
  }
  phi[(NN+1):N] = phi_y;
  theta[(NN+1):N] = theta_y;
}
model {
  vector[ND] d = dist_sphere(phi[dInd1], phi[dInd2], theta[dInd1], theta[dInd2]); // pairwise geographic distances of all but geocode estimates
  matrix[N,N] L_y = gp_L(d, alpha_gp, rho, sigma_gp); 
  sigma_d_raw ~ std_normal();
  beta_s ~ std_normal();
  beta_mut ~ std_normal();
  rho_raw ~ inv_gamma(5, 5);
  alpha_gp_raw ~ std_normal();
  sigma_gp_raw ~ std_normal();
  normalizing ~ gamma(kdonut, kdonut); // keep cartesian parameters near surface of unit sphere
  target += -sum(2 * log(normalizing)); // unit vector jacobian
  loc[,1:N] ~ fisher(loc_mu, kappa); // actual tip locations constrained within geocoding uncertainty
  loc[,self] ~ spherical_cauchy(loc[,ancestor], link_scale(sigma_d_prior * sigma_d_raw * stime)); // ancestral node locations shrink toward location of direct descendants
  y ~ multi_normal_cholesky(alpha_y, L_y); // gaussian process fit to observed covariates and predicts values at tip/node locations
  present ~ bernoulli_logit(alpha_s + beta_s * y); // likelihood of sampling a virus at a given location as function of covariates
  mut ~ poisson_log_glm(append_col(ltime, y[NY + self]), alpha_mut, beta_mut); // number of mutations as a function of time and covariates
}
