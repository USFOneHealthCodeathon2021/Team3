// estimate ancestral node locations on a sphere assuming presence of large outlier movements (e.g. flights) and noisy geocoded tip estimates
// predict covariate values at estimated node locations
// estimate association between covariates and both likelihood of sampling locations (e.g. not likely to occur in low population density such as oceans) and number of mutations leading to nodes (the more interesting part)
// should generalize to multiple covariates so for instance important confounders like pop density can be accounted for 
functions {
  real wrapped_cauchy_lpdf(vector d, vector scale) {
    return(sum(log(sinh(scale))
               - log(cosh(scale) - cos(d)))
           - rows(d) * log(2 * pi()));
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
  vector<lower=-pi(), upper=pi()>[NT] phi_mu; // tip longitude geocoded estimates in radians
  vector<lower=0, upper=pi()>[NT] theta_mu; // tip latitude geocoded estimates in radians
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
  int dInd1[ND]; // indices to map pairwise distances
  int dInd2[ND];
  int dInd[N,N];
  int dEdgeInd[NN-1]; // indices to select distances between adjacent nodes
  int iter = 1;
  for(j in 1:(N-1)) {
    for(i in (j+1):N) {
      dInd1[iter] = i;
      dInd2[iter] = j;
      dInd[i,j] = iter;
      dInd[j,i] = iter;
      iter += 1;
    }
  }
  for(n in 1:(NN-1)) {
    dEdgeInd[n] = dInd[self[n], ancestor[n]];
  }
}
parameters {
  real<lower=0> sigma_d_raw; // normalized dispersal rate
  unit_vector[3] loc[NN]; // tip and node locations
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
  vector[N] y = append_row(y_obs, y_pred); // observed and estimated values of covariates
  vector[N] phi; // tip/node longitudes in radians
  vector[N] theta; // tip/node latitudes in radians
  real<lower=0> rho = rho_raw * rho_prior; // gp length scale
  real<lower=0> alpha_gp = alpha_gp_raw * scale_gp_prior; // covariate variance explained by GP
  real<lower=0> sigma_gp = sigma_gp_raw * scale_gp_prior; // covariate residual variance
  for(n in 1:(NN)) {
    phi[n] = atan2(loc[n,2], loc[n,1]);
    theta[n] = atan2(sqrt(dot_self(loc[n,1:2])), loc[n,3]); 
  }
  phi[(NN+1):N] = phi_y;
  theta[(NN+1):N] = theta_y;
}
model {
  vector[ND] d = dist_sphere(phi[dInd1], phi[dInd2], theta[dInd1], theta[dInd2]); pairwise geographic distances of all but geocode estimates
  vector[NT] d_err = dist_sphere(phi[1:NT], phi_mu, theta[1:NT], theta_mu); // geographic distances between geocoded and estimated tip locations
  matrix[N,N] L_y = gp_L(d, alpha_gp, rho, sigma_gp); 
  sigma_d_raw ~ std_normal();
  beta_s ~ std_normal();
  beta_mut ~ std_normal();
  rho_raw ~ inv_gamma(5, 5);
  alpha_gp_raw ~ std_normal();
  sigma_gp_raw ~ std_normal();
  d_err ~ von_mises(0, kappa); // actual tip locations constrained within geocoding uncertainty
  d[dEdgeInd] ~ wrapped_cauchy(sigma_d_prior * sigma_d_raw * stime); // ancestral node locations shrink toward location of direct descendants
  y ~ multi_normal_cholesky(alpha_y, L_y); // gaussian process fit to observed covariates and predicts values at tip/node locations
  present ~ bernoulli_logit(alpha_s + beta_s * y); // likelihood of sampling a virus at a given location as function of covariates
  mut ~ poisson_log_glm(append_col(ltime, y[NY+self]), alpha_mut, beta_mut); // number of mutations as a function of time and covariates
}
