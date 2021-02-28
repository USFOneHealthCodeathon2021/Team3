// estimate ancestral node locations on a sphere assuming presence of large outlier movements (e.g. flights) and noisy geocoded tip estimates
// predict covariate values at estimated node locations
// estimate association between covariates and both likelihood of sampling locations (e.g. not likely to occur in low population density such as oceans) and number of mutations leading to nodes (the more interesting part)
// should generalize to multiple covariates so for instance important confounders like pop density can be accounted for 
// donut prior based on https://discourse.mc-stan.org/t/divergence-treedepth-issues-with-unit-vector/8059/3
// Fisher distribution inspired by same forum
// spherical cauchy in: https://arxiv.org/pdf/1510.07679.pdf
functions {
  real spherical_cauchy_dots_lpdf(vector dots, vector scale) {
    real ss = square(scale);
    return(3 * sum(log1m(ss) - log1p(ss - 2 * scale .* dots)));
  } // dropped constants
  real fisher_lpdf(matrix y, matrix mu, vector kappa) {
    return(sum(log(kappa)
               - log(2 * pi())
               - log_diff_exp(kappa, -kappa)
               + kappa .* columns_dot_product(y,mu)'));
  }
  vector link_scale(vector scale) {
    return(scale ./ (scale + 1));
  } // not sure about priors with this link
  dist_sphere_dots(vector dots) {
    int ND = rows(dots);
    vector[ND] d;
    for(n in 1:ND) {
      d[n] = atan2(sqrt(1-square(dots[n])), dots[n]);
    }
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
  matrix gp_L(vector d, int N, real alpha, real rho, real sigma) {
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
  matrix[3,NY] loc_y; // covariate locations, cartesian unit
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
  int dInd1[ND]; // indices to map pairwise distances
  int dInd2[ND];
  int dInd[N,N];	
  int dEdgeInd[NN-1]; // indices to select distances between adjacent nodes
  for(n in 1:N) kdonut[n] = fmax(kappa[n],100);
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
  matrix[3,N] loc; // tip and node locations
  vector[N] y = append_row(y_pred, y_obs); // observed and estimated values of covariates
  real<lower=0> rho = rho_raw * rho_prior; // gp length scale
  real<lower=0> alpha_gp = alpha_gp_raw * scale_gp_prior; // covariate variance explained by GP
  real<lower=0> sigma_gp = sigma_gp_raw * scale_gp_prior; // covariate residual variance
  for(n in 1:(NN)) {
    normalizing[n] = sqrt(dot_self(locvec[n]));
    loc[,n] = locvec[n] / normalizing[n];
  }
  loc[,(NN+1):N] = loc_y;
}
model {
  vector[ND] dots = columns_dot_product(loc[,dInd1], loc[,dInd2])';
  matrix[N,N] L_y = gp_L(dist_sphere_dots(dots), N, alpha_gp, rho, sigma_gp); 
  sigma_d_raw ~ std_normal();
  beta_s ~ std_normal();
  beta_mut ~ std_normal();
  rho_raw ~ inv_gamma(5, 5);
  alpha_gp_raw ~ std_normal();
  sigma_gp_raw ~ std_normal();
  normalizing ~ gamma(kdonut, kdonut); // keep cartesian parameters near surface of unit sphere
  target += -sum(2 * log(normalizing)); // unit vector jacobian
  loc[,1:NT] ~ fisher(loc_mu, kappa); // actual tip locations constrained within geocoding uncertainty
  dots[dEdgeInd] ~ spherical_cauchy_dots(link_scale(sigma_d_prior * sigma_d_raw * stime)); // ancestral node locations shrink toward location of direct descendants
  y ~ multi_normal_cholesky(alpha_y, L_y); // gaussian process fit to observed covariates and predicts values at tip/node locations
  present ~ bernoulli_logit(alpha_s + beta_s * y); // likelihood of sampling a virus at a given location as function of covariates
  mut ~ poisson_log_glm(append_col(ltime, y[self]), alpha_mut, beta_mut); // number of mutations as a function of time and covariates
}
