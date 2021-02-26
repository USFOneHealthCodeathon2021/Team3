// estimate ancestral node locations on a sphere assuming presence of large outlier movements
functions {
  real wrapped_cauchy_lpdf(real y, real mu, real gamma) {
    return(log(sinh(gamma))
           - log(cosh(gamma) - cos(y-mu))
           - log(2 * pi()));
  }
}
data {
  int<lower=2> N; // number of tips
  int<lower=1> NI; //number of internal nodes
  vector<lower=-pi(), upper=pi()>[N] phi_tips; //longitude in radians
  vector<lower=-pi()/2, upper=pi()/2>[N] theta_tips; //latitude in radians
  vector<lower=0>[N+NI-1] time; //temporal distances between adjacent tips/nodes
  int self[N+NI-1]; // index for each tip/node in a single vector
  int ancestor[N+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_prior;
}
parameters {
  real<lower=0> sigma_raw; //dispersal rate
  unit_vector[3] loc_anc[NI]; //ancestral node location
}
transformed parameters {
  vector[N+NI] phi; //longitude in radians
  vector[N+NI] theta; //latitude in radians
  vector[N+NI-1] d; //geographic distances between adjacent tips/nodes
  phi[1:N] = phi_tips;
  theta[1:N] = theta_tips;
  for(n in 1:NI) {
    phi[N+n] = atan2(loc_anc[n,2], loc_anc[n,1]);
    theta[N+n] = atan2(sqrt(dot_self(loc_anc[n,1:2])), loc_anc[n,3]);
  }
  for(n in 1:(N+NI-1)) {
    real deltaPhi = phi[self[n]] - phi[ancestor[n]];
    d[n] = atan2(sqrt(square(cos(theta[self[n]]) * sin(deltaPhi))
                      + square(cos(theta[ancestor[n]]) * sin(theta[self[n]])
                               - sin(theta[ancestor[n]]) * cos(theta[self[n]]) * cos(deltaPhi))),
                 sin(theta[ancestor[n]]) * sin(theta[self[n]])
                 + cos(theta[ancestor[n]]) * cos(theta[self[n]]) * cos(deltaPhi));
  }
}
model {
  sigma_raw ~ std_normal();
  for(n in 1:(N+NI-1))
    d[n] ~ wrapped_cauchy(0, sqrt(time[n]) * sigma_prior * sigma_raw);
}
