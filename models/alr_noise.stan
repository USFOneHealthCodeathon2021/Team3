// estimate ancestral node locations on a sphere assuming presence of large outlier movements
functions {
  real wrapped_cauchy_lpdf(real y, real mu, real gamma) {
    return(log(sinh(gamma))
           - log(cosh(gamma) - cos(y-mu))
           - log(2 * pi()));
  }
  real dist_sphere(phi1, phi2, theta1, theta2) {
    real deltaPhi = phi1 - phi2;
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
  vector<lower=-pi(), upper=pi()>[N] phi_mu; // tip longitude estimate in radians
  vector<lower=-pi()/2, upper=pi()/2>[N] theta_mu; // tip latitude estimate in radians
  vector<lower=0>[N] kappa; // precision of geocoding estimate (1/variance)
  vector<lower=0>[N+NI-1] time; //temporal distances between adjacent tips/nodes
  int self[N+NI-1]; // index for each tip/node in a single vector
  int ancestor[N+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_prior;
}
parameters {
  real<lower=0> sigma_raw; //dispersal rate
  unit_vector[3] loc[N+NI]; //ancestral node location
}
transformed parameters {
  vector[N+NI] phi; //longitude in radians
  vector[N+NI] theta; //latitude in radians
  vector[N] e; //error distances between geocoded and estimated node locations
  vector[N+NI-1] d; //geographic distances between adjacent tips/nodes
  for(n in 1:(N+NI)) {
    phi[n] = atan2(loc[n,2], loc[n,1]);
    theta[n] = atan2(sqrt(dot_self(loc[n,1:2])), loc[n,3]);
  }
  for(n in 1:N) {
    e[n] = dist_sphere(phi[n], phi_mu[n], theta[n], theta_mu[n]);
  }
  for(n in 1:(N+NI-1)) {
    d[n] = dist_sphere(phi[ancestor[n]], phi[self[n]], theta[ancestor[n]], theta[self[n]]);
  }
}
model {
  sigma_raw ~ std_normal();
  for(n in 1:N)
    e[n] ~ von_mises(0, kappa[n]);
  for(n in 1:(N+NI-1))
    d[n] ~ wrapped_cauchy(0, sqrt(time[n]) * sigma_prior * sigma_raw);
}
