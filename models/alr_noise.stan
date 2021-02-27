// estimate ancestral node locations on a sphere assuming presence of large outlier movements
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
}
data {
  int<lower=2> N; // number of tips
  int<lower=1> NI; // number of internal nodes
  vector<lower=-pi(), upper=pi()>[N] phi_mu; // tip longitude geocoded estimates in radians
  vector<lower=0, upper=pi()>[N] theta_mu; // tip latitude geocoded estimates in radians
  vector<lower=0>[N] kappa; // precision of geocoding estimate (1/variance)
  vector<lower=0>[N+NI-1] time; // time between adjacent tips/nodes
  int self[N+NI-1]; // index for each tip/node in a single vector
  int ancestor[N+NI-1]; // index for each tip/node's ancestor in a single vector
  real<lower=0> sigma_prior; // prior dispersal rate estimate
}
parameters {
  real<lower=0> sigma_raw; // normalized dispersal rate
  unit_vector[3] loc[N+NI]; // tip/node location
}
transformed parameters {
  vector[N+NI] phi; // tip/node longitudes in radians
  vector[N+NI] theta; // tip/node latitudes in radians
  vector[N] e; // geographic distances between geocoded and estimated tip locations
  vector[N+NI-1] d; // geographic distances between adjacent tips/nodes
  for(n in 1:(N+NI)) {
    phi[n] = atan2(loc[n,2], loc[n,1]);
    theta[n] = atan2(sqrt(dot_self(loc[n,1:2])), loc[n,3]);
  }
  e = dist_sphere(phi[1:N], phi_mu, theta[1:N], theta_mu);
  d = dist_sphere(phi[ancestor], phi[self], theta[ancestor], theta[self]);
}
model {
  sigma_raw ~ std_normal();
  e ~ von_mises(0, kappa);
  d ~ wrapped_cauchy(sqrt(time) * sigma_prior * sigma_raw);
}
