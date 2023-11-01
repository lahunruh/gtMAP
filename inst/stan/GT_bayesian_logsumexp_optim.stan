data {
  int<lower=0> N;
  int<lower=0> T;
  vector[T] Ct;
  matrix[T, N] M;
  int<lower=0> n;
  real<lower=0> sigma;
  real alpha;
  real beta;
  real a;
  real b;
  real m;
  real s;
  int w[N];
  real<lower=0, upper=1> p;
  real ctbase;
  real ctmult;
}

transformed data {
  real x_r[1] = rep_array(0.0,1);
}

parameters {
  vector<lower=0>[N] Ct_ind;
  //real<lower=0, upper=1> p;
}

transformed parameters {
  vector[T] Ct_pred;

  Ct_pred = rep_vector(0, T);

  for (i in 1:T) {
    for (j in 1:N) {
      if (M[i,j] > 0) {
        if (w[j] > 0) {
          Ct_pred[i] = Ct_pred[i] + M[i,j] * exp2(ctbase - Ct_ind[j]) * ctmult;
        }
      }
    }
  }

  //Ct_pred = log2(Ct_pred);
  for (i in 1:T) {
    if (Ct_pred[i] > 0) {//negative_infinity() ==
        Ct_pred[i] = ctbase - log2(Ct_pred[i]/(ctmult * sum(M[i,]))); // added 37 - and /(100 *
    }
  }
}

model {
  int sumW;

  sumW = sum(w);

  for (i in 1:N) {
    if (w[i] != 0) {
      target += normal_lpdf(Ct_ind[i] | m, s);//Ct_lpdf(Ct_ind[i] | a,b,m,s,x_r);//
    }
  }

  //sumW ~ binomial(121, p); // num_elements(w), p);
  //w ~ bernoulli(p); // deactivate if prevalence not wanted
  //p ~ beta(alpha, beta);
  Ct ~ normal(Ct_pred, sigma);

}

generated quantities {
  real total;
  vector[T] Ct_pred_out;

  Ct_pred_out = rep_vector(0, T);

  for (i in 1:T) {
    for (j in 1:N) {
      if (M[i,j] > 0) {
        if (w[j] > 0) {
          Ct_pred_out[i] = Ct_pred_out[i] + M[i,j] * exp2(ctbase - Ct_ind[j]) * ctmult;
        }
      }
    }
  }

  for (i in 1:T) {
    if (Ct_pred_out[i] > 0) {//negative_infinity() ==
        Ct_pred_out[i] = ctbase - log2(Ct_pred_out[i]/(ctmult*sum(M[i,]) ));//0;
    }
  }

  total = 0;
  for (i in 1:T) {
    total += abs(Ct[i] - Ct_pred_out[i]); // changed from squared
  }
}
