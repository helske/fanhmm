data {
  int<lower=1> N; // number of individuals
  int<lower=2> T; // number of time points
  int<lower=2> M; // number of observed symbols
  int<lower=1> S; // number of hidden states
  array[T, N] int<lower=0, upper=M> obs; // observations
  real<lower=0> lambda;
  matrix[S, S - 1] Qs;
  matrix[M, M - 1] Qm;
}
parameters {
  // parameters for the initial state probability vector
  vector[S - 1] eta_pi;
  // parameters for the transition probabilities
  matrix[S - 1, S] eta_A;
  // parameters for the emissions probabilities
  matrix[M - 1, S] eta_B;
}
transformed parameters {
  vector[S] Pi = softmax(Qs * eta_pi);
  matrix[S, S] A;
  matrix[S, M] log_B;
  for (s in 1:S) {
    A[s, ] = softmax(Qs * eta_A[, s])';
    log_B[s, ] = log_softmax(Qm * eta_B[, s])';
  }
  real log_lik = 0;
  {
    matrix[S, T] log_omega;
    for (i in 1:N) {
      for (t in 1:T) {
        log_omega[, t] = log_B[, obs[t, i]];
      }
      log_lik += hmm_marginal(log_omega, A, Pi);
    }
  }
}
model {
  // priors matching the regularization in ML case
  target += -0.5 * lambda * sum(square(eta_pi));
  target += -0.5 * lambda * sum(square(to_vector(eta_A)));
  target += -0.5 * lambda * sum(square(to_vector(eta_B)));
  target += log_lik;
}

generated quantities {
  real max_A = max(A); // invariant to label-switching
}
