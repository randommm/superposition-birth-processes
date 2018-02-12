long_time_sample <- function(lambda, mu, rho, p) {

  yule_start <- function(lambda_, d) {
    self = list()
    self$lambda_ = lambda_
    self$upsilon = rexp(d, lambda_)
    return(self)
  }

  yule_get_next <- function(self) {
    current_min_index = which.min(self$upsilon)
    current_min = self$upsilon[current_min_index]

    self$upsilon[current_min_index] = current_min + rexp(1, self$lambda_)
    self$upsilon[length(self$upsilon) + 1] = current_min + rexp(1, self$lambda_)

    self$current_min = current_min
    return(self)
  }

  m <- rpois(1, rho)
  d <- rbinom(1, m, p)

  if (d == 0) {
    return(Inf)
  }

  v <- d

  process <- yule_start(lambda, d)
  new_t <- 0
  y <- rexp(d, mu)

  while (TRUE) {
    process <- yule_get_next(process)
    new_t <- process$current_min

    #A previous cell got promoted
    if (min(y) < new_t)
      return(min(y))

    new_y <- rexp(1, mu)

    #This cell got promoted
    if (new_y < new_t)
      return(new_t)

    y <- c(y, new_y)
    v <- v+1
  }
}

require(rstan)

stan_code1 <- "
data {
  int n_obs;
  int n_cen;
  real t_obs[n_obs];
  real t_cen[n_cen];
  real x1_obs[n_obs];
  real x1_cen[n_cen];
}
parameters {
  real<lower=0> lambda;
  real<upper=0> beta0;
  real beta1;
  real<lower=0> mu;
}
model {
  for (i in 1:n_cen) {
    real t;
    t = t_cen[i];

    target +=
      (exp(mu*t) - 1) *
      exp(beta0 + beta1 * x1_cen[i] + lambda*t) /
      (exp(lambda*t) - exp(t*(lambda + mu)) - 1);
  }

  for (i in 1:n_obs) {
    real log_theta_obs;
    real t;

    log_theta_obs = beta0 + beta1 * x1_obs[i];
    t = t_obs[i];

    target +=
      (exp(mu*t) - 1) *
      exp(log_theta_obs + lambda*t) /
      (exp(lambda*t) - exp(t*(lambda + mu)) - 1);

    target +=
      log(-lambda + mu*(exp(lambda*t) - 1) + (lambda + mu)*
      (-exp(lambda*t) + exp(t*(lambda + mu)) + 1));

    target +=
      log_theta_obs;

    target +=
      - 2 * log(-exp(lambda*t) + exp(t*(lambda + mu)) + 1);
  }

  lambda ~ cauchy(0, 2.5);
  mu ~ cauchy(0, 2.5);
}
generated quantities {
  real p;
  p = exp(beta0);
}"

lambda <- 0.7
mu <- 1 / 4.5
p <- 0.7

beta <- 0.45
n <- 10000
x <- rnorm(n, -1, 2)

t <- Vectorize(long_time_sample)(lambda, mu, exp(beta * x), p)

cen_index <- (t==Inf)
obs_index <- (t!=Inf)

n_obs <- sum(obs_index)
n_cen <- sum(cen_index)
t_obs <- t[obs_index]
t_cen <- t[cen_index]
x1_obs <- x[obs_index]
x1_cen <- x[cen_index]

t_cen <- abs(rnorm(n_cen, 0, 25))

compiled_model1 <- rstan::stan_model(model_code=stan_code1)

print(rstan::optimizing(compiled_model1))

fit <- rstan::sampling(compiled_model1, iter=1e4, chains=2,
                       init_r=.1, cores=3,
                       control=list(adapt_delta=0.99),
                       refresh=100)
