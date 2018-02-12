# Primeiramente definimos a função geradora de variaveis aleatorias no R:

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

    # A previous cell got promoted
    if (min(y) < new_t)
      return(min(y))

    new_y <- rexp(1, mu)

    # This cell got promoted
    if (new_y < new_t)
      return(new_t)

    y <- c(y, new_y)
    v <- v+1
  }
}

# Em seguida escolhemos os valores dos parametros:

lambda <- c(0.5, 1.0, 1.5, 2.0)
mu <- 0.4
p <- 0.7
rho <- 1.25

# lambda <- c(1.0, 3.0, 5.0, 7.0, 10, 15)
# mu <- 0.4
# p <- 0.25
# rho <- 3

observations <- observations_cured <- observations_not_cured <- list()
for (i in 1:length(lambda)) {
  # Feito isso, geremos algumas observações:
  observations[[i]] <- replicate(10000, long_time_sample(lambda[i], mu, rho, p))

  # Agora vamos calcular quantos pacientes foram curados
  # (ou seja, sobreviveram até o infinito).
  observations_cured[[i]] <- observations[[i]][observations[[i]] == Inf]
  observations_not_cured[[i]] <- observations[[i]][observations[[i]] < Inf]

  # Com isso podemos calcular a taxa de cura empirica:
  print("Taxa de cura empirica:")
  print(
  length(observations_cured[[i]]) /
  (length(observations_cured[[i]]) + length(observations_not_cured[[i]]))
  )

  print("Taxa de cura verdadeira:")
  print(dpois(0, rho*p))
}

colors <- c("blue", "red", "forestgreen", "purple", "yellow", "black")

# Agora vamos plotar a cumulativa empirica para os pacientes em risco
pdf("ecdf.pdf",width=7,height=5)
plot(ecdf(observations_not_cured[[1]]), lty=1, lwd=1,
  col=colors[1])
for (i in 2:length(lambda))
  lines(ecdf(observations_not_cured[[i]]), lty=i, lwd=i,
  col=colors[i])

legend("bottomright", legend=lambda, cex=0.8, col=colors, inset=.1,
      lty=1:4, lwd=1:4, bty="n")
dev.off()
#Por fim, a sobrevivencia empirica (i.e.: Kaplan-Meier)
require(survival)
pdf("surv_sim_1.pdf",width=7,height=5)
plot(survfit(Surv(observations_not_cured[[1]])~1, conf.type="none"), lty=1, lwd=2.5,
     col=colors[1], xlim=c(0, 3.0), ylab="Survival function")
for (i in 2:length(lambda))
  lines(survfit(Surv(observations_not_cured[[i]])~1, conf.type="none"),
  lty=i, lwd=2.5, col=colors[i])

legend("topright", legend=lambda, cex=0.8, col=colors, inset=.1,
       lty=1:length(lambda), lwd=2.5, bty="n")
dev.off()
