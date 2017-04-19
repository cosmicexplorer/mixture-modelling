vt_tmpl <- vector("double", 1)

calc_expectation <- function(X, mu, sd, alpha, pdf, like_fun = log) {
    ## TODO: have configurable pdf, like_fun with good defaults (and a
    ## standardized set of arguments they'll accept)
    data_size <- dim(X)
    n <- data_size[1]
    d <- data_size[2]
    k <- dim(mu)[1]
    ## TODO: allow deviations along:
    ## 1. not just the coordinate axes
    ## 2. fewer than d axes
    stopifnot(
        all(dim(mu) == c(k, d) &
            dim(sd) == c(k, d) &
            length(alpha) == k))
    stopifnot(isTRUE(all.equal(sum(alpha), 1)))

    prob <- vapply(1:k, function (i) pdf(X, mu[i,], sd[i, ]), vt_tmpl) %*% alpha

    sum_prob <- sum(prob)
    all_post <- prob / sum_prob

    transformed_likelihood <- like_fun(sum_prob)
    sum_transformed_likelihood <- sum(transformed_likelihood)

    list(post = all_post,
         likelihood_tr = sum_transformed_likelihood)
}

maximize_expectation <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])

  comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
  comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)

  comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
  comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n

  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)

  list("mu" = c(comp1.mu, comp2.mu),
       "var" = c(comp1.var, comp2.var),
       "alpha" = c(comp1.alpha, comp2.alpha))
}
