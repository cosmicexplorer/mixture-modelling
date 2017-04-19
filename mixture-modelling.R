## X is n x d
## mu is k x d
## sd is k x d
## alpha is vector(k)
## pdf is a function accepting (n x d, vector(d), vector(d)) and returning n x k
calc_expectation <- function(X, mu, sd, alpha, pdf, like_fun = log) {
    ## TODO: have configurable pdf, like_fun with good defaults (and a
    ## standardized set of arguments they'll accept)
    data_dim <- dim(X)
    n <- data_dim[1]
    d <- data_dim[2]
    k <- dim(mu)[1]
    ## TODO: allow deviations along:
    ## 1. not just the coordinate axes
    ## 2. fewer than d axes
    stopifnot(
        all(dim(mu) == c(k, d) &
            dim(sd) == c(k, d) &
            length(alpha) == k))
    stopifnot(isTRUE(all.equal(sum(alpha), 1)))

    ## n x k
    prob <- alpha * vapply(1:k, function (i) pdf(X, mu[i,], sd[i,]),
                           vector("double", n))

    ## vector(n)
    sum_prob <- rowSums(prob)
    ## n x k
    all_post <- prob / sum_prob

    ## vector(n)
    transformed_likelihood <- like_fun(sum_prob)
    ## scalar
    sum_transformed_likelihood <- sum(transformed_likelihood)

    list(posterior = all_post,          # n x k
         likelihood_tr = sum_transformed_likelihood) # scalar
}

## X is n x d
## single_mu is vector(d)
## returns n x d
euclidean <- function(X, single_mu) {
    stopifnot(dim(X)[2] == length(single_mu))
    ## subtract from each data point
    ## same as t(t(x) - single_mu)^2
    sweep(X, 2, single_mu)^2
}

## X is n x d
## posterior is n x k
## distance is a function of (n x d, vector(d)) -> n x d
maximize_expectation <- function(X, posterior, distance = euclidean) {
    data_dim <- dim(X)
    n <- data_dim[1]
    d <- data_dim[2]
    post_dim <- dim(posterior)
    k <- post_dim[2]
    stopifnot(post_dim[1] == n)

    ## vector(k)
    sum_post <- colSums(posterior)

    ## k x d
    unscaled_mu <- t(posterior) %*% X
    ## k x d
    mu <- sum_post^(-1) * unscaled_mu

    ## k x d
    var <- sum_post^(-1) * vapply(1:k, function (i) {
        ## n x d
        scaled <- posterior[,i] * distance(X, mu[i,])
        ## 1 x d
        colSum(scaled)
    }, vector("double", d))

    alpha

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
