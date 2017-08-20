grcmle <-
function(counts, scheme, method = c("Poisson", "ZIP"),
  do.plot = T, init.guess = NULL, optimizing.algorithm.index = 2,
  lambda.extend.ratio = 3, conf.level = 0.95){
  method <- match.arg(method)
  if(!all(counts == as.integer(counts))){
    warning("non-integer count found. p-value may not be reliable.")
  }
  if(min(counts) < 5){
    warning("one or more count numbers are less than 5. p-value may not be reliable.")
  }
  #require(nloptr)
  method <- method[1]
  stopifnot(method %in% c("Poisson", "ZIP"))
  if(!is.null(init.guess)){
    if(method == "Poisson") stopifnot(length(init.guess) == 1)
    if(method == "ZIP") stopifnot(length(init.guess) == 2)
  }
  ## COMPILING SCHEME
  scheme <- sort(as.integer(scheme[is.finite(scheme)]), decreasing = F)
  scheme.cpld <- lapply(1:length(scheme), function(i){
    if (i == length(scheme)) return(scheme[i] - 1L)
    return(scheme[i]:(scheme[i+1] - 1))
  })
  scheme <- c(scheme, Inf)
  ## MAKING FUNCTION TO OPTIMIZE
  num.groups <- length(scheme.cpld)
  if(num.groups <= 1) stop("insufficient groups")
  if((num.groups <= 2) & (method == "ZIP")) stop("insufficient groups")
  stopifnot(length(counts) == num.groups)
  optim.algorithm <- c(
    "NLOPT_GN_DIRECT_L", "NLOPT_GN_DIRECT",
    "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL",
    "NLOPT_GN_DIRECT_L_NOSCAL", "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
    "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L"
  )[optimizing.algorithm.index]
  if(method == "Poisson"){
    fn1 <- function(lambda){
      blocks <- c(mapply(function(sgmt) sum(dpois(sgmt, lambda)),
        scheme.cpld[-num.groups]), ppois(scheme.cpld[[num.groups]], lambda, F))
      output <- counts * log(blocks)
      output[!is.finite(output)] <- 0
      return(sum(output))
    }
    lower <- 0
    upper <- scheme[num.groups] * lambda.extend.ratio
    ## FINDING INITIAL GUESS
    x.test <- c(init.guess, seq(lower + 1e-5, upper, length.out = 100))
    y.test <- mapply(fn1, x.test)
    plot.x <- x.test; plot.y <- y.test
    x.test <- x.test[which(y.test == max(y.test))[1]]; rm(y.test)
    ## DO OPTIMIZATION
    x.mle <- nloptr(x0 = x.test, eval_f = function(x)-fn1(x), lb = lower + 1e-5,
      ub = upper, opts = list(algorithm = optim.algorithm,
      xtol_rel = 1e-5, maxeval = 3000))$solution
    if(do.plot){
      plot(plot.x, plot.y, type = "l", lwd = 1.5, col = "black",
        xlab="lambda", ylab="log likelihood")
      lines(c(x.mle, x.mle), range(plot.y), col="red", lwd=2)
    }
    blocks <- c(mapply(function(sgmt) sum(dpois(sgmt, x.mle)),
      scheme.cpld[-num.groups]), ppois(scheme.cpld[[num.groups]], x.mle, F))
    Q <- sum((counts - blocks * sum(counts))^2 / (blocks * sum(counts)))
    df <- num.groups - 2
    p.value <- 1 - pchisq(Q, df = df)
    ## Fisher information
    Fisher.info.Poisson <- sum(diff(dpois(scheme - 1, lambda = x.mle)) ^ 2 /
      diff(-ppois(scheme - 0.5, x.mle, F)))
    CI.lambda <- qnorm(p=0.5 + c(-conf.level, conf.level) / 2, mean = x.mle,
      sd = 1 / sqrt(sum(counts) * Fisher.info.Poisson), lower.tail = TRUE)
    CI.p <- NULL
    std.err <- 1 / sqrt(sum(counts) * Fisher.info.Poisson)
  }else{
    fn <- function(p, lambda){
      blocks <- p * c(mapply(function(sgmt) sum(dpois(sgmt, lambda)),
        scheme.cpld[-num.groups]), ppois(scheme.cpld[[num.groups]], lambda, F))
      blocks[1] <- 1 - p + blocks[1]
      output <- counts * log(blocks)
      output[!is.finite(output)] <- 0
      return(sum(output))
    }
    lower <- c(0,0)
    upper <- c(1, scheme[num.groups] * lambda.extend.ratio)
    ## FINDING INITIAL GUESS
    x1.test <- c(init.guess[1], seq(lower[1] + 1e-5, upper[1], length.out = 50))
    x2.test <- c(init.guess[2], seq(lower[2] + 1e-5, upper[2], length.out = 50))
    y.test <- outer(x1.test, x2.test, Vectorize(fn))
    x.grid <- x1.test; y.grid <- x2.test; val.grid <- y.test
    x.ind <- which(y.test == max(y.test), arr.ind = T)[1, ]
    ## DO OPTIMIZATION
    x.mle <- nloptr(x0 = c(x1.test[x.ind[1]], x2.test[x.ind[2]]),
      eval_f = function(x)-fn(x[1], x[2]), lb = lower + 1e-5, ub = upper,
        opts = list(algorithm = optim.algorithm,
      xtol_rel = 1e-5, maxeval = 3000))$solution
    if(do.plot){
      image(x.grid, y.grid, val.grid, xlab="p", ylab="lambda")
      lines(x.mle[1], x.mle[2], pch=19,type="p", col="black")
    }
    blocks <- x.mle[1] * c(mapply(function(sgmt) sum(dpois(sgmt, x.mle[2])),
        scheme.cpld[-num.groups]), ppois(scheme.cpld[[num.groups]], x.mle[2], F))
    blocks[1] <- 1 - x.mle[1] + blocks[1]
    Q <- sum((counts - blocks * sum(counts))^2 / (blocks * sum(counts)))
    df <- num.groups - 3
    p.value <- 1 - pchisq(Q, df = df)
    ## Fisher information
    Fisher.info.Poisson <- sum(diff(dpois(scheme - 1, lambda = x.mle[2])) ^ 2 /
      diff(-ppois(scheme - 0.5, x.mle[2], F)))
    p <- x.mle[1]
    theta1 <- ppois(scheme[2] - 0.5, lambda = x.mle[2], lower.tail = TRUE)
    F11 <- p * (p-1) * dpois(scheme[2]-1, x.mle[2])^2 / theta1 / (1-p+p*theta1) +
      p * Fisher.info.Poisson
    F22 <- (1-theta1)/p/(1-p+p*theta1)
    CI.lambda <- qnorm(p=0.5 + c(-conf.level, conf.level) / 2, mean = x.mle[2],
      sd = 1 / sqrt(sum(counts) * F11), lower.tail = TRUE)
    CI.p <- qnorm(p=0.5 + c(-conf.level, conf.level) / 2, mean = x.mle[1],
      sd = 1 / sqrt(sum(counts) * F22), lower.tail = TRUE)
    CI.p <- pmin(1, pmax(0, CI.p))
    std.err <- c(1 / sqrt(sum(counts) * F22), 1 / sqrt(sum(counts) * F11))
  }
  CI.lambda <- pmax(0, CI.lambda)
  return(list(mle = x.mle, p.value = p.value, df = df, CI.lambda = CI.lambda,
    CI.p = CI.p, conf.level = conf.level, std.err = std.err))
}
