find.scheme <-
function(N,
  densityFUN, lambda.lwr, lambda.upr, p.lwr, p.upr,
  probs, lambdas, ps,
  is.0.isolated=TRUE, model = c("Poisson", "ZIP"),
  matSc = c("A", "D", "E"), M = "auto")
{
  # ASSUME: CONDITIONS ALREADY VALIDATED
  N <- as.integer(N)
  model <- match.arg(model)
  matSc <- match.arg(matSc)

  # DECIDE PROBABILITY TYPE
  temp <- missing(densityFUN) + missing(probs)
  if(temp == 0){
    stop("please specify 'densityFUN' or 'probs'")
  } else if(temp == 2) {
    stop("specify 'densityFUN' or 'probs' but NOT BOTH")
  }

  # DISCRETE DISTRIBUTIONS
  if(missing(densityFUN)) {
    stopifnot(length(probs) == length(lambdas))
    probs <- probs / sum(probs)
  }

  # AUTO M MODE - CHECKPOINT
  if(M == "auto") {
    not.Finished <- TRUE
    # begin to ROUGHLY estimate the 95% quantile to make the starting point
    P.BYDESIGN <- 0.95
    if(missing(probs)){
      quantile95 <- P.BYDESIGN * lambda.upr + (1 - P.BYDESIGN) * lambda.lwr
    } else {
      ord.lambdas <- order(lambdas, decreasing = TRUE)
      lambda.ind <- ord.lambdas[which(cumsum(probs[ord.lambdas]) > (1 - P.BYDESIGN))[1]]
      quantile95 <- lambdas[lambda.ind]
      rm(ord.lambdas, lambda.ind)
    }
    M.try <- max(as.integer((quantile95 + (N/2) + 0.5)), N - 1)
    times.tried <- 1
    # BEGIN TESTING DIFFERENT M'S
    while(not.Finished){
      # if it really takes too long to run, the user needs to stop the
      # code manually
      output <- find.scheme.internal(M = M.try,
        N = N, matSc = matSc, densityFUN = densityFUN, lambda.lwr = lambda.lwr,
        lambda.upr = lambda.upr, p.lwr = p.lwr, p.upr = p.upr,
        probs = probs, lambdas = lambdas, ps = ps,
        is.0.isolated = is.0.isolated, model = model)
      if(output$succeed) {
        not.Finished <- FALSE
      } else {
        M.try <- M.try + 2
        times.tried <- times.tried + 1
      }
    }
    output$M <- M.try
    output$times.tried <- times.tried
    return(output)
  } else {
    M <- as.integer(M)
    output <- find.scheme.internal(M = M,
      N = N, matSc = matSc, densityFUN = densityFUN, lambda.lwr = lambda.lwr,
      lambda.upr = lambda.upr, p.lwr = p.lwr, p.upr = p.upr,
      probs = probs, lambdas = lambdas, ps = ps,
      is.0.isolated = is.0.isolated, model = model)
    output$M <- M
    output$times.tried <- 1
    return(output)
  }
  return(NULL)
}
