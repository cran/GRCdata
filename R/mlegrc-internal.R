.gen.output <-
function(best.scheme, succeed = TRUE){
  return(list(
    best.scheme.compact   = .parse.scheme.compact(best.scheme),
    best.scheme.loose     = .parse.scheme.loose(best.scheme),
    best.scheme.innerCode = best.scheme,
    succeed               = succeed))
}
.parse.scheme <-
function(scheme){
  return(unlist(lapply(1:(length(scheme) - 1), function(i){
    if (scheme[i] == scheme[i + 1] - 1){
      output <- scheme[i]
    } else {
      output <- paste(scheme[i], "~", scheme[i + 1] - 1, sep = "")
    }
    return(output)
  })))
}
.parse.scheme.compact <-
function(scheme){
  return(sub("~Inf","+",paste(.parse.scheme(scheme), collapse="/", sep = "")))
}
.parse.scheme.loose <-
function(scheme){
  return(paste("{", .parse.scheme(scheme), "}", sep = ""))
}
scoringFun <-
function(model = c("A", "D", "E"), batch = FALSE){
  model <- match.arg(model)
  if(!batch){
    if("A" == model) {
      function(mat)det(mat) / sum(diag(mat))
    } else if ("D" == model) {
      function(mat)det(mat)
    } else if ("E" == model) {
      function(mat)min(eigen(mat, symmetric = T, only.values = T)$values)
    }
  } else {
    if("A" == model) {
      function(m) (m[1,] * m[3,] - m[2,]^2) / (m[1,] + m[3,])
    } else if ("D" == model) {
      function(m) (m[1,] * m[3,] - m[2,]^2)
    } else if ("E" == model) {
      function(m) {
        m.det <- (m[1,] * m[3,] - m[2,]^2)
	m.trc <- (m[1,] + m[3,])
	0.5 * (m.trc - sqrt(pmax(0, m.trc^2 - 4 * m.det)))
      }
    }
  }
}
prob2averageFctl <-
function(probs, lambdas, ps, model=c("Poisson", "ZIP")){
  p <- length(probs)
  stopifnot(p == length(lambdas))
  model <- match.arg(model)
  if("ZIP" == model){
    stopifnot(p == length(ps))
    function(fn){ # fn(lambda,p)
      return(sum(probs * mapply(fn, lambdas, ps)))
    }
  } else if("Poisson" == model) {
    function(fn){ # fn(lambda)
      return(sum(probs * mapply(fn, lambdas)))
    }
  }
}
density2averageFctl <-
function(densityFn, lambda.lwr, lambda.upr, p.lwr, p.upr,
  model=c("Poisson", "ZIP"), tol=1e-9){
  model <- match.arg(model)
  if("ZIP" == model){
    function(fn){ # fn(lambda,p)
      integrand <- function(x)fn(x[1], x[2]) * densityFn(x[1], x[2])
      integral <- adaptIntegrate(integrand, lowerLimit=c(lambda.lwr, p.lwr),
        upperLimit=c(lambda.upr, p.upr), tol=tol, doChecking = T)
      if(0 != integral$returnCode){
        warning(paste("numerical integration warning code:", integral$returnCode))
      }
      return(integral$integral)
    }
  } else if("Poisson" == model) {
    function(fn){ # fn(lambda)
      integrand <- Vectorize(function(x)fn(x) * densityFn(x))
      integral <- integrate(integrand, lower=lambda.lwr, upper=lambda.upr,
        rel.tol=tol)
      if("OK" != integral$message){
        warning(paste("numerical integration warning:", integral$message))
      }
      return(integral$value)
    }
  }
}
find.scheme.internal <-
function(M, N, matSc = c("A", "D", "E"),
  densityFUN, lambda.lwr, lambda.upr, p.lwr, p.upr,
  probs, lambdas, ps,
  is.0.isolated=TRUE, model = c("Poisson", "ZIP"))
{
  # DECIDE PROBABILITY TYPE
  averageFctl <- if(missing(probs)){
    density2averageFctl(densityFn = densityFUN,
      lambda.lwr = lambda.lwr, lambda.upr = lambda.upr,
      p.lwr = p.lwr, p.upr = p.upr, model = model)
  } else {
    prob2averageFctl(probs = probs, lambdas = lambdas, ps = ps, model = model)
  }

  # TRIVIAL CASES
  if("Poisson" == model){
    if (N <= 1L) stop("Number of groups less than 2 for Poisson model. Invalid input.")
  } else {
    if(N <= 2L) stop("Number of groups less than 3 for ZIP model. Invalid input.")
  }
  if(M <= N - 2L) {warning("M too small, changed to N-1"); M <- N-1}
  if(is.0.isolated) if(N == 2) return(.gen.output(c(0, 1, Inf)))

  # MAKING R MATRIX
  R <- matrix(NaN, nrow = M+1L, ncol = M+2L)
  R[1,M+2] <- 0
  for(i in 1:(M+1)) for(j in i:(M+1)) {
    R[i,j] <- averageFctl(function(lambda, p){
      diff(dpois(c(i-2L, j-1L), lambda))^2 / sum(dpois((i-1L):(j-1L), lambda))
    })
  }
  for(i in 2:(M+1)) {
    R[i, M+2] <- averageFctl(function(lambda, p){
      temp <- dpois(i-2L, lambda)^2
      if(temp == 0) return(0)
      temp / ppois(i-2L, lambda, lower.tail = F)
    })
  }

  # MAKING W
  W <- averageFctl(function(lambda, p){
    a <- dpois((-1L):M, lambda)
    return(1/lambda - sum(diff(a)^2 / a[-1]))
  })

  # MAKING U MATRIX
  if("ZIP" == model){
    M.0 <- if(is.0.isolated) 1L else M-N+3L
    U <- matrix(0, nrow = 3L, ncol = M-N+3L)
    for(i in 1:M.0){
      U[1,i] <- averageFctl(function(lambda, p){
        mu <- ppois(q = i - 0.5, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
        mu.prime <- -dpois(x = i-1L, lambda = lambda, log = FALSE)
        return(p * (p-1) * (mu.prime^2) / (mu * (1 - p + p * mu)))
      })
      U[2,i] <- averageFctl(function(lambda, p){
        mu <- ppois(q = i - 0.5, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
        mu.prime <- -dpois(x = i-1, lambda = lambda, log = FALSE)
        return( -mu.prime / (1 - p + p * mu))
      })
      U[3,i] <- averageFctl(function(lambda, p){
        mu <- ppois(q = i - 0.5, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
        return((1 - mu) / (p * (1 - p + p * mu)))
      })
    }
    p0 <- averageFctl(function(lambda, p)p)
  }

  # MAKING SEARCHING COMBINATIONS; N index M+3 := Inf SINCE "-1" LATER WHEN USED
  scheme.search <- if(is.0.isolated){
    rbind(1L, 2L, if(M == 2) 3 else combn(x = 3:(M+1), m = N - 2L), M + 3L)
  }else{
    rbind(1L, if(M == 1) 2 else combn(x = 2:(M+1), m = N - 1L), M + 3L)
  } # EACH *COLUMN* IS A CODE VECTOR

  # BEGIN: SEARCHING
  performance <- double(ncol(scheme.search))
  for(i in 1:(nrow(scheme.search) - 1)) {
    performance <- performance + R[cbind(scheme.search[i, ], scheme.search[i+1, ] - 1)]
  }
  if("ZIP" == model) { # MODIFYING PERFORMANCE SCORE ACCORDINGLY
    J <- matrix(0, nrow = 3, ncol = ncol(scheme.search))
    J[1, ] <- U[1, scheme.search[2, ] - 1] + p0 * performance
    J[2, ] <- U[2, scheme.search[2, ] - 1]
    J[3, ] <- U[3, scheme.search[2, ] - 1]
    performance <- scoringFun(model = matSc, batch = TRUE)(J)
    rm(J)
  }
  best.perf <- max(performance)
  best.indx <- which(performance == best.perf)
  if(length(best.indx) > 1){
    warning(paste("more than one best index found:", paste(best.indx, collapse = " ")))
    best.indx <- best.indx[1]
  }
  #write(x = paste(performance), file = "find.schemeSh.txt", ncolumns = 1)
  rm(performance)


  # MAKING VALIDATING COMBINATIONS
  scheme.validate<- if(is.0.isolated){
    rbind(1L, 2L, if(M == 2) 3 else combn(x = 3:(M+1), m = N - 3L), M + 2L)
  }else{
    rbind(1L, if(M == 1) 2 else combn(x = 2:(M+1), m = N - 2L), M + 2L)
  }

  # BEGIN: VALIDATING
  performance <- rep(W, ncol(scheme.validate))
  for(i in 1:(nrow(scheme.validate) - 1)) {
    performance <- performance + R[cbind(scheme.validate[i, ], scheme.validate[i+1, ] - 1)]
  }
  if("ZIP" == model) {
    J <- matrix(0, nrow = 3, ncol = ncol(scheme.validate))
    J[1, ] <- U[1, scheme.validate[2, ] - 1] + p0 * performance
    J[2, ] <- U[2, scheme.validate[2, ] - 1]
    J[3, ] <- U[3, scheme.validate[2, ] - 1]
    performance <- scoringFun(model = matSc, batch = TRUE)(J)
  }
  #browser()
  succeed <- (max(performance) <= best.perf)
  #write(x = paste(performance), file = "find.schemeVr.txt", ncolumns = 1)

  # GENERATE OUTPUT
  best.scheme <- scheme.search[, best.indx, drop = TRUE] - 1
  best.scheme[length(best.scheme)] <- Inf
  return(.gen.output(best.scheme = best.scheme, succeed = succeed))
}
