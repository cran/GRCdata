\name{find.scheme}
\alias{find.scheme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Search for the global optimal grouping scheme of grouped and/or right-censored count data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given the prior distribution (or values) of parameters,
and the total/maximum number of groups (N) allowed for grouping schemes,
this function finds the global optimal grouping scheme that makes
the sampling process most informative.

}
\usage{
find.scheme(N,
  densityFUN, lambda.lwr, lambda.upr, p.lwr, p.upr,
  probs, lambdas, ps,
  is.0.isolated = TRUE, model = c("Poisson", "ZIP"),
  matSc = c("A", "D", "E"), M = "auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
(maximum) number of groups allowed for all grouping schemes. A non-integral value will be coerced to an integer.
}
  \item{densityFUN, lambda.lwr, lambda.upr, p.lwr, p.upr}{
prior information of parameters in a continuous form.
These parameters denote the prior probability density function (optional),
the lower bound of \eqn{\lambda} (for Poisson models), the higher bound of \eqn{\lambda} (for Poisson models),
the lower bound of \code{p} (optional for ZIP models),
the higher bound of \code{p} (optional for ZIP models), respectively.
}
  \item{probs, lambdas, ps}{
prior information of the parameters in a discrete form.
These parameters are vectors denoting the mass probabilities, the corresponding values of \eqn{\lambda} and \code{p} (optional), respectively.
}
  \item{is.0.isolated}{
%%     ~~Describe \code{is.0.isolated} here~~
a logical value indicating whether zero is contained and only contained in a single group.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
underlying Poisson models to be used for optimal designs: \code{Poisson} or \code{ZIP}.
The default value is \code{Poisson}.
}
  \item{matSc}{
%%     ~~Describe \code{matSc} here~~
A character indicating types of optimality functions of the Fisher information (matrix).
It must be one from the three letters: \code{A}, \code{D}, and \code{E}.
In particular,
if \code{J} is the 2-by-2 Fisher information matrix, then

\code{"A"} or A-optimality maximizes \eqn{1/\mathrm{tr}(J^{-1})};

\code{"D"} or D-optimality maximizes \eqn{\mathrm{det}(J)};

\code{"E"} or E-optimality maximizes the minimum eigenvalue of \eqn{J}.
}
  \item{M}{
%%     ~~Describe \code{M} here~~
a sufficiently large integer needed to facilitate the search, or a character \code{"auto"}.
Theoretically, it could be the lowest integer contained in the last right-censored
group of the global optimal grouping scheme. A non-integral value will be coerced to an integer.
If \code{M} is set to be \code{"auto"}, the algorithm takes longer time to converge
because it will automatically determine \code{M} and return the global optimal
grouping scheme;
The default value of \code{M} is \code{"auto"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function tries to find the N-group scheme maximizing Fisher information (matrix).
If \code{model} is specified as Poisson, \code{p.lwr} or \code{p.upr} will be ignored.
When the prior distribution is discrete, \code{lambdas} specify discrete values that \eqn{\lambda} may take,
and \code{probs} specify probabilities associated with \code{p}.
In the ZIP model, \code{lambdas} and \code{ps} specify discrete values that \eqn{\lambda} and p may take, respectively.
\code{probs} denotes joint mass probabilities associated with (\eqn{\lambda}, p). The values of
(p.lwr, p.upr) cannot be (0, 1) as the algorithm will not converge. Instead, approximate values, such as
(0.000001, 0.999999), can be used.

A sufficiently large integer \code{M} should be provided by the user so that infinitely many grouping schemes
could be handled by the search algorithm. \code{M} is in theory the lowest integer to be contained in the last
right-censored group of the global optimal grouping scheme. In practice, the choice of \code{M} should be slightly higher than
its theoretical value because the search algorithem is designed in a way that it prevents any acceptance of a false optimal
solution at the cost of tolerating false rejection of the correct optimal grouping scheme. This idea is implemented by
a logical indicator \code{succeed} in the output. Its value will be \code{TRUE} if the real optimal grouping scheme is identified.
Otherwise, a \code{FALSE} output means that \code{M} is not large enough to gurantee that the grouping scheme yielded by
the search algorithm is the global optimal grouping scheme. Researchers then need to select a larger \code{M} and repeat this
process until the logical indicator \code{succeed} becomes \code{TRUE}. Alternatively, users may use the \code{"auto"} option so that
this iterative process will be automatically implemented.

}
\value{
The returned value is a list with components.
\item{best.scheme.compact, best.scheme.loose, best.scheme.innerCode}{the same optimal grouping scheme is printed in various forms.}
\item{succeed}{see Details. This is a logical variable. The global optimal grouping scheme is obtained if it is \code{TRUE};
a larger \code{M} needs to be selected for a successful search if it is \code{FALSE}.}
}
\references{
%% ~put references to the literature/web site here ~
Qiang Fu, Xin Guo and Kenneth C. Land. Conditionally accepted. "Optimizing Count Responses in Surveys: A Machine-Learning Approach." Sociological Methods & Research.

}
\author{
%%  ~~who you are~~
Xin Guo <x.guo@polyu.edu.hk>, Qiang Fu <qiang.fu@ubc.ca>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# Example 1 ####################################
# M=7, N=3, 0 is not required to be contained
# in a separate group of grouping schemes.
# Poisson model, lambda takes 4 and 5 and each value has a probability of 0.5.
find.scheme(probs = c(0.5, 0.5), lambdas = c(4,5),
  M = 7, N = 3, is.0.isolated = FALSE, model = "Poisson")

# Example 2 ####################################
# N=3, 0 is required to be contained in a separate group of grouping schemes.
# Poisson model, lambda takes 4 and 5 and each value has a probability of 0.5.
# M is not given, so it will be selected automatically.
find.scheme(probs = c(0.5, 0.5), lambdas = c(4,5),
  N = 3, is.0.isolated = TRUE, model = "Poisson")

# Example 3 ####################################
# M=7, N=3, 0 is not required to be contained in a separate group.
# ZIP model, (lambda, p) take (4, 0.3) and (5, 0.4)
# with their probabilities denoted by c(0.5, 0.5)
\donttest{
find.scheme(probs = c(0.5, 0.5), lambdas = c(4,5), ps = c(0.3, 0.5),
  M = 7, N = 3, is.0.isolated = FALSE, model = "ZIP")
}

# Example 4 ####################################
# N=3, 0 is not required to be contained in a separate group.
# Poisson model, lambda takes a normal distribution truncated to [1, 10]
# M is not given, so it will be selected automatically.
\donttest{
find.scheme(densityFUN = function(lambda)
  dnorm(lambda, mean = 3, sd = 1),
  lambda.lwr = 1, lambda.upr = 10,
  N = 3, is.0.isolated = FALSE, model = "Poisson")
}

# Example 5 ####################################
# M=7, N=3, 0 is required to be contained in a separate group.
# Poisson model, lambda takes a normal distribution truncated to [1, 10]
\donttest{
find.scheme(densityFUN = function(lambda)
  dnorm(lambda, mean = 3, sd = 1),
  lambda.lwr = 1, lambda.upr = 10,
  M = 7, N = 3, is.0.isolated = TRUE, model = "Poisson")
}

# Example 6 ####################################
# N=3, 0 is required to be contained in a separate group.
# Poisson model, lambda takes an uniform distribution on [1, 10]
# M is not given, so it will be selected automatically.
find.scheme(densityFUN = function(lambda)
  dunif(lambda, min = 1, max = 10),
  lambda.lwr = 1, lambda.upr = 10,
  N = 3, is.0.isolated = TRUE, model = "Poisson")

# Example 7 #################################
# M=7, N=3, 0 is required to be contained in a separate group.
# ZIP model, (lambda, p) has an uniform distribution with
# lambda on [1,10] and p on [0.1, 0.9]
\donttest{
find.scheme(densityFUN = function(...) 1,
  lambda.lwr = 1, lambda.upr = 10, p.lwr = 0.0001, p.upr = 0.9999,
  M = 7, N = 3, is.0.isolated = TRUE, model = "ZIP")
}

# Example 8 ####################################
# M=7, N=3, 0 is required to be contained in a separate group.
# ZIP model, (lambda, p) has a normal distribution centered
# at (5.5, 0.5) with a covariance matrix showing their correlation
#    /                 \
#    |  11/3      3    |
#    |    3     11/3   |
#    \                 /.
# This normal distribution is also truncated to
# [1, 10] X [0.1, 0.9]
# Note: this example may take several minutes to converge,
# depending on your computer configuration.
\donttest{
dsty <- function(lambda, p){
  vec <- c(lambda - 5.5, p - 0.5)
  mat <- matrix(c(11/3,3,3,11/3), nrow = 2, ncol = 2)
  pw <- -0.5 * sum(vec * solve(mat, vec))
  return(exp(pw))
}
find.scheme(densityFUN = dsty,
  lambda.lwr = 1, lambda.upr = 10, p.lwr = 0.1, p.upr = 0.9,
  M = 7, N = 3, is.0.isolated = TRUE, model = "ZIP")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
