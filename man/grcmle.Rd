\name{grcmle}
\alias{grcmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Maximum likelihood estimation of Poisson or ZIP parameters at the aggregate level.
}
\description{
This function infers Poisson or zero-inflated Poisson (ZIP) parameters from grouped and right-censored count data,
and conducts a chi-squared goodness-of-fit test.
A grouped and right-censored scheme may look like

\code{  0, 1, 2--4, 5--8, 9+.}

For grouped and right-censored count data collected in a survey, such as frequency of alcohol drinking, number of births or occurrence of crimes,
the response category designed as the example above means never, once, 2 to 4 times, 5 to 8 times, 9 times and more.
The frequency distribution from a sample corresponding to the example above may look like

\code{  3, 15, 168, 155, 15}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grcmle(counts, scheme, method = c("Poisson", "ZIP"),
  do.plot = T, init.guess = NULL,
  optimizing.algorithm.index = 2, lambda.extend.ratio = 3,
  conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
specifies the frequency distribution of the grouped and right-censored count data. For the example above, one may input

\code{  counts = c(3, 15, 168, 155, 15)}.
}
  \item{scheme}{
specifies the grouping scheme. It should be a vector of integers containing the starting
point (or the lowest integer) of each group. For example, to input the scheme above

\code{  0, 1, 2--4, 5--8, 9+,}

one may use

\code{  scheme = c(0, 1, 2, 5, 9)}.
}
  \item{method}{
a string parameter specifies which statistical model to use. Currently there are two options
\code{"Poisson"} and \code{"ZIP"}. The default value is \code{"Poisson"}.
It can be abbreviated.
%%     ~~Describe \code{method} here~~
}
  \item{do.plot}{
a logical variable indicating whether or not to plot the log likelihood.
The default is \code{T}.
}
  \item{init.guess}{
the initial value used for the optimization procedure of the likelihood estimation. The default value
is \code{NULL}, which instructs the function \code{grcmle} to select the initial value automatically.
}
  \item{optimizing.algorithm.index}{
%%     ~~Describe \code{optimizing.algorithm.index} here~~
defines which optimization algorithm to use. Currently the possible values are \code{1,2,3,4,5,6,7} and \code{8},
representing the following algorithms, respectively:

\code{  NLOPT_GN_DIRECT_L}

\code{  NLOPT_GN_DIRECT}

\code{  NLOPT_GN_DIRECT_L_RAND}

\code{  NLOPT_GN_DIRECT_NOSCAL}

\code{  NLOPT_GN_DIRECT_L_NOSCAL}

\code{  NLOPT_GN_DIRECT_L_RAND_NOSCAL}

\code{  NLOPT_GN_ORIG_DIRECT}

\code{  NLOPT_GN_ORIG_DIRECT_L}

For details of these algorithms, please see the manual of the \R package \code{"nloptr"}.
The default value is 2.
}
  \item{lambda.extend.ratio}{
specifies the searching interval of possible \eqn{\lambda} as \eqn{[0, nr]}, where \eqn{n} is the left end (i.e., the lowest integer)
of the last right-censored group, and \eqn{r} is \code{lambda.extend.ratio}. By default,
we set

\code{lambda.extend.ratio=3}.
}
  \item{conf.level}{
    confidence level of the confidence interval(s) for the parameter(s) inferred
  }
}
\details{
Maximum likelihood estimation is used for the inference.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The returned value is a list containing
\item{mle}{the parameter(s) inferred. For Poisson model, it is the estimate of \eqn{\lambda}.
For ZIP model, it shows a vector of length 2: the first element is the estimate of \eqn{p} and
the second element is the estimate of \eqn{\lambda}.}
\item{p.value}{the p-value of the chi-squared test of goodness-of-fit.}
\item{df}{the degree(s) of freedom of the chi-squared test of goodness-of-fit.}
\item{CI.lambda}{the confidence interval of \eqn{\lambda} obtained by normal approximation}
\item{CI.p}{the confidence interval of \eqn{p} obtained by normal approximation}
\item{conf.level}{the confidence level}
\item{std.err}{the standard error of \eqn{\lambda}
  or the standard errors of \eqn{(p, \lambda)}, if a ZIP model is specified}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% ~put references to the literature/web site here ~
%}
\author{
%%  ~~who you are~~
Authors: Xin Guo <x.guo@polyu.edu.hk>, Qiang Fu <qiang.fu@ubc.ca>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
grcmle(counts=c(6, 15, 168, 155, 15), scheme = c(0, 1, 2, 5, 9))
\donttest{
grcmle(counts=c(6, 15, 168, 155, 15), scheme = c(0, 1, 2, 5, 9), method = "ZIP")
}
}
% \keyword{ Maximum Likelihood Estimate }
%  \keyword{ Poisson }
%  \keyword{ ZIP }% __ONLY ONE__ keyword per line
