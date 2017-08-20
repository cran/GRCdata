\name{GRCdata-package}
\alias{GRCdata-package}
\alias{GRCdata}
\docType{package}
\title{
Parameter inference and optimal designs for grouped and/or right-censored count data
}
\description{
This package consists of two main functions:
The first function uses a given grouped and/or right-censored grouping scheme and empirical data to infer parameters, and implements chi-square goodness-of-fit tests;
The second function searches for the global optimal grouping scheme of grouped and/or right-censored count responses in surveys.

This R package is designed to implement methods and algorithms developed in the following papers and please cite these articles at your convenience:

  Qiang Fu, Xin Guo and Kenneth C. Land. Forthcoming. "A Poisson-Multinomial Mixture Approach to Grouped and Right-Censored Counts." Communications in Statistics -- Theory and Methods.
  DOI: 10.1080/03610926.2017.1303736
  (mainly about the first function for aggregate-level parameter inference)
  
  Qiang Fu, Xin Guo and Kenneth C. Land. Conditionally accepted. "Optimizing Count Responses in Surveys: A Machine-Learning Approach." Sociological Methods & Research. 
  (mainly about the second function for finding optimal grouping schemes)

To install the package "GRCdata_1.0.tar.gz", one may place this
file in the working directory/folder of R, and type

\code{install.packages("GRCdata", repos = NULL, type = "source")}

To check the current working directory of R, one may type

\code{getwd()}

To see the source code, one could extract the package ``GRCdata_1.0.tar.gz''.
There would be two directories/folders: \code{man} and \code{R}. The source code
is under the \code{R} directory/folder.
}
\details{
\tabular{ll}{
Package: \tab GRCdata\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab July 28, 2017\cr
License: \tab GPLv3\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Authors: Xin Guo <x.guo@polyu.edu.hk>, Qiang Fu <qiang.fu@ubc.ca>

Maintainers: Xin Guo <x.guo@polyu.edu.hk>
}
\references{
  Qiang Fu, Xin Guo and Kenneth C. Land. Conditionally accepted. "Optimizing Count Responses in Surveys: A Machine-Learning Approach." Sociological Methods & Research. 

  Qiang Fu, Xin Guo and Kenneth C. Land. Forthcoming. "A Poisson-Multinomial Mixture Approach to Grouped and Right-Censored Counts." Communications in Statistics -- Theory and Methods.
  DOI: 10.1080/03610926.2017.1303736  
}
