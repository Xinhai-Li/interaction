% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FPIR1twoway.R
\name{FPIR1twoway}
\alias{FPIR1twoway}
\title{Fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M *x2^N}
\usage{
FPIR1twoway(y, x1, x2)
}
\arguments{
\item{y}{the dependent variable in regression y ~ x1 + x2 + x1^M *x2^N}

\item{x1}{one of the independent variables in regression y ~ x1 + x2 + x1^M *x2^N}

\item{x2}{another independent variable in regression y ~ x1 + x2 + x1^M *x2^N}
}
\value{
return a data.frame with columns of M, N, R square value of the whole model,
proportion of the variance explained by the interaction term, coefficients of all terms
(the interception, x1, x2, and the interaction).
}
\description{
This function estimates the optimum values of exponents (M and N) in FPIR,
 return model parameters such as regression coefficients, sum of squares of the interaction term, R square, AIC.
}
\examples{

 x1 = runif(100); x2 = runif(100);
 y = x1 * x2^10; y = y + rnorm(length(y), 0, 0.01) # generate a dataset
 results = FPIR1twoway (y, x1, x2)
 results[[1]] # parameters for 10000 models
 results[[2]] # parameters for the best model
 results[[3]] # regression coefficients
 results[[4]] # adjusted R square

}
\author{
Xinhai Li (Xinhai_li_edu@126.com)
}
