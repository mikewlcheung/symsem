[![Build Status](https://travis-ci.org/mikewlcheung/symsem.svg?branch=master)](https://travis-ci.com/github/mikewlcheung/symsem)
[![cran version](http://www.r-pkg.org/badges/version/symSEM)](https://cran.r-project.org/package=symSEM)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/symSEM)](http://cranlogs.r-pkg.org/badges/symSEM)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/symSEM)](http://cranlogs.r-pkg.org/badges/grand-total/symSEM)
[![Rdoc](http://www.rdocumentation.org/badges/version/symSEM)](http://www.rdocumentation.org/packages/symSEM)

# symSEM
Symbolic Computation for Structural Equation Models

# Example

## Compute a symbolic model-implied covariance and correlation matrices

```
library(symSEM)

## A regression model
model <- "y ~ b1*x1 + b2*x2
          ## Covariance between x1 and x2
          x1 ~~ x2
          ## Means
          y ~ b0*1
          x1 ~ m1*1
          x2 ~ m2*1"

## Convert it into a RAM speculation
RAM <- metaSEM::lavaan2RAM(model)

## Implied covariance matrix and mean structure
impliedS(RAM, corr=FALSE)

## Implied correlation matrix
impliedS(RAM, corr=TRUE)
```
