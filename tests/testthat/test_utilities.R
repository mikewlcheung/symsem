context("Checking utility functions")

test_that("sym() works correctly", {
  ## Test one element
  A1 <- "log(sin(x))"
  expect_identical(sym(A1), Ryacas::ysym("Ln(Sin(x))"))

  ## Test two elements
  A2 <- c("log(x)", "cos(y)")
  expect_identical(sym(A2), Ryacas::ysym(c("Ln(x)", "Cos(y)")))

  ## Test a column vector
  A3 <- matrix(c("tan(a1)", "acos(a2)"), ncol=1)
  expect_identical(sym(A3), 
                   Ryacas::ysym(matrix(c("Tan(a1)", "ArcCos(a2)"), 
                   ncol=1)))

  ## Test a matrix
  A4 <- matrix(paste0("a", 1:6), ncol=2)
  expect_identical(sym(A4), Ryacas::ysym(A4))              
})

test_that("as.matrix() works correctly", {

  A1 <- matrix(c(1, 2, 3, "a", "b", "c"), ncol=2, nrow=3)
  A2  <- Ryacas::ysym(A1)
  expect_identical(A1, as.matrix(A2))

  A3 <- matrix(c(1, 2, 3, "exp(a)", "log(b)", "sqrt(c)"), ncol=2, nrow=3)
  A4 <- matrix(c(1, 2, 3, "Exp(a)", "Ln(b)", "Sqrt(c)"), ncol=2, nrow=3)
  expect_identical(A3, as.matrix(Ryacas::ysym(A4)))
})

test_that("impliedS() works correctly", {

    model1 <- "y ~ c*x + b*m
               m ~ a_a*x
               x ~~ 1*x
               y ~ b_0*1
               m ~ m0*1
               x ~ x0*1"
    
    RAM1 <- metaSEM::lavaan2RAM(model1)

    S1a <- impliedS(RAM1, corr=TRUE)
    S1b <- impliedS(RAM1, corr=FALSE)

    R1a <- list(Sigma = structure(c("1", "a_a*c+b", "c+b*a_a",
                                    "b+a_a*c", "1", "a_a",
                                    "c+b*a_a", "a_a", "1"),
                                  .Dim = c(3L, 3L),
                                  .Dimnames = list(c("y", "m", "x"),
                                                   c("y", "m", "x"))),
                mu = structure(c(0, 0, 0),
                               .Dim = c(1L, 3L),
                               .Dimnames = list("1", c("y", "m", "x"))),
                corr = TRUE)
    R1b <- list(Sigma = structure(c("yWITHy+mWITHm*b^2+(c+b*a_a)^2", "mWITHm*b+a_a*(c+b*a_a)", "c+b*a_a",
                                    "b*mWITHm+(c+b*a_a)*a_a", "mWITHm+a_a^2", "a_a",
                                    "c+b*a_a", "a_a", "1"),
                                  .Dim = c(3L, 3L),
                                  .Dimnames = list(c("y", "m", "x"),
                                                   c("y", "m", "x"))),
                mu = structure(c("b_0+m0*b+x0*(c+b*a_a)",
                                 "m0+x0*a_a","x0"),
                               .Dim = c(1L, 3L),
                               .Dimnames = list("1", c("y", "m", "x"))),
                corr = FALSE)

    expect_identical(S1a, R1a)
    expect_identical(S1b, R1b)
})
