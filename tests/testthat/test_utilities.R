## Skip the tests in cran as it is difficult to install the required environments, such as python and sympy.
skip_on_cran()

test_that("impliedS() works correctly", {

    model1 <- "y ~ c*x + b*m
               m ~ a_a*x
               x ~~ 1*x
               y ~ b_0*1
               m ~ m0*1
               x ~ x0*1
               ## Adhoc parameter names to test special characters (E, I) in sympy
               m ~~ EvarmE*m
               y ~~ IvaryI*y"
    
    RAM1 <- metaSEM::lavaan2RAM(model1)

    S1a <- impliedS(RAM1, corr=TRUE)
    S1b <- impliedS(RAM1, corr=FALSE)

    R1a <- list(Sigma = structure(c("1", "a_a*c + b", "a_a*b + c",
                                    "a_a*c + b", "1", "a_a",
                                    "a_a*b + c", "a_a", "1"),
                                  .Dim = c(3L, 3L),
                                  .Dimnames = list(c("y", "m", "x"),
                                                   c("y", "m", "x"))),
                Mu = structure(c(0, 0, 0),
                               .Dim = c(1L, 3L),
                               .Dimnames = list("1", c("y", "m", "x"))),
                corr = TRUE)
    R1b <- list(Sigma = structure(c("IvaryI + EvarmE*b^2 + (a_a*b + c)^2", 
                                    "a_a*(a_a*b + c) + EvarmE*b", "a_a*b + c",
                                    "a_a*(a_a*b + c) + EvarmE*b", 
                                    "a_a^2 + EvarmE", "a_a", "a_a*b + c", "a_a", "1"),
                                  dim = c(3L, 3L),
                                  dimnames = list(c("y", "m", "x"), c("y", "m", "x"))),
                Mu = structure(c("b*m0 + b_0 + x0*(a_a*b + c)", "a_a*x0 + m0", "x0"),
                               dim = c(1L, 3L), dimnames = list("1", c("y", "m", "x"))),
                corr = FALSE)

    expect_identical(S1a$Sigma, R1a$Sigma)
    expect_identical(S1b$Sigma, R1b$Sigma)
    expect_identical(S1b$Mu, R1b$Mu)
})
