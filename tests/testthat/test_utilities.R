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

    R1a <- list(Sigma = structure(c("1", "a_a*c + b", "a_a*b + c",
                                    "a_a*c + b", "1", "a_a",
                                    "a_a*b + c", "a_a", "1"),
                                  .Dim = c(3L, 3L),
                                  .Dimnames = list(c("y", "m", "x"),
                                                   c("y", "m", "x"))),
                mu = structure(c(0, 0, 0),
                               .Dim = c(1L, 3L),
                               .Dimnames = list("1", c("y", "m", "x"))),
                corr = TRUE)
    R1b <- list(Sigma = structure(c("b^2*mWITHm + yWITHy + (a_a*b + c)^2", 
                                    "a_a*(a_a*b + c) + b*mWITHm", "a_a*b + c", "a_a*(a_a*b + c) + b*mWITHm", 
                                    "a_a^2 + mWITHm", "a_a", "a_a*b + c", "a_a", "1"),
                                  dim = c(3L, 3L),
                                  dimnames = list(c("y", "m", "x"), c("y", "m", "x"))),
                mu = structure(c("b*m0 + b_0 + x0*(a_a*b + c)", "a_a*x0 + m0", "x0"),
                               dim = c(1L, 3L), dimnames = list("1", c("y", "m", "x"))),
                corr = FALSE)

    expect_identical(S1a, R1a)
    expect_identical(S1b, R1b)
})
