test_that("lmmprobe runs on simple random intercept data", {
    skip_on_cran() # Skip on CRAN to avoid timeouts if it's slow

    set.seed(123)
    n <- 20
    m <- 4
    p <- 10
    N <- n * m

    # Predictors
    Z <- matrix(rnorm(N * p), ncol = p)
    colnames(Z) <- paste0("z", 1:p)

    # Coefficients
    beta_true <- rep(0, p)
    beta_true[1:2] <- c(1, -1)

    # Random effects
    sigma_b <- 0.5
    b <- rnorm(n, 0, sigma_b)
    ID <- rep(1:n, each = m)
    b_vec <- b[ID]

    # Error
    sigma_e <- 1
    Y <- as.vector(Z %*% beta_true + b_vec + rnorm(N, 0, sigma_e))

    # Inputs for lmmprobe
    # V should be matrix of IDs for random intercept case (per docs/examples)
    V <- matrix(ID, ncol = 1)
    ID_data <- as.numeric(ID)

    # Run lmmprobe
    # We use minimal cpu to avoid heavy parallel overhead in test
    # Note: lmmprobe might require specific backend, we'll try cpus=1
    # If cpus=1 fails (some snowfall implementations assume >=2 or specific cluster), we might need 2.
    # But usually 1 is fine or 2. We'll try 1.

    fit <- lmmprobe(
        Y = matrix(Y, ncol = 1),
        Z = Z,
        V = V,
        ID_data = ID_data,
        # Pass test data same as train for simplicity in this basic test
        Y_test = matrix(Y, ncol = 1),
        Z_test = Z,
        V_test = V,
        ID_test = ID_data,
        alpha = 0.05,
        maxit = 5, # Very few iterations for speed
        ep = 1.0, # Loose convergence to stop early
        cpus = 1,
        B = 3
    )

    expect_type(fit, "list")
    expect_true("beta" %in% names(fit))
    expect_true("gamma" %in% names(fit))
    expect_equal(length(fit$beta), p)

    # Check if it picked up something (loose check)
    # beta[1] should likely be positive, beta[2] negative
    # But with 5 iterations and defaults it might not be converged.
    # Just checking structure for now.
})
