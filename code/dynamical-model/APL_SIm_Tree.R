# Nick's APL_Sim_Tree script translated from Matlab to R and edited for iterations

dyn_mod_fun <- function(iter, gen, init_cond, parameters) {
  # Function for simulating dynamics with annuals, perennials, and litter, including tree litter effects.
  
  s <- parameters[["s"]][iter, ]  # [sA, sP, pS, pP]
  sA <- s$sA; sP <- s$sP; pS <- s$pS; pP <- s$pP
  
  y <- parameters[["y"]][iter, ]  # [yA, yP, f]
  yA <- y$yA; yP <- y$yP; f <- y$f
  
  g <- parameters[["g"]][iter, ]  # [gA, gP]
  gA <- g$gA; gP <- g$gP
  
  e <- as.numeric(parameters[["e"]][iter, ])  # [eA, eP]
  
  decay <- parameters[["decay"]][iter, ]  # [bA, bP, d, bT, delta]
  bA <- decay$bA; bP <- decay$bP; d <- decay$d; bT <- decay$bT; delta <- decay$delta
  
  alpha <- parameters[["alpha"]][iter, ]  # [alphaA, alphaP, gamma]
  alphaA <- alpha$alphaA; alphaP <- alpha$alphaP; gamma <- alpha$gamma
  
  beta <- as.numeric(parameters[[7]][iter, ])  # [betaA, betaP]
  
  # Initialize vectors
  N_A <- numeric(gen) # number of generations
  L <- numeric(gen)
  N_P <- matrix(0, nrow = 2, ncol = gen)  # Rows: [perennial seeds, perennial adults]
  
  # Set initial conditions
  N_A[1] <- init_cond[1]
  L[1] <- init_cond[2]
  N_P[,1] <- init_cond[3:4] # perennial seeds, perennial adults
  
  # initialize germination and competition
  E <- matrix(0, nrow = 2, ncol = gen) # Rows: [annual seeds, perennial seeds]
  CA <- numeric(gen)
  CP <- numeric(gen)
  
  E[,1] <- e / (1 + beta * L[1]) # first year establishment
  CA[1] <- 1 + alphaAA * E[1,1] * gA * N_A[1] + # first year competitive effects
    alphaAS * gP * E[2,1] * N_P[1,1] +
    alphaAP * N_P[2,1]
  CP[1] <- 1 + alphaPA * E[1,1] * gA * N_A[1] + # first year competitive effects
    alphaPS * gP * E[2,1] * N_P[1,1] +
    alphaPP * N_P[2,1]
  
  # Time loop
  for (t in 2:gen) {
    # Annual Dynamics
    N_A[t] <- N_A[t-1] * (sA * (1 - gA) + gA * E[1,t-1] * yA / CA[t-1])
    
    # Litter Dynamics
    L[t] <- bA * N_A[t-1] * E[1,t-1] * gA +
      bP * (gP * E[2,t-1] * N_P[1,t-1] * (pS * delta + 1 - pS) +
              N_P[2,t-1] * (pP * delta + 1 - pP)) +
      (1 - d) * L[t-1] + bT
    
    # Perennial Dynamics
    M <- matrix(c(sP * (1 - gP) + E[2,t-1] * gP * yP * f / CP[t-1],
                  yP / CP[t-1],
                  gP * E[2,t-1] * pS,
                  pP), nrow = 2, byrow = TRUE)
    
    N_P[,t] <- M %*% N_P[,t-1]
    
    # Update E and C
    E[,t] <- e / (1 + beta * L[t])
    C[t] <- 1 + alphaA * E[1,t] * gA * N_A[t] +
      alphaP * gamma * gP * E[2,t] * N_P[1,t] +
      alphaP * N_P[2,t]
  }
  
  # Combine outputs
  sys <- rbind(N_A, L, N_P)  # Rows: 1. Annual seeds, 2. Litter, 3. Perennial seeds, 4. Perennial adults
  
  # convert to tibble
  out <- tibble(
    generation = 1:gen,
    annual_seeds = sys[1, ],
    litter = sys[2, ],
    perennial_seeds = sys[3, ],
    perennial_adults = sys[4, ]
  )
  
  return(out)
}
