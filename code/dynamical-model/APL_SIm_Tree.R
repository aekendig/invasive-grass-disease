# Nick's script translated from Matlab to R

APL_Sim_Tree <- function(gen, init_cond, parameters) {
  # Function for simulating dynamics with annuals, perennials, and litter, including tree litter effects.
  
  # Unpack parameters
  s <- parameters[[1]]  # [sA, sP, pS, pP]
  sA <- s[1]; sP <- s[2]; pS <- s[3]; pP <- s[4]
  
  y <- parameters[[2]]  # [yA, yP, f]
  yA <- y[1]; yP <- y[2]; f <- y[3]
  
  g <- parameters[[3]]  # [gA, gP]
  gA <- g[1]; gP <- g[2]
  
  e <- parameters[[4]]  # [eA, eP]
  
  decay <- parameters[[5]]  # [bA, bP, d, bT, delta]
  bA <- decay[1]; bP <- decay[2]; d <- decay[3]; bT <- decay[4]; delta <- decay[5]
  
  alpha <- parameters[[6]]  # [alphaA, alphaP, gamma]
  alphaA <- alpha[1]; alphaP <- alpha[2]; gamma <- alpha[3]
  
  beta <- parameters[[7]]  # [betaA, betaP]
  
  # Initialize vectors
  N_A <- numeric(gen) # number of generations
  L <- numeric(gen)
  N_P <- matrix(0, nrow = 2, ncol = gen)  # Rows: [perennial seeds, perennial adults]
  
  # Set initial conditions
  N_A[1] <- init_cond[1]
  L[1] <- init_cond[2]
  N_P[,1] <- init_cond[3:4]
  
  # initialize germination and competition
  E <- matrix(0, nrow = 2, ncol = gen)
  C <- numeric(gen)
  
  E[,1] <- e / (1 + beta * L[1]) # first year establishment
  C[1] <- 1 + alphaA * E[1,1] * gA * N_A[1] + # first year competitive effects
    alphaP * gamma * gP * E[2,1] * N_P[1,1] +
    alphaP * N_P[2,1]
  
  # Time loop
  for (t in 2:gen) {
    # Annual Dynamics
    N_A[t] <- N_A[t-1] * (sA * (1 - gA) + gA * E[1,t-1] * yA / C[t-1])
    
    # Litter Dynamics
    L[t] <- bA * N_A[t-1] * E[1,t-1] * gA +
      bP * (gP * E[2,t-1] * N_P[1,t-1] * (pS * delta + 1 - pS) +
              N_P[2,t-1] * (pP * delta + 1 - pP)) +
      (1 - d) * L[t-1] + bT
    
    # Perennial Dynamics
    M <- matrix(c(sP * (1 - gP) + E[2,t-1] * gP * yP * f / C[t-1],
                  yP / C[t-1],
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
  return(sys)
}
