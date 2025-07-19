# Nick's APL_Sim_Tree script translated from Matlab to R and edited for iterations

dyn_mod_fun <- function(iter, gen, init_cond, parameters, init_C = NULL) {
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
  
  alpha <- parameters[["alpha"]][iter, ]  # [alphaAA, alphaAP, alphaAS, alphaPA, alphaPP, alphaPS]
  alphaAA <- alpha$alphaAA; alphaAP <- alpha$alphaAP; alphaAS <- alpha$alphaAS
  alphaPA <- alpha$alphaPA; alphaPP <- alpha$alphaPP; alphaPS <- alpha$alphaPS
  
  beta <- as.numeric(parameters[["beta"]][iter, ])  # [betaA, betaP]
  
  # Initialize vectors
  N_A <- numeric(gen) # number of generations
  L <- numeric(gen)
  N_P <- matrix(0, nrow = 2, ncol = gen)  # Rows: [perennial seeds, perennial adults]
  
  # Set initial conditions
  N_A[1] <- init_cond[1]
  L[1] <- init_cond[2]
  N_P[,1] <- init_cond[3:4] # perennial seeds, perennial adults
  
  # initialize establishment
  E <- matrix(0, nrow = 2, ncol = gen) # Rows: [annual seeds, perennial seeds]
  E[,1] <- e / (1 + beta * L[1]) # first year establishment
  
  # initialize competition
  CA <- numeric(gen)
  CP <- numeric(gen)
  
  # first year competitive effects
  if(is.null(init_C)){
    
    # based on initial conditions
    CA[1] <- 1 + alphaAA * E[1,1] * gA * N_A[1] + 
      alphaAS * gP * E[2,1] * N_P[1,1] +
      alphaAP * N_P[2,1]
    CP[1] <- 1 + alphaPA * E[1,1] * gA * N_A[1] + 
      alphaPS * gP * E[2,1] * N_P[1,1] +
      alphaPP * N_P[2,1]
    
  } else {
    
    # based on provided values
    CA[1] <- init_C[1] # experienced by annual
    CP[1] <- init_C[2] # experienced by perennial
    
  }
  
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
    CA[t] <- 1 + alphaAA * E[1,t] * gA * N_A[t] +
      alphaAS * gP * E[2,t] * N_P[1,t] +
      alphaAP * N_P[2,t]
    CP[t] <- 1 + alphaPA * E[1,t] * gA * N_A[t] +
      alphaPS * gP * E[2,t] * N_P[1,t] +
      alphaPP * N_P[2,t]
  }
  
  # convert to tibble
  out <- tibble(
    .draw = parameters[["draws"]][iter, ]$.draw,
    generation = 1:gen,
    annual_seeds = N_A,
    litter = L,
    perennial_seeds = N_P[1, ],
    perennial_adults = N_P[2, ],
    annual_competition = CA,
    perennial_competition = CP
  )
  
  return(out)
}
