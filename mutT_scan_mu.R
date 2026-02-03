# Simulation model of mutation-biased adaptation
# Author: A. Couce
# Description: This program simulates adaptation in a haploid asexual population under a Wright-Fisher
# model for two resistance pathways (TEM-1 and PBP-3), with configurable mutation bias and population size.
# This example explores how the scaling of mutation-driven adaptation with population size depends on the basal
# mutation rate, tested across five values. To test the effect of other parameters, adjust them accordingly.
# Reference: Barber, J; Couce, A. (2025). Mutation-biased adaptation is consequential even in large
# bacterial populations. bioRxiv; doi: https://doi.org/10.1101/2025.06.16.655099 

# Clear previous session
rm(list = ls())

# Start timer
clock <- Sys.time()

# Values of mu to scan
mu_vals <- c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

# Parameters (fixed across scans unless changed below)
replicate <- 100       # number of replicates
m <- 100               # mutation bias coefficient (mutT mutator: 100x for PBP-3)
assy <- 1              # is PBP peak higher?
fmax <- 10             # maximum fitness
nrow <- 9              # 1 + 8 PBP mutations
ncol <- 5              # 1 + 4 TEM mutations

# Fitness landscape
fmat <- matrix(0, nrow = nrow, ncol = ncol)
initial_value <- 1
row_increment <- (fmax * assy - 1) / (nrow - 1)
col_increment <- (fmax - 1) / (ncol - 1)
epis <- 0

# Epistasis (default = 0): implemented by subtracting a fraction of the multiplicative effect between PBP3 and TEM mutations from their additive fitness.

for (i in 1:nrow) {
  for (j in 1:ncol) {
    fmat[i, j] <- initial_value + (i - 1) * row_increment + (j - 1) * col_increment -
      epis*((i - 1) * row_increment * (j - 1) * col_increment)
  }
}

# To control peak height, the additive effects of PBP3 and TEM-1 mutations are truncated at a defined fitness maximum.
fmat[fmat > fmax * assy] <- fmax * assy
fmat[fmat < 0] <- 0

# Visualizing landscape
cols<-c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')
paler_cols <- adjustcolor(cols, alpha.f = 0.7)
image(t(fmat[nrow(fmat):1,]), col=paler_cols)
abline(v=c(0.125, 0.375, 0.625, 0.875))
abline(h=c(0.0625,0.0625*3,0.0625*5,0.0625*7,0.0625*9,0.0625*11,0.0625*13,0.0625*15 ))


# ALGORITHM
# External loop over desired parameter values (here, basal mutation rate)
for (mu in mu_vals) {

  # Plot label for exploratory purposes
  main <- paste0('T100x_mu_', format(mu, scientific = TRUE))

  # Storage matrix for results
  reps <- matrix(0, nrow = replicate, ncol = 17)
  c <- 0

  # Iterate over population sizes
  for (n in c(1e5 * sqrt(sqrt(10))^(0:16))) {
    c <- c + 1
    Nmax <- n
    cat("Simulating N =", Nmax, "with mu =", mu, "\n")

    for (r in 1:replicate) {
      if (r %% 20 == 0) cat("  Replicate", r, "\n")

      pop <- matrix(0, nrow = nrow, ncol = ncol)
      pop[1, 1] <- Nmax
      t <- 1
      # Continue while the most-fit genotype has not reached max. fitness OR not all individuals are in the last row/column 
      # (epistasis may cause populations to get stuck on suboptimal peaks, never fulfilling the first condition)
      while (round(fmat[which(pop == max(pop), arr.ind = TRUE)]) < fmax && !( sum(pop[nrow,],pop[,ncol])==sum(pop) ) ){ 
        t <- t + 1

        # Calculate mean fitness
        w_prom <- sum(pop * fmat) / sum(pop)

        # Reproduction step (weighted sampling)
        for (i in 1:nrow) {
          for (j in 1:ncol) {
            expected_freq <- (pop[i, j] * fmat[i, j]) / (sum(pop) * w_prom)
            pop[i, j] <- rpois(1, Nmax * expected_freq)
          }
        }

        # Mutation step (here, bias applied to PBP pathway)
        for (i in 1:(nrow - 1)) {
          for (j in 1:(ncol - 1)) {
            mb_r <- rpois(1, lambda = pop[i, j] * mu * m)	# PBP pathway
            mb_c <- rpois(1, lambda = pop[i, j] * mu)	    	# TEM pathway

            pop[i, j] <- pop[i, j] - mb_r - mb_c
            pop[i + 1, j] <- pop[i + 1, j] + mb_r
            pop[i, j + 1] <- pop[i, j + 1] + mb_c
          }
        }
      }

      # Quantify adaptation route (TEM% vs total)
      tem <- sum(pop %*% (0:(ncol - 1)))
      pbp <- sum(t(pop) %*% (0:(nrow - 1)))
      reps[r, c] <- tem / (tem + pbp)
    }
  }

  # Visualizing results
  # Smoothing function
  mav <- function(x,n){filter(x,rep(1/n,n), method = c("convolution"), sides=2)}

  x11(width=5.3,height=5.5)
  plot(c(1e5 * sqrt(sqrt(10))^(0:16)), apply(reps,2,function(x) quantile(x,0.5)), log='x', ylim=c(0,1),
  type='o', lwd=2, col='deepskyblue2', xlab='Ne', main=main)
  lines(c(1e5 * sqrt(sqrt(10))^(0:16)), apply(reps,2,function(x) quantile(x,0.75)), ylim=c(0,1),
  type='o', lwd=2, col='deepskyblue2', xlab='Ne')
  lines(c(1e5 * sqrt(sqrt(10))^(0:16)), apply(reps,2,function(x) quantile(x,0.25)), ylim=c(0,1),
  type='o', lwd=2, col='deepskyblue2', xlab='Ne')

  Ne<-c(1e5 * sqrt(sqrt(10))^(0:16))

  # Create and add neutral expectation with 90% confidence intervals
  s_med<-as.matrix(mav(apply(reps,2,median),3))
  s_med[c(1,17)]<-apply(reps[,c(1,17)],2,median)

  s_upp<-as.matrix(mav(apply(reps,2,function(x) quantile(x,0.75)),3))
  s_upp[c(1,17)]<-apply(reps[,c(1,17)],2,function(x) quantile(x,0.75))

  s_low<-as.matrix(mav(apply(reps,2,function(x) quantile(x,0.25)),3))
  s_low[c(1,17)]<-apply(reps[,c(1,17)],2,function(x) quantile(x,0.25))

  # Create and add shaded area depicting interquartile range
  polygon(x = c(Ne,rev(Ne)), y =c(s_upp,rev(s_low)), col= adjustcolor("black", alpha.f = 0.15), border=NA)

  # Print key results for quick monitoring
  cat("\n# KEY OUTPUT (mu =", mu, ")\n")
  print(cbind(Ne, s_med, s_upp, s_low))

  # Write to file for further analyses
  write.table(cbind(mu, Ne, s_med, s_upp, s_low), file = "example.tsv", sep = "\t", row.names = FALSE, col.names = c('mu', 'N', 'med', 'UQ', 'LQ'), quote = FALSE, append = TRUE)
}

# Print elapsed time
cat("Elapsed time:\n")
print(Sys.time() - clock)

