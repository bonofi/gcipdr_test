# Required Packages
if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")

library(psych)
library(Matrix)
library(mvtnorm)

set.seed(42)
n_obs <- 500

# 1. GENERATE MOCK EXPERIMENTAL SEED DATA
# Continuous variable (e.g., income baseline)
X_cont <- rnorm(n_obs, mean = 50, sd = 10)
# Latent variable causing continuous/binary linkage
latent_link <- 0.6 * scale(X_cont) + rnorm(n_obs, sd = sqrt(1 - 0.6^2))
# Two correlated binary indicators mapped via thresholds
Y_bin1 <- as.numeric(latent_link > 0.2)
Y_bin2 <- as.numeric((0.5 * latent_link + rnorm(n_obs)) > 0.5)

df_empirical <- data.frame(X_cont = X_cont, Y_bin1 = Y_bin1, Y_bin2 = Y_bin2)

# 2. TRANSFORM CONTINUOUS DATA VIA VAN DER WAERDEN RANK TRANSFORMATION
r_cont <- rank(df_empirical$X_cont)
Z_cont <- qnorm(r_cont / (n_obs + 1))
df_mixed_normal_scores <- data.frame(Z_cont = Z_cont, Y_bin1 = Y_bin1, Y_bin2 = Y_bin2)

# 3. COMPUTE MIXED LATENT CORRELATION MATRIX (Pearson + Biserial + Tetrachoric)
# psych::mixedCor automatically determines pairs based on column classes/unique values
mixed_matrix <- psych::mixedCor(df_mixed_normal_scores, c=1, d=2:3)
Sigma_raw <- mixed_matrix$rho

cat("--- Raw Mixed Correlation Matrix ---\n")
print(Sigma_raw)

# 4. ENFORCE POSITIVE DEFINITENESS (Higham's Alternating Projections Method)
# Hybrid matrices frequently yield negative eigenvalues; nearPD fixes this safely
Sigma_pd <- as.matrix(Matrix::nearPD(Sigma_raw, corr = TRUE)$mat)

cat("\n--- Smooth Near-PD Correlation Matrix ---\n")
print(Sigma_pd)

# 5. EXECUTE NORTA SIMULATION ENGINE
n_sim <- 10000
# Step A: Sample from Multivariate Standard Normal
Z_sim <- mvtnorm::rmvnorm(n_sim, mean = rep(0, 3), sigma = Sigma_pd)

# Step B: Compute Component-Wise Uniform Marginals via Normal CDF
U_sim <- phi(Z_sim) # psych::phi operates identically to pnorm

# Step C: Inverse Mapping to Target Profiles
# Target 1: Parametric approach for continuous variable (Privacy Protected)
fit_mean <- mean(df_empirical$X_cont)
fit_sd   <- sd(df_empirical$X_cont)
X_simulated_cont <- qnorm(U_sim[, 1], mean = fit_mean, sd = fit_sd)

# Target 2 & 3: Latent thresholds derived directly from historical probabilities
p_bin1 <- mean(df_empirical$Y_bin1)
p_bin2 <- mean(df_empirical$Y_bin2)

X_simulated_bin1 <- as.numeric(U_sim[, 2] > (1 - p_bin1))
X_simulated_bin2 <- as.numeric(U_sim[, 3] > (1 - p_bin2))

# Final Synthesized Output
df_synthetic <- data.frame(
  X_cont = X_simulated_cont,
  Y_bin1 = X_simulated_bin1,
  Y_bin2 = X_simulated_bin2
)

cat("\n--- Synthetic Generation Complete ---\n")
print(head(df_synthetic))
