set.seed(1234) # reproducibility

# Simulating the independent variable (arbitrarily
# chosen as ~Normal and the error ~ N(0, sigma^2)
sigma = 1/rgamma(n = 1, shape = 0.5, rate = 1)
x = rnorm(n = 1000, mean = 3, sd = 1.5)
e = rnorm(n = 1000, mean = 0, sd = sqrt(sigma))
V = matrix(sigma * c(10, 0, 0, 10), ncol = 2, nrow = 2)

# Priors of the parameters - intercept and slope
betas = MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = V)

# The regression model
y = betas[1] + (betas[2] * x) + e
da = data.frame(y, x)

# Plotting our data
da %>%
  ggplot(aes(x, y)) +
  geom_point(colour = 'skyblue', size = 2, alpha = 0.75) +
  geom_smooth(method = 'lm', colour = 'grey', linetype = 'dotted') +
  theme_bw() +
  labs(x = 'Independent variable', y = 'Dependent variable')
## `geom_smooth()` using formula 'y ~ x'


# Bayesian model ---------------------------------------------------
# Priors were already assigned as -------------
# Bs ~ Normal(m, V) = Normal(0, sigma^2 * 10)
# Sigma ~ IG(a, b) = IG(0.5, 1)
# Posterior distribution ----------------------
# Betas | Sigma, y ~ N(m*, V*)
# Sigma | y ~ IG(a*, b*), where:
#
# m* = (V^-1 + B'B)^-1 * (V^-1*m + B'Y)
# V* = (V^-1 + B'B)^-1
# a* = a + n/2
# b* = b + (m' V^-1 m + Y'Y - (m*)'(V*)^-1m*)/2
#----------------------------------------------
v_star = solve(solve(V) + t(mm) %*% mm)
m_star = v_star %*% (solve(V) %*% c(0, 0) + t(mm) %*% y)
a_star = 0.5 + n/2
b_star = 1 + (t(c(0, 0)) %*% solve(V) %*% c(0, 0) +
                (t(y) %*% y) - t(m_star) %*% solve(v_star) %*% m_star)/2

# Sampling from the posteriors -----------------------------------------------------
sim = 100000
gamma = rgamma(sim, shape = a_star,
               rate = b_star)
# For the variance
sigma_sim = 1/gamma
# Consider that you have a random variable
# Y ~ Normal(mu, variance), if you multiply it
# by a constant 'c' the variance is multiplied
# by c squared, aka, cY ~ Normal(mu, variance * c^2)
# For the random error, we used that, if Y ~ Normal(0, v*),
# sqrt(sigma) * Y ~ Normal(0, v*sigma), which is our
# target distribution
err = sqrt(sigma_sim)*MASS::mvrnorm(n = sim, mu = c(0, 0), v_star)
# consider now that we are adding the random
# error ~ Normal(0, v*sigma) to a constant
# (the estimated mean for the betas), which will lead us to have
# beta ~ Normal(m_star, v*sigma), as we just added a
# constant to the location of the distribution
params = data.frame(par = c(rep(c(m_star), each = sim) +
                              c(err[,1], err[,2]), sigma_sim))
params$groups = as.factor(rep(c(1:3), each = sim))
params$groups_label = factor(params$groups, labels =
                               c('beta[0]', 'beta[1]', 'sigma^2'))
params_prior = c(betas, sigma)

vline = function(group){
  geom_vline(data = dplyr::filter(params,
                                  groups == group),
             aes(xintercept = params_prior[group]), linetype = 'dotted')
}

params %>%
  ggplot(aes(par)) +
  geom_density(colour = 'skyblue', size = 1.2, alpha = 0.75) +
  #1:3 %>% purrr::map(vline) +
  facet_wrap(~groups_label, scales = 'free',
             labeller = label_parsed) +
  theme_bw() +
  labs(x = 'Parameters', y = 'Density',
       title = 'Posterior densities of the parameters')


seed = 1234
n = 1e4 # number of samples
sigma = 1 # the variance of the target distribution (Normal distribution)
mu = 0 # the mean of the target distribution (Normal distribution)
f_x = function(x){ # the PDF of the target distribution (Normal distribution)
  (1/(sigma * sqrt(2 * pi))) * exp(-(1/2) * ((x - mu)/sigma)^2)
}
range = seq(-sigma * 5, sigma * 5, by = 1/n) # since 99% of values are within 3 sd, we would expect this range to be the range of the target distribution.

# The envelope distribution is a Uniform distribution U~(-3, 3):
Y = runif(length(range), -sigma * 5, sigma * 5) # Uniform distribution that should

# Check if the envelope density requires scaling: 
plot(range, f_x(range), type = 'l')
lines(range, dunif(range, -sigma * 5, sigma * 5))

# The value by which we could scale the Uniform distribution can be approximated by:
scale = max(f_x(range) / dunif(range, -sigma * 5, sigma * 5))
scale
## [1] 3.989423
# Check if the envelope density is scaled correctly: 
plot(range, f_x(range), type = 'l')
lines(range, scale * dunif(range, -sigma * 5, sigma * 5))


# Let's try sampling using both scaled and unscaled values and compare results:
## 1. unscaled proposal distribution:
set.seed(seed)
# Sampling from a standard Uniform distribution for the accept reject rule: 
U = runif(length(range))
accept = c() # to store samples from the normal distribution:
trial = 1 # used in algorithm loop.

# Rejection sampling algorithm:
while (trial <= length(range)) {
  test_u = U[trial] # grab a value from the standard Uniform
  test_x = f_x(Y[trial])
  if (test_u <= test_x)
    accept = rbind(accept, Y[trial])
  trial = trial + 1
}

## 2. scaled envelope distribution:
set.seed(seed)
# Sampling from a standard Uniform distribution for the accept reject rule: 
U = runif(length(range))
accept2 = c() # to store samples from the normal distribution:
trial = 1 # used in algorithm loop.

# Rejection sampling algorithm:
while (trial <= length(range)) {
  test_u = U[trial] # grab a value from the standard Uniform
  test_x = f_x(Y[trial]) / (scale * dunif(Y[trial], -sigma * 5, sigma * 5)) # accounting for scaling
  if (test_u <= test_x)
    accept2 = rbind(accept2, Y[trial])
  trial = trial + 1
}
hist(accept)
