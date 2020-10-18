library("rstan")
library("bridgesampling")
library("psych")

models_dir <- "stan_models"

data <- data.frame(
  study = c( "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9"),
  t     = c(2.51, 2.39, 2.55, 2.03, 2.23, 2.41, 1.31, 1.92, 2.96),
  N     = c( 100,  150,   97,   99,  100,  150,  200,  100,   50)
)


WFE_null_model  <- rstan::stan_model(file.path(models_dir,"WFE_null.stan"))
FEt_null_loglik <- function(data){
  out <- list()
  out$logml  <- sum(dt(data$t, data$N - 2, log = T))
  class(out) <- "bridge"
  return(out)
}

REt_null_model <- rstan::stan_model(file.path(models_dir,"REt_null_nc.stan"))
WRE_null_model <- rstan::stan_model(file.path(models_dir,"WRE_null_nc.stan"))

FEt_model    <- rstan::stan_model(file.path(models_dir,"FEt.stan"))
WFE_model    <- rstan::stan_model(file.path(models_dir,"WFE.stan")) 

REt_model    <- rstan::stan_model(file.path(models_dir,"REt_nc.stan"))
WRE_model    <- rstan::stan_model(file.path(models_dir,"WRE_nc.stan")) 

#### inference ####

# with prior specified in Bem, Utts & Johnson's 2011 reply:
# centered at zero and 90 quantilie of abs values at .50
#get_prior_sd <- function(x, sd){
#  y <- numeric(1)
#  y <- qnorm(.05, 0, x, lower.tail = F) - sd
#  y
#  }
#round(nleqslv::nleqslv(1, get_prior_sd, sd = .50)$x,3) # = .304

# get data in stan format
dl <- list(
  t      = data$t,
  N      = data$N,
  K      = nrow(data),
  crit_t = qt(.025, data$N-2, lower.tail = F),
  test_type = 2,
  
  # set priors
  prior_mu1     = 0,
  prior_mu2     = .304,
  prior_tau1    = 0,
  prior_tau2    = 1,
  prior_weight1 = 1,
  prior_weight2 = 1,
  
  # type of the weight function
  weight_type = 1,
  
  # some additional settings for the sampler
  tolerance = .001
)

### normal models
m.REt_null<- rstan::sampling(REt_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 1, control = list(adapt_delta = .90))

log_lik.m.REt_null <- bridgesampling::bridge_sampler(m.REt_null)
log_lik.m.FEt_null <- FEt_null_loglik(data)

m.FEt  <- rstan::sampling(FEt_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m.REt  <- rstan::sampling(REt_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

log_lik.m.FEt <- bridgesampling::bridge_sampler(m.FEt)
log_lik.m.REt <- bridgesampling::bridge_sampler(m.REt)

# adding selection models
# type 1 - (t/crit_t)^w
dl$weight_type   <- 1
dl$prior_weight1 <- 1
dl$prior_weight2 <- 5

m1.WFE     <- rstan::sampling(WFE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m1.WRE     <- rstan::sampling(WRE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

m1.WFE_null<- rstan::sampling(WFE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m1.WRE_null<- rstan::sampling(WRE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

log_lik.m1.WFE_null <- bridgesampling::bridge_sampler(m1.WFE_null)
log_lik.m1.WRE_null <- bridgesampling::bridge_sampler(m1.WRE_null)

log_lik.m1.WFE <- bridgesampling::bridge_sampler(m1.WFE)
log_lik.m1.WRE <- bridgesampling::bridge_sampler(m1.WRE)


# type 2 - step functions
dl$weight_type   <- 2
dl$prior_weight1 <- 1
dl$prior_weight2 <- 1

m2.WFE     <- rstan::sampling(WFE_model,      data = dl, iter = 7000, warmup = 5000, chains = 3, cores = 3, control = list(adapt_delta = .99))
m2.WRE     <- rstan::sampling(WRE_model,      data = dl, iter = 7000, warmup = 5000, chains = 3, cores = 3, control = list(adapt_delta = .99))

m2.WFE_null<- rstan::sampling(WFE_null_model, data = dl, iter = 7000, warmup = 5000, chains = 3, cores = 3, control = list(adapt_delta = .99))
m2.WRE_null<- rstan::sampling(WRE_null_model, data = dl, iter = 7000, warmup = 5000, chains = 3, cores = 3, control = list(adapt_delta = .99))

log_lik.m2.WFE_null <- bridgesampling::bridge_sampler(m2.WFE_null)
log_lik.m2.WRE_null <- bridgesampling::bridge_sampler(m2.WRE_null)

log_lik.m2.WFE <- bridgesampling::bridge_sampler(m2.WFE)
log_lik.m2.WRE <- bridgesampling::bridge_sampler(m2.WRE)


# type 3 - beta*pt^w
dl$weight_type   <- 3
dl$prior_weight1 <- 10
dl$prior_weight2 <- 10

m3.WFE     <- rstan::sampling(WFE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m3.WRE     <- rstan::sampling(WRE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

m3.WFE_null<- rstan::sampling(WFE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m3.WRE_null<- rstan::sampling(WRE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

log_lik.m3.WFE_null <- bridgesampling::bridge_sampler(m3.WFE_null)
log_lik.m3.WRE_null <- bridgesampling::bridge_sampler(m3.WRE_null)

log_lik.m3.WFE <- bridgesampling::bridge_sampler(m3.WFE)
log_lik.m3.WRE <- bridgesampling::bridge_sampler(m3.WRE)


# type 4
dl$weight_type   <- 4
dl$prior_weight1 <- 2
dl$prior_weight2 <- 5

m4.WFE     <- rstan::sampling(WFE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m4.WRE     <- rstan::sampling(WRE_model,      data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

m4.WFE_null<- rstan::sampling(WFE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))
m4.WRE_null<- rstan::sampling(WRE_null_model, data = dl, iter = 5000, warmup = 3000, chains = 3, cores = 3, control = list(adapt_delta = .95))

log_lik.m4.WFE_null <- bridgesampling::bridge_sampler(m4.WFE_null)
log_lik.m4.WRE_null <- bridgesampling::bridge_sampler(m4.WRE_null)

log_lik.m4.WFE <- bridgesampling::bridge_sampler(m4.WFE)
log_lik.m4.WRE <- bridgesampling::bridge_sampler(m4.WRE)


### classical BFs
# combined BF reported in Bem's reply 13,669
bridgesampling::bayes_factor(log_lik.m.FEt, log_lik.m.FEt_null) # only fixed-effect models
bridgesampling::bayes_factor(log_lik.m.REt, log_lik.m.REt_null) # only random-effect models
only_unweighted <- bridgesampling::post_prob(log_lik.m.FEt, log_lik.m.REt,
                                             log_lik.m.FEt_null, log_lik.m.REt_null)
sum(only_unweighted[1:2])/sum(only_unweighted[3:4])             # model averaged meta-analysis

### using all models
prior_odds <- rep(c(1, rep(1/4, 4)), 4)
prior_prob <- prior_odds/sum(prior_odds)
marginal_weights <- bridgesampling::post_prob(
  log_lik.m.FEt_null, log_lik.m1.WFE_null, log_lik.m2.WFE_null, log_lik.m3.WFE_null, log_lik.m4.WFE_null,
  log_lik.m.REt_null, log_lik.m1.WRE_null, log_lik.m2.WRE_null, log_lik.m3.WRE_null, log_lik.m4.WRE_null,
  log_lik.m.FEt,      log_lik.m1.WFE,      log_lik.m2.WFE,      log_lik.m3.WFE,      log_lik.m4.WFE,
  log_lik.m.REt,      log_lik.m1.WRE,      log_lik.m2.WRE,      log_lik.m3.WRE,      log_lik.m4.WRE,
  prior_prob = prior_prob
)

marginal_probs <- matrix(
  c(marginal_weights[1],  sum(marginal_weights[2:5]),
    marginal_weights[6],  sum(marginal_weights[7:10]),
    marginal_weights[11], sum(marginal_weights[12:15]),
    marginal_weights[16], sum(marginal_weights[17:20])),
  nrow = 2
)
rownames(marginal_probs) <- c("normal", "weighted")
colnames(marginal_probs) <- c("FE null", "RE null", "FE", "RE")
round(marginal_probs, 3)

# BF for weighted vs unweighted models
sum(marginal_probs[2,])/sum(marginal_probs[1,])

# BF for existence of an effect
sum(marginal_probs[1:2,3:4])/sum(marginal_probs[1:2,1:2])

# assuming weighted fixed-effect model
conditional_WFE <- bridgesampling::post_prob(log_lik.m1.WFE, log_lik.m2.WFE, log_lik.m3.WFE, log_lik.m4.WFE)
round(conditional_WFE,3)

# assuming weighted random-effect model
conditional_WRE <- bridgesampling::post_prob(log_lik.m1.WRE, log_lik.m2.WRE, log_lik.m3.WRE, log_lik.m4.WRE)
round(conditional_WRE,3)

rethinking::precis(m3.WFE, prob = .95)
rethinking::precis(m3.WRE, prob = .95)

samples.m3.WFE <- rstan::extract(m3.WFE)

x <- seq(0,5, .01)
plot(NA, type = "n", xlim = c(0,5), ylim = c(0,1), xlab = "z-scores", ylab = "P(published)", main = "", las = 1)
temp_beta <- sample(samples.m3.WFE$beta, 300)
for(i in 1:length(temp_beta))lines(x, exp(-(temp_beta[i])*(pnorm(abs(x), lower.tail = F)*2)^2), lwd = .1, col = "grey30")
lines(x, exp(-(329.83)*(pnorm(abs(x), lower.tail = F)*2)^2), lwd = 3)

samples.m.FEt <- rstan::extract(m.FEt)
plot(density(samples.m3.WFE$mu), main = "Bem 2011", xlab = "mu", lwd = 2, ylim = c(0, 15), las = 1)
lines(density(samples.m.FEt$mu))
legend("topleft",bty = "n", legend = c("weighted FE", "classical FE"), lwd = c(2,1))


# all-model marginal mean
marginal_weights <- bridgesampling::post_prob(
  log_lik.m.FEt_null, log_lik.m1.WFE_null, log_lik.m2.WFE_null, log_lik.m3.WFE_null, log_lik.m4.WFE_null,
  log_lik.m.REt_null, log_lik.m1.WRE_null, log_lik.m2.WRE_null, log_lik.m3.WRE_null, log_lik.m4.WRE_null,
  log_lik.m.FEt,      log_lik.m1.WFE,      log_lik.m2.WFE,      log_lik.m3.WFE,      log_lik.m4.WFE,
  log_lik.m.REt,      log_lik.m1.WRE,      log_lik.m2.WRE,      log_lik.m3.WRE,      log_lik.m4.WRE,
  prior_prob = prior_prob
)
all_mu <- list(
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  rstan::extract(m.FEt)$mu, rstan::extract(m1.WFE)$mu, rstan::extract(m2.WFE)$mu, rstan::extract(m3.WFE)$mu, rstan::extract(m4.WFE)$mu,
  rstan::extract(m.REt)$mu, rstan::extract(m1.WRE)$mu, rstan::extract(m2.WRE)$mu, rstan::extract(m3.WRE)$mu, rstan::extract(m4.WRE)$mu
)
marginal_mu <- replicate(10000, {
  
  i <- sample(1:length(marginal_weights), 1, prob = marginal_weights)
  
  temp_mu <- all_mu[[i]]
  sample(temp_mu, 1)
})

plot(density(samples.m3.WFE$mu), main = "Bem 2011", xlab = "mu", lwd = 2, ylim = c(0, 15), las = 1)
lines(density(samples.m.FEt$mu))
lines(density(marginal_mu), lwd = 2, col = "blue")
legend("topleft",bty = "n", legend = c("classical FE", "weighted FE", "marginal"), 
       lwd = c(1,2,2), col = c("black","black", "blue"))


