# Nick Clark (https://github.com/nicholasjclark) version
# Basically copy everything that sjwild did for data prep :)
library(tidyverse)
library(tidybayes)
library(brms)
library(cmdstanr)
library(here)
library(ggplot2); theme_set(theme_classic())

# Function to calculate rate per 100k
per_100k <- function(x, y){
  return(x / y * 100000)
}


# Load data
d <- read.csv("Data/wisqars_self_harm_10_to_19.csv",
              na.strings = "--")
d <- d[1:42,]


# Clean data for brms and cmdstanr
d$Estimated.Number <- as.numeric(str_remove_all(d$Estimated.Number, ","))
d$Population <- as.numeric(str_remove_all(d$Population, ","))
d$Standard.Error <- as.numeric(str_remove_all(d$Standard.Error, ","))
d$Lower.95..CI <- as.numeric(str_remove_all(d$Lower.95..CI, ","))
d$Upper.95..CI <- as.numeric(str_remove_all(d$Upper.95..CI, ","))

d$Year <- as.numeric(d$Year)
d$mean <- per_100k(d$Estimated.Number, d$Population)
d$se <- per_100k(d$Standard.Error, d$Population)
d$lb <- per_100k(d$Lower.95..CI, d$Population)
d$ub <- per_100k(d$Upper.95..CI, d$Population)
d$Year <- d$Year - min(d$Year)


# Subset data and model with a Hilbert space squared exponential GP in brms
d_1014 <- d[d$Age.Group == "10 to 14" & !is.na(d$mean),]
mod <- brm(mean | se(se, sigma = TRUE) ~ 
             1 + gp(Year, k = 12, c = 5/4, scale = FALSE),
           data = d_1014,
           chains = 4,
           cores = 4,
           iter = 4000,
           warmup = 2000,
           backend = "cmdstanr")

# View conditional GP trend
conditional_effects(mod)

# Calculate posterior predictions and credible intervals
preds <- as.data.frame(predict(mod))
preds2 <- as.data.frame(predict(mod, probs = c(0.1, 0.9)))
preds3 <- as.data.frame(predict(mod, probs = c(0.25, 0.75)))

# Prep data for plot
preds$Year <- preds2$Year <- preds3$Year <-d_1014$Year
colnames(preds) <- colnames(preds2) <- colnames(preds3) <- 
  c("mean", "se", "lb", "ub", "Year")
preds$group <- "GP95"
preds2$group <- "GP80"
preds3$group <- "GP50"
d_1014$group <- "Raw"

d_plot <- rbind(d_1014[, c("mean", "lb", "ub", "group", "Year")],
                preds[, c("mean", "lb", "ub", "group", "Year")],
                preds2[, c("mean", "lb", "ub", "group", "Year")],
                preds3[, c("mean", "lb", "ub", "group", "Year")])

ggplot(d_plot,
            mapping = aes(x = Year,
                          y = mean,
                          ymin = lb,
                          ymax = ub,
                          lower = lb,
                          upper = ub,
                          middle = mean)) + 
  geom_ribbon(data = d_plot %>%
                dplyr::filter(group == 'GP95'),
              fill = 'darkred',
              alpha = 0.2) +
  geom_ribbon(data = d_plot %>%
                dplyr::filter(group == 'GP80'),
              fill = 'darkred',
              alpha = 0.3) +
  geom_ribbon(data = d_plot %>%
                dplyr::filter(group == 'GP50'),
              fill = 'darkred',
              alpha = 0.4) +
  geom_boxplot(data = d_plot %>%
                 dplyr::filter(group == 'Raw'),
               aes(group = as.factor(Year)),
               stat = 'identity',
               fill = scales::alpha('grey', 0.3),
               width = 0.25) +
  labs(title = "Comparing raw values with a GP",
       x = "Year",
       y = "Rate per 100k") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("Images/brms_gp.png",
       height = 800, width = 1200, units = "px")

# Try a version in mvgam
library(mvgam)
library(marginaleffects)

# This requires the creation of pseudo-series, so it has potential
# to be far less efficient than brms. On the other hand, it should work
# for any response distribution that mvgam supports :\
samp_dat <- do.call(rbind, lapply(1:25, function(x){
  data.frame(yrep = rnorm(NROW(d_1014),
                          mean = d_1014$mean,
                          sd = d_1014$se),
             year = 2000 + d_1014$Year,
             time = rev(1:NROW(d_1014)),
             series = paste0('yrep_', x))
})) %>%
  dplyr::mutate(series = as.factor(series)) %>%
  dplyr::arrange(series, time)

# Plot the samples
ggplot(samp_dat, aes(x = year, y = yrep)) +
  geom_point(alpha = 0.3) +
  labs(y = 'Rate per 100k',
       x = 'Year')

# Set up the trend_map to force all observed 'series' to share the
# same latent process model
trend_map <- data.frame(series = unique(samp_dat$series),
                        trend = 1)

# Fit a model with shared observation variance among the series
mod_mvgam <- mvgam(
  # Observation formula is just the constant
  formula = yrep ~ 1,
  
  # Process formula is where the dynamics occur
  trend_formula = ~
    # No intercept but the GP alpha param should be able
    # to capture this average variability
    gp(year, k = 12, c = 5/4, scale = FALSE),
  
  # Currently mvgam requires autoregressive dynamics
  # for most State Space models; this needs to be relaxed
  # in future versions, but for now we can just use a 
  # suppressed AR1
  trend_model = AR(),
  
  # Suppress the AR1 coefficient and assume the error primarily
  # comes from observation error
  priors = c(prior(normal(0, 0.05),
                   class = ar1,
                   lb = -0.1,
                   ub = 0.1),
             prior(exponential(2),
                   class = sigma)),
  
  # Force all pseudoseries to share the same latent dynamics
  trend_map = trend_map,
  
  # Lognormal observations with a shared observation variance, 
  # which assumes that sampling uncertainty increases proportionally 
  # to the mean (big / dangerous assumption?)
  family = lognormal(),
  share_obs_params = TRUE,
  data = samp_dat,
  
  # MCMC settings
  burnin = 750,
  samples = 1000,
  
  # Don't compute residuals for all 25 series
  residuals = FALSE)

# Model summary shows good convergence overall, though it took 
# much longer than the brms version because of the extra 
# 'observations' and the pesky (un-needed) AR1
summary(mod_mvgam)

# No major issues in parameter estimation
mcmc_plot(mod_mvgam, 
          variable = c('ar1','sigma'), 
          regex = TRUE,
          type = 'hist')

mcmc_plot(mod_mvgam, 
          variable = c('alpha_gp','rho_gp'), 
          regex = TRUE,
          type = 'trace')

# Plot some realisations of the GP
plot(mod_mvgam, type = 'smooths', 
     trend_effects = TRUE,
     realisations = TRUE,
     n_realisations = 25)

# Plot posterior predictions and include the pseudo-observations
plot_predictions(mod_mvgam, condition = 'year',
                 points = 0.4)

# Plot again but show different levels of prediction uncertainty
newdata <- datagrid(model = mod_mvgam,
                    year = seq(2000, 2020, by = 0.5))
preds <- predictions(mod_mvgam,
                     newdata = newdata)
ggplot(preds, aes(year, yrep)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              fill = 'darkred',
              alpha = 0.2) +
  geom_ribbon(data = predictions(mod_mvgam, 
                                 conf_level = 0.8,
                                 newdata = newdata),
              aes(ymin = conf.low,
                  ymax = conf.high),
              fill = 'darkred',
              alpha = 0.3) +
  geom_ribbon(data = predictions(mod_mvgam, 
                                 conf_level = 0.5,
                                 newdata = newdata),
              aes(ymin = conf.low,
                  ymax = conf.high),
              fill = 'darkred',
              alpha = 0.4) +
  geom_point(data = samp_dat, alpha = 0.4,
             size = 0.5) +
  labs(y = 'Rate per 100k',
       x = 'Year',
       title = 'Lognormal State-Space model')

ggsave("Images/mvgam_gp.png",
       height = 800, width = 1200, units = "px")
