library(tidyverse)
library(tidybayes)
library(brms)
library(cmdstanr)
library(here)

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


# subset data and model with spline in brms
# Use defaults for s() and priors because I am lazy
d_1014 <- d[d$Age.Group == "10 to 14" & !is.na(d$mean),]
mod <- brm(mean | se(se, sigma = TRUE) ~ 1 + s(Year),
           data = d_1014,
           chains = 4,
           cores = 4,
           iter = 4000,
           warmup = 3000,
           control = list(adapt_delta = .99),
           backend = "cmdstanr")


# predict mu
fits <- data.frame(fitted(mod))


# prep data for plot
fits$Year <- d_1014$Year
colnames(fits) <- c("mean", "se", "lb", "ub", "Year")
fits$group = "Spline"
d_1014$group = "Raw"

d_plot <- rbind(d_1014[, c("mean", "lb", "ub", "group", "Year")],
                fits[, c("mean", "lb", "ub", "group", "Year")])


# plot
p <- ggplot(d_plot,
            mapping = aes(x = Year,
                          y = mean,
                          ymin = lb,
                          ymax = ub,
                          group = group,
                          fill = group)) + 
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  scale_fill_manual(values = c("grey50", "red"),
                    name = NULL) +
  labs(title = "Comparing raw values with a spline",
       x = "Year",
       y = "Rate per 100k") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")


# save image
ggsave("Images/brms_spline.png", plot = p,
       height = 800, width = 1200, units = "px")



#### Fit local linear trend model ####  



