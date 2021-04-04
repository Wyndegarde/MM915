########### Required packages and Data ##############
library(tidyverse)
library(ggplot2)
library(patchwork) 
library(gstat)
library(latticeExtra)
library(sp)

require(spdep)
require(sf)
require(rgdal)
require(rgeos)
library(shapefiles)
library(R2WinBUGS)
library(rjags)
library(GGally)

load("london.RData")
########### Question 1 ################


# Create a subset containing 2005

london@data <- london@data[london@data$year == 2005, ]
plot(london)
# Create new column: SIR

london$SIR <- london$Observed / london$Expected

# Because the Poisson regression uses a log link, we want to look at the Log SIR

london$logSIR <- log(london$SIR)


# Create SpPlot to explore SIR
london_simple <- st_as_sf(london)


spplot(london, 'SIR', col.regions = rev(heat.colors(200)))

plot_box <- tibble(
  xmin = st_bbox(london)[1],
  xmax = st_bbox(london)[3],
  ymin = st_bbox(london)[2],
  ymax = st_bbox(london)[4],
  xrange = xmax - xmin,
  yrange = ymax - ymin
)

ggplot(london_simple) +
  geom_sf(aes(fill = logSIR, colour = logSIR)) +
  scale_fill_viridis_c(option = "C") +
  scale_colour_viridis_c(option = "C") +
  labs(title = "Spatial plot of London", subtitle  = "Using levels of Log(SIR) as measurement") +
  annotate(
    "rect",
    ymin = plot_box$ymin + plot_box$yrange * 0.54,
    ymax = plot_box$ymin + plot_box$yrange * 0.75,
    xmin = plot_box$xmin + plot_box$xrange * 0.5,
    xmax = plot_box$xmin + plot_box$xrange * 0.65,
    colour = "black",
    alpha = 0.05,
    size = 0.75
  ) +
  annotate(
    "curve",
    xend = plot_box$xmin + plot_box$xrange * 0.65,
    yend = plot_box$ymin + plot_box$yrange * 0.75 ,
    x = plot_box$xmin + plot_box$xrange * 0.7,
    y = plot_box$ymin + plot_box$yrange * 0.98,
    arrow = arrow(length = unit(0.03, "npc")),
    curvature = 0.3,
    size = 0.75
  ) +
  annotate(
    "rect",
    ymin = plot_box$ymin + plot_box$yrange * 0.48,
    ymax = plot_box$ymin + plot_box$yrange * 0.63,
    xmin = plot_box$xmin + plot_box$xrange * 0.34,
    xmax = plot_box$xmin + plot_box$xrange * 0.455,
    colour = "black",
    alpha = 0.05,
    size = 0.75
  ) +
  annotate(
    "curve",
    xend = plot_box$xmin + plot_box$xrange * 0.455,
    yend = plot_box$ymin + plot_box$yrange * 0.63 ,
    x = plot_box$xmin + plot_box$xrange * 0.7,
    y = plot_box$ymin + plot_box$yrange * 0.98,
    arrow = arrow(length = unit(0.03, "npc")),
    curvature = 0.3,
    size = 0.75
  ) +
  annotate(
    "text",
    x = plot_box$xmin + plot_box$xrange * 0.86,
    y = plot_box$ymin + plot_box$yrange * 0.97,
    label = "These two areas show signs \n of possible spatial dependence."
  ) +
  annotate(
    "text",
    x = plot_box$xmin + plot_box$xrange * 0.635,
    y = plot_box$ymin + plot_box$yrange * 0.73,
    label = "1",
    size = 5
  ) +
  annotate(
    "text",
    x = plot_box$xmin + plot_box$xrange * 0.44,
    y = plot_box$ymin + plot_box$yrange * 0.61,
    label = "2",
    size = 5
  ) +
  theme_void()

par(mar = c(5, 4, 4, 2), mfrow  = c(1, 1))
hist(london$SIR,
     main = "Respiratory Disease SIR Across London",
     xlab = "Standardised incidence ratio")
H

spplot(london, "logSIR")



ggplot(london_simple, aes(JSA, logSIR)) + geom_point(colour = 'black') +
  labs(title = "Log Standard Incidence Ratio \nplotted against three key variables") +
  xlab("Job Seekers Allowance") -> p8

ggplot(london_simple, aes(PM25, logSIR)) + geom_point(colour = 'black') +
  xlab("Average fine particulate matter concentration") -> p9

ggplot(london_simple, aes(Price, logSIR)) + geom_point(colour = 'black') +
  xlab("Annual average sale price of homes") -> p10

layout = '
AABB
#CC#
'
p8 + p9 + p10 + plot_layout(design = layout)

ggpairs(london@data[, c("PM25", 'JSA', 'Price', 'logSIR')])
# Looking at these ggpairs plots, there does not appear to be any worries about multi-collinearity
# But also no strong relationships at all really.


########### Question 2 ##############################

# Use to avoid running BUGS code
load("LondonRegressionModel.RData")

Data <-
  list(
    Y = london$Observed,
    E = london$Expected,
    N = nrow(london@data),
    JSA = london@data$JSA - mean(london@data$JSA),
    PM25 = london@data$PM25 - mean(london@data$PM25),
    Price = london@data$Price - mean(london@data$Price)
  )

Init = list(list(
  beta0 = 0,
  beta1 = 0,
  beta2 = 0,
  beta3 = 0
),
list(
  beta0 = 0,
  beta1 = 0,
  beta2 = 0,
  beta3 = 0
))
# We set initial value because we are taking logs within our calculation
# Typically set initial value somewhere around the mean of your prior distribution

london.sim <- bugs(
  Data,
  parameters.to.save = c("beta0", "beta1", "beta2", "beta3"),
  model.file = "poisson-regression.txt",
  inits = Init,
  n.iter = 7500,
  n.burnin = 5000,
  n.thin = 1,
  n.chains = 2,
  bugs.directory = "C:/Users/liamh/OneDrive/Documents/WinBUGS14"
)

print(london.sim, digits.summary = 3)

london.samples <- as.data.frame(london.sim$sims.list)

ggplot(london.samples) + geom_line(aes(1:nrow(london.samples), beta0)) +
  labs(title = "Traceplots of saved parameters",
       subtitle = "Convergence checks",
       x = "Index") -> TP0

ggplot(london.samples) + geom_line(aes(1:nrow(london.samples), beta1)) +
  labs(x = "Index") -> TP1

ggplot(london.samples) + geom_line(aes(1:nrow(london.samples), beta2)) +
  labs(x = "Index")  -> TP2

ggplot(london.samples) + geom_line(aes(1:nrow(london.samples), beta3)) +
  labs(x = "Index")  -> TP3

(TP0 + TP1) / (TP2 + TP3)


geweke.plot(london.sim)

geweke.diag(london.samples)

# Pearson Residuals

beta_samples <-
  sapply(london.sim$sims.list[grepl("beta", names(london.sim$sims.list))], mean)

london.sim$summary[1, 1]

# For the purposes of this I'll assume 2 covariates
# you can extend this to any number of covariates
lp <-
  (beta_samples[1]) + (beta_samples[2] * london@data$JSA) + 
  (beta_samples[3] *london@data$PM25) + (beta_samples[4] * london@data$Price)

mu <- london@data$Expected * exp(lp)

london_res <- (london@data$Expected - mu) / sqrt(mu)

var(london_res) # Suggests overdispersion - making CARS model appropriate (state why in paper)

ggplot() + geom_point(aes(mu, london_res)) + xlab("Mu") + ylab("Residuals") +
  labs(title = "Fitted values vs residuals plot", subtitle = "Checking for Over/under dispersion")

# next do the Moran's I

moran.test(london_res, nb2listw(poly2nb(london)), alternative = "two.sided")

# Significant result so There is spatial dependence within the dataset.

save(london.sim, file = "LondonRegressionModel.RData")


x <- 1
########### Question 3 - CARS model ###############

# We use a poisson CARS model if we have tried to fit a normal poisson regression model
# And the assumptions have not been met.

# Get our adjacency matrix
W <- nb2mat(poly2nb(london), style = "B") # This is Queen adjacency

# Extract the index of neighbouring regions
inds <- lapply(1:nrow(W), function(i)
  which(W[i,] == 1))

# convert previous result into a vector
Adj <- Reduce("c", inds)

# calculate the number of neighbours for each region
Num.Adj <- rowSums(W)

# Calculate the total number of neighbouring pairs (should equal length of Adj)
SumNumNeigh <- sum(Num.Adj)


#Use to load CARS model
load("LondonCARSModel.RData")

# set up data - this includes my spatial information

Data.CARS <-
  list(
    Y = london$Observed,
    E = london$Expected,
    N = nrow(london@data),
    JSA = london@data$JSA - mean(london@data$JSA),
    PM25 = london@data$PM25 - mean(london@data$PM25),
    Price = london@data$Price - mean(london@data$Price),
    Adj = Adj,
    Num = Num.Adj,
    SumNumNeigh = SumNumNeigh
  )

Init.CARS <-
  list(
    list(
      beta0 = 0,
      beta1 = 0,
      beta2 = 0,
      beta3 = 0,
      phi = rep(1, nrow(london))
    ),
    list(
      beta0 = 0,
      beta1 = 0,
      beta2 = 0,
      beta3 = 0,
      phi = rep(1, nrow(london))
    )
  )
# In BUGS Code, Tau.v Was 0.1,0.1 previously.
# Tau.v tells you about the level of spatial dependence in the data.



london.CARS <- bugs(
  Data.CARS,
  parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "phi"),
  model.file = "CARSmodel.txt",
  inits = Init.CARS,
  n.iter = 23000,
  n.burnin = 5000,
  n.thin = 3,
  n.chains = 2,
  bugs.directory = "C:/Users/liamh/OneDrive/Documents/WinBUGS14"
)
# Set this so I am aiming to get 2000 effective samples, with 6000 total samples.
# 5000 burn in.
# So, k = 6000/2000 = 3. iterations = K * total samples +  burn-in = 3 * 6000 + 5000 = 23000

print(london.CARS, digits.summary = 3)

min(london.CARS$summary[, 8])
max(london.CARS$summary[, 8])

london.CARS.samples <- as.data.frame(london.CARS$sims.list)

ggplot(london.CARS.samples) + geom_line(aes(1:nrow(london.CARS.samples), beta0)) +
  labs(title = "Traceplots of saved parameters",
       subtitle = "Convergence checks (CARS model)",
       x = "Index") -> CARSTP0

ggplot(london.CARS.samples) + geom_line(aes(1:nrow(london.CARS.samples), beta1)) +
  labs(x = "Index") -> CARSTP1

ggplot(london.CARS.samples) + geom_line(aes(1:nrow(london.CARS.samples), beta2)) +
  labs(x = "Index")  -> CARSTP2

ggplot(london.CARS.samples) + geom_line(aes(1:nrow(london.CARS.samples), beta3)) +
  labs(x = "Index")  -> CARSTP3

(CARSTP0 + CARSTP1) / (CARSTP2 + CARSTP3)


geweke.plot(london.CARS)

geweke.diag(london.CARS.samples)


# Pearson Residuals

beta_samplesCARS <-
  sapply(london.CARS$sims.list[grepl("beta", names(london.CARS$sims.list))], mean)

london.CARS$summary[1, 1]

lpCARS <-
  beta_samplesCARS[1] + (beta_samplesCARS[2] * london@data$JSA) + 
  (beta_samplesCARS[3] *london@data$PM25) + (beta_samplesCARS[4] * london@data$Price)


muCARS <- london@data$Expected * exp(lpCARS)

london_resCARS <- (london@data$Observed - muCARS) / sqrt(muCARS)

var(london_resCARS)
length(london_resCARS)
ggplot() + geom_point(aes(muCARS, london_resCARS)) + xlab("Mu") + ylab("Residuals") +
  labs(title = "Plot of fitted values against residuals", subtitle  = "Testing mean = variance assumption")

# next do the Moran's I

moran.test(london_resCARS, nb2listw(poly2nb(london)), alternative = "two.sided")

save(london.CARS, file = "LondonCARSModel.RData")

########### Question 4 #######################

# Poison Reg
pois_reg_DIC <- london.sim$DIC
pois_reg_var <- var(london_res) # Suggests overdispersion - making CARS model appropriate (state why in paper)
pois_reg_Moran <- moran.test(london_res, nb2listw(poly2nb(london)), alternative="two.sided")
pois_reg_moran_stat <- pois_reg_Moran$statistic

pois_reg_PM25_mean <- exp(london.sim$summary[3,1])
pois_reg_PM25_lower <- exp(london.sim$summary[3,3])
pois_reg_PM25_upper <- exp(london.sim$summary[3,7])

# Poison_CAR

pois_CAR_DIC <- london.CARS$DIC
pois_CAR_var <- var(london_resCARS) 
pois_CAR_Moran <- moran.test(london_resCARS, nb2listw(poly2nb(london)), alternative="two.sided")
pois_CAR_moran_stat <- pois_CAR_Moran$statistic

pois_CAR_PM25_mean <- exp(london.CARS$summary[3,1])
pois_CAR_PM25_lower <- exp(london.CARS$summary[3,3])
pois_CAR_PM25_upper <- exp(london.CARS$summary[3,7])

model_comparison <- data.frame(
  Poisson_regression = 
    c(pois_reg_DIC, pois_reg_var, pois_reg_moran_stat, pois_reg_PM25_lower, pois_reg_PM25_mean, pois_reg_PM25_upper),
  Poisson_CARS_regression =
    c(pois_CAR_DIC,pois_CAR_var,pois_CAR_moran_stat, pois_CAR_PM25_lower, pois_CAR_PM25_mean, pois_CAR_PM25_upper),
  row.names = c("DIC","Variance","Moran Statistic","PM25 lower CI","PM25 mean","PM25 upper CI"))

colnames(model_comparison) <- c("Regression model","CARS model")

model_comparison


pois_CAR_mean <- exp(london.CARS$summary[1:4,1])
pois_CAR_lower <- exp(london.CARS$summary[1:4,3])
pois_CAR_upper <- exp(london.CARS$summary[1:4,7])
########### Extra Code ##############


par(mar = c(0, 0, 0, 0), mfrow = c(1, 2))
nb_queen <- poly2nb(london)
nb_rook <- poly2nb(london, queen = F)
plot(london,
     border = "grey60",
     axes = F,
     bty = NULL)
plot(
  nb_queen,
  coordinates(london),
  pch = 19,
  cex = 0.6,
  add = TRUE
)
plot(london,
     border = "grey60",
     axes = F,
     bty = NULL)
plot(
  nb_rook,
  coordinates(london),
  pch = 19,
  cex = 0.6,
  add = TRUE
)
# Don't run again
london.CARS <- bugs(
  Data.CARS,
  parameters.to.save = c("beta0", "beta1", "beta2", "beta3","phi"),
  model.file = "CARSmodel.txt",
  inits = Init.CARS,
  n.iter = 7500,
  n.burnin = 5000,
  n.thin = 1,
  n.chains = 2,
  bugs.directory = "C:/Users/liamh/OneDrive/Documents/WinBUGS14") 
# This model resulted in very few effective samples. So introduced thinning


write.csv(model_comparison, "modelComparison.csv")



########### BUGS Files ################

# Poisson Regression

model {
  for(i in 1:N)
  {
    Y[i]    ~ dpois(mu[i])
    log(mu[i])   <- log(E[i]) + beta0 + JSA[i]*beta1 + PM25[i]*beta2 + Price[i]*beta3
  }
  beta0   ~ dnorm(0, 0.01)
  beta1   ~ dnorm(0, 0.01)
  beta2   ~ dnorm(0, 0.01)
  beta3   ~ dnorm(0, 0.01)
}

# CARS model

model {
  for(i in 1:N)
  {
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(E[i]) + beta0 + JSA[i]*beta1 + PM25[i]*beta2 + Price[i]*beta3 + phi[i]
  }
  
  phi[1:N] ~ car.normal(Adj[], weights[], Num[], tau.v)
  for(k in 1:SumNumNeigh){
    weights[k] <- 1
  }
  
  # priors
  tau.v   ~ dgamma(0.5, 0.0005)
  beta0   ~ dnorm(0, 0.001)
  beta1   ~ dnorm(0, 0.001)
  beta2   ~ dnorm(0, 0.001)
  beta3   ~ dnorm(0, 0.001)
}