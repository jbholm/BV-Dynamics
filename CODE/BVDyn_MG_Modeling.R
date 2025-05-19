################################################################################
#Title: Mixed effect model of microbiome vs. MET study phase with metagenomic data
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20241022
#See Readme file for details
################################################################################

# Clear workspace and restart R
#rm(list=ls())  
#.rs.restartR()

# Load libraries
library(devtools) 
library(magrittr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(splines)
library(lme4)
library(MuMIn)
library(ggplot2)
library(lattice)
library(broom)
library(performance)
library(car)


# Set seed
set.seed(54321)

# Set working 
workingDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/RDS_Files/"
setwd(workingDir)
outputDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/Figures/"

# Load files
MET.plot.data <- readRDS("BVDyn_MG_Modeling_Figure_2.RDS")
all.species <- readRDS("BVDyn_MG_Modeling_Figure_3.RDS")



### Modeling of DL, L. iners, and BVT with mixed effects regression
poly.model <- function(data, formula) {
  # Fit the full model
  model <- lmer(formula, data = data)
  
  # Calculate AIC, BIC, and Log-Likelihood
  aic <- AIC(model)
  bic <- BIC(model)
  log.lik <- logLik(model)
  
  # Calculate R²
  r.squared <- r.squaredGLMM(model)
  marginal.r2 <- r.squared[1]
  conditional.r2 <- r.squared[2]
  
  # Calculate ICC
  icc <- icc(model)$ICC_adjusted
  
  # Residual diagnostics
  residuals <- resid(model)
  fitted <- fitted(model)
  shapiro.test <- shapiro.test(residuals)
  
  # Bootstrapping
  bootstrap.results <- bootMer(model, function(fit) fixef(fit), nsim = 100)
  
  # Create a summary
  summary <- list(AIC = aic,
                  BIC = bic,
                  LogLik = logLik,
                  Marginal_R2 = marginal.r2,
                  Conditional_R2 = conditional.r2,
                  ICC = icc,
                  Shapiro_Wilk_Test = shapiro.test,
                  Bootstrap_Results = bootstrap.results
  )
  
  return(summary)
}

lacto.eq.MET <- poly.model(formula = Lacto.combined ~ poly(TP, 3) + (1 | SID_2), data = MET.plot.data)
iners.eq.MET <- poly.model(formula = Lactobacillus_iners ~ poly(TP, 3) + (1 | SID_2), data = MET.plot.data)
bv.eq.MET <- poly.model(formula = BV.combined ~ poly(TP, 3) + (1 | SID_2), data = MET.plot.data)

# Modeling of qPCR with mixed effects regression requires commenting out ICC calculation 
#(lines 59 & 75 because of singularity in the data)
qpcr.eq.MET <- poly.model(formula = log2.qPCR ~ poly(TP, 3) + (1 | SID_2), data = MET.plot.data)


### Plot regression of DL, BVT, L. iners, and qPCR
## Use lmer and add random effect to regression line
# Define a function to fit the lmer model and compute predictions with confidence intervals
predict.lmer.ci <- function(data, response) {
  formula <- as.formula(paste(response, "~ poly(TP, 3) + (1 | SID_2)"))
  model.plot <- lmer(formula, data = data)
  
  # Create new data for prediction
  new.data <- data.frame(TP = seq(min(data$TP) - 0.1, max(data$TP) + 0.1, length.out = 100))
  new.data$SID_2 <- data$SID_2[1]  # Dummy grouping variable
  
  # Predict with confidence intervals
  preds <- predict(model.plot, newdata = new.data, re.form = NA, se.fit = TRUE)
  new.data$fit <- preds$fit
  new.data$se <- preds$se.fit
  new.data <- new.data %>% mutate(lwr = fit - 1.96 * se,
                                  upr = fit + 1.96 * se,
                                  response = response)
  
  return(new.data)
}

# Apply the function to each response variable and bind the results
responses <- c("Lacto.combined", "Lactobacillus_iners", "BV.combined", "log2.qPCR")
predictions <- bind_rows(lapply(responses, function(resp) predict.lmer.ci(MET.plot.data, resp)))

## Set plot characteristics
Colors <- c("Lactobacillus_iners" = "orange", "BV.combined" = "cornflowerblue", "Lacto.combined" = "red", "log2.qPCR" = "black")
vertical.lines <- c(0.5, 4.5, 7.5) 

# Plot
theme_set(theme_bw())
MET.poly <- ggplot(MET.plot.data, aes(x = TP)) +
                  geom_jitter(aes(y = Lacto.combined, color = "Lacto.combined"), show.legend = TRUE, width = 0.1) +
                  geom_jitter(aes(y = Lactobacillus_iners, color = "Lactobacillus_iners"), show.legend = TRUE, width = 0.2) +
                  geom_jitter(aes(y = BV.combined, color = "BV.combined"), show.legend = TRUE, width = 0.1) +
                  geom_jitter(aes(y = log2.qPCR, color = "log2.qPCR"), show.legend = TRUE, width = 0.1) +
                  geom_line(data = predictions, aes(y = fit, color = response), size = 1) +
                  geom_ribbon(data = predictions, aes(ymin = lwr, ymax = upr, fill = response), alpha = 0.2, show.legend = FALSE) +
                  geom_vline(xintercept = vertical.lines, linetype = "dotted", color = "black") +
                  scale_fill_manual(values = Colors,
                                    guide = "none") +
                  scale_x_continuous(name = "",
                                     limits = c(-3.5, 9.5),
                                     breaks = seq(-3, 9, 1)) +
                  scale_y_continuous(name = "log2[Abundance + 1]",
                                     limits = c(5, 28.8),
                                     breaks = seq(5, 27.5, 2.5), 
                                     expand = expansion(mult = c(0, 0.05))) +
                  coord_cartesian(clip = "off") +
                  annotate(geom = "text", x = -2 , y = -Inf, vjust = 4, label = "Prior to SBV", size = 3.5) +
                  annotate(geom = "text", x = 0 , y = -Inf, vjust = 4, label = "BVDX", size = 3.5) +
                  annotate(geom = "text", x = 2.5 , y = -Inf, vjust = 4, label = "Early MET", size = 3.5) +
                  annotate(geom = "text", x = 6 , y = -Inf, vjust = 4, label = "Late MET", size = 3.5) +
                  annotate(geom = "text", x = 8.5 , y = -Inf, vjust = 4, label = "After MET", size = 3.5) +
                  theme(aspect.ratio = 5/6,
                        legend.title = element_text(size = 9),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.minor.x = element_blank()) 
MET.poly

ggsave(MET.poly, filename = paste(outputDir, "Figure_2.pdf", sep = ""), 
       height = 8, width = 11)



### Model all species and plot species of interest
# Subtract each observation from BVDX
desired.order <- c(0, -3,-2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9)
all.species.chng <- all.species %>%
                                    group_by(SID_2) %>%
                                    arrange(match(TP, desired.order), .by_group = TRUE) %>%
                                    mutate(across(where(is.numeric), ~ . - first(.))) %>%
                                    slice(-1)
all.species.chng$TP <- as.numeric(all.species.chng$TP)

## Create new function for mixed effect model applied to list of species
poly.model.2 <- function(data, formula) {
  # Fit the model
  model <- lmer(formula, data = data)
  
  # Get the summary of the model
  summary.model <- summary(model)
  
  # Extract coefficients
  coefficients <- fixef(model)
  
  # Calculate marginal and conditional R²
  r.squared <- r.squaredGLMM(model)
  marginal.r2 <- r.squared[1]
  conditional.r2 <- r.squared[2]
  
  # Calculate AIC, BIC, and Log-Likelihood
  aic <- AIC(model)
  bic <- BIC(model)
  logLik <- logLik(model)
  
  # Calculate ICC
  icc <- icc(model)$ICC_conditional
  
  # Shapiro-Wilk test for residuals
  residuals <- resid(model)
  shapiro.test.result <- shapiro.test(residuals)
  
  # Bootstrapping
  has_random_effects <- length(ranef(model)) > 0
  
  bootstrap.results <- tryCatch({
    if (has_random_effects) {
      bootMer(model, function(fit) fixef(fit), nsim = 100)
    } else {
      boot(model, function(fit) coef(fit), nsim = 100)
    }
  }, error = function(e) {
    message("Bootstrapping failed: ", e$message)
    NA
  })
  
  # Extract coefficients
  intercept <- coefficients[1]
  beta1 <- coefficients[2]
  beta2 <- coefficients[3]
  beta3 <- coefficients[4]
  
  # Construct the regression equation
  eq <- paste("y = ", round(beta3, 2), "x³ + ", round(beta2, 2), "x² + ", round(beta1, 2), "x + ", round(intercept, 2))
  
  return(list(eq = eq,
              marginal.r2 = marginal.r2,
              conditional.r2 = conditional.r2,
              aic = aic,
              bic = bic,
              logLik = logLik,
              icc = icc,
              shapiro.test = shapiro.test.result,
              bootstrapping.test = bootstrap.results))
}

## Run function for all species
species.list <- colnames(all.species.chng)[8:174]

# Iterate over each specie and save to dataframe
eq.list.chng <- list()

for (i in species.list) {
  species.model <- poly.model.2(data = all.species.chng,
                                formula = as.formula(paste(i, "~ poly(TP, 3) + (1 | SID_2)")))
  
  # Store the results in the list
  eq.list.chng[[i]] <- list(Model = species.model$eq,
                            Marginal.R2 = species.model$marginal.r2,
                            Conditional.R2 = species.model$conditional.r2,
                            AIC = species.model$aic,
                            BIC = species.model$bic,
                            LogLik = species.model$logLik,
                            ICC = species.model$icc,
                            Shapiro.Wilk.Test = species.model$shapiro.test,
                            Bootstrap.Results = species.model$bootstrapping.test)
}


# Create function to extract specific results from the list and round numeric values
extract.results <- function(results.list, stat) {
  sapply(results.list, function(x) {
    if (is.null(x[[stat]])) return(NA)
    
    result <- x[[stat]]
    
    if (is.list(result)) {
      # Handle specific list cases
      if (stat == "Shapiro.Wilk.Test") {
        result <- result$p.value
      } else if (stat == "Bootstrap.Results") {
        result <- as.numeric(result$t0[1])  
      }}
    
    # Round if result is numeric
    if (is.numeric(result)) {
      return(round(result, 3))
    }
    
    # Return non-numeric result as is
    return(result)})
}


# Create a dataframe from the results list
eq.chng.df <- data.frame(Species = names(eq.list.chng),
                         Model = extract.results(eq.list.chng, "Model"),
                         Marginal_R2 = extract.results(eq.list.chng, "Marginal.R2"),
                         Conditional_R2 = extract.results(eq.list.chng, "Conditional.R2"),
                         AIC = extract.results(eq.list.chng, "AIC"),
                         BIC = extract.results(eq.list.chng, "BIC"),
                         LogLik = extract.results(eq.list.chng, "LogLik"),
                         ICC = extract.results(eq.list.chng, "ICC"),
                         Shapiro_Wilk_Test = extract.results(eq.list.chng, "Shapiro.Wilk.Test.Res"),
                         Bootstrap_Results = extract.results(eq.list.chng, "Bootstrap.Results"),
                         stringsAsFactors = FALSE)

write.csv(eq.chng.df, file = paste(outputDir, "Species_vs_METDay_Poly_Regression_Eqs_NewCounts.csv", sep = ""), fileEncoding = "UTF-8")


### Plot species of interest with (Y = X - BVDX)

species.to.plot <- c("L. crispatus", "L. jensenii", "L. gasseri", "L. mulieris", 
                     "G. vaginalis", "P. amnii", "P. spNov1", "P. spNov2", "P. spNov3", 
                     "P. sp000758925", "B. sp001552935", "D. sp001553355")

all.species.chng.df <- all.species.chng[all.species.chng$Species %in% species.to.plot, ]
all.species.chng.df$Species <- factor(all.species.chng.df$Species, levels = species.to.plot)

# Set plot characteristics
vertical.lines <- c(0.5, 4.5, 7.5) 

a <- xyplot(Count ~ TP | Species, 
            data = all.species.chng.df, as.table = TRUE,
            aspect = 0.9:1,
            groups = SID_2,
            panel = function(x, y, subscripts, groups, ...) {
              panel.xyplot(x, y, ...)
              current.groups <- groups[subscripts]
              fm <- lmer(y ~ poly(x, 3) + (1 | current.groups))
              x.new <- seq(min(x), max(x), length.out = 100)
              
              # Predict with standard errors
              pred <- predict(fm, newdata = data.frame(x = x.new, current.groups = unique(current.groups)[1]), re.form = NA, se.fit = TRUE)
              
              # Calculate confidence intervals
              upper <- pred$fit + 1.96 * pred$se.fit
              lower <- pred$fit - 1.96 * pred$se.fit
              
              # Plot confidence intervals
              panel.lines(x.new, upper, col = "red", lty = 2)
              panel.lines(x.new, lower, col = "red", lty = 2)
              
              # Plot the fitted line
              y.new <- pred$fit
              panel.lines(x.new, y.new, col = "red", lwd = 2)
              
              # Get the panel number
              pn <- panel.number()
              # Convert panel number to row and column indices
              row <- (pn - 1) %/% 3 + 1
              col <- (pn - 1) %% 3 + 1
              # Calculate R² (marginal and conditional R² for fixed effects and random effects)
              r.squared <- r.squaredGLMM(fm)
              marginal.r2 <- r.squared[1]
              conditional.r2 <- r.squared[2]
              
              # Display marginal and conditional R²
              min.y <- min(y) 
              #panel.text(0, min.y, labels = sprintf("Marg. R² = %.2f", marginal.r2), cex = 0.6)
              panel.text(0, min.y + 0.4, labels = sprintf("R² = %.2f", conditional.r2), cex = 0.8)
              panel.abline(v = vertical.lines, col = "black", lty = 2)
            },
            scales = "free",
            jitter.x = TRUE,
            pch = 21,
            cex = .5,
            #min.y <- min(y) - 0.5,
            #ylim = c(min.y, NA),
            xlim = c(-3, 9),
            ylab = list(label = "log2[Abundance + 1]", fontsize = 15),
            xlab = list(label = "Day of MET Treatment", fontsize = 15),
            par.strip.text = list(cex = 0.8, fontface = c("italic", "bold")),
            layout = c(4, 3)
)
a

trellis.device(device = "pdf", file = paste(outputDir, "Figure_3.pdf", sep = ""), height = 9, width = 11.7)
print(a)
dev.off()



