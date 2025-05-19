################################################################################
#Title: Mixed effect model of cytokine vs. MET study phase
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


# Set seed
set.seed(54321)

# Set working 
workingDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/RDS_Files/"
setwd(workingDir)
outputDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/Figures/"

# Load files
immune.4.analysis <- readRDS("BVDyn_Cytokine_Modeling_Figure_4.RDS")


### Make scatter plots with polynomial regression ###
cyto.plot.data <- melt(immune.4.analysis)
colnames(cyto.plot.data) <- c("TP", "UID", "SID_2", "STATUS", "TREATMENT", "NUGENT_CLASS", "Cytokine", "Count")
cyto.plot.data <- cyto.plot.data %>% drop_na(Count)

# Set factors
stat.order <- c("Prior to SBV", "BVDX", "Early MET", "Late MET", "After MET")
cyto.plot.data$STATUS <- factor(cyto.plot.data$STATUS, levels = stat.order)

cyto.order <- c("sEcad", "IL1a", "IL17A", "IL1b", "IL6", "IL8", "IP10", "IFNa2a", "MIG", "MIP1b", "MIP3a", "MMP9")
cyto.plot.data$Cytokine <- factor(cyto.plot.data$Cytokine, levels = cyto.order)

cyto.plot.data$TP <- as.numeric(cyto.plot.data$TP)

## Set plot characteristics
vertical.lines <- c(0.5, 4.5, 7.5) 

## Plot cytokines
a <- xyplot(Count ~ TP | Cytokine, 
            data = cyto.plot.data, as.table = TRUE,
            aspect = 1:1,
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
              #panel.text(0, min.y, labels = sprintf("Marginal R² = %.2f", marginal.r2), cex = 0.6)
              panel.text(0, min.y + 0.2, labels = sprintf("R² = %.2f", conditional.r2), cex = 0.8)
              
              panel.abline(v = vertical.lines, col = "black", lty = 2)
            },
            scales = "free",
            jitter.x = TRUE,
            pch = 21, 
            cex = .5,
            xlim = c(-3, 9),
            ylab = list(label = "log2[Abundance (pg/mL) + 1]", fontsize = 15),
            xlab = list(label = "Day of MET Treatment", fontsize = 15),
            par.strip.text = list(cex = 0.8),
            layout = c(4, 3),
            index.cond = list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
)
a

trellis.device(device = "pdf", file = paste(outputDir, "Figure_4.pdf", sep = ""), height = 8.3, width = 11.7)
print(a)
dev.off()



### Model all cytokines
## Create new function for mixed effect model applied to list of cytokines
poly.model <- function(data, formula) {
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

## Run function for all cyto
immune.4.analysis$TP <- as.numeric(immune.4.analysis$TP)
# Iterate over each cytokines and save to dataframe
eq.list.chng <- list()

for (i in cyto.order) {
  cyto.model <- poly.model(data = immune.4.analysis,
                                formula = as.formula(paste(i, "~ poly(TP, 3) + (1 | SID_2)")))
  
  # Store the results in the list
  eq.list.chng[[i]] <- list(Model = cyto.model$eq,
                            Marginal.R2 = cyto.model$marginal.r2,
                            Conditional.R2 = cyto.model$conditional.r2,
                            AIC = cyto.model$aic,
                            BIC = cyto.model$bic,
                            LogLik = cyto.model$logLik,
                            ICC = cyto.model$icc,
                            Shapiro.Wilk.Test = cyto.model$shapiro.test,
                            Bootstrap.Results = cyto.model$bootstrapping.test)
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
eq.chng.df <- data.frame(Cytokine = names(eq.list.chng),
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

write.csv(eq.chng.df, file = paste(outputDir, "Cytokine_vs_METDay_Poly_Regression_Eqs_NewCounts.csv", sep = ""), fileEncoding = "UTF-8")

