################################################################################
#Title: Finalizing BV Dynamincs Study Figure - MET
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20241015
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
library(ggplot2)
library(tidytext)
library(ggnewscale)

# Set seed
set.seed(54321)

# Set working 
setwd("~/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/Code")
workingDir <- "../RDS_Files/"
setwd(workingDir)
outputDir <- "../Figures/"

# Load files
all.days <- readRDS("BVDyn_Study_Figure_1.RDS")

### Make scatter plots ###
study.phase <- c("Prior to SBV" = "#66C2A5", "BVDX" = "#D73027", "Early MET" = "#8DA0CB", 
                 "Late MET" = "#08519C", "After MET" = "#FC8D62")
status.shapes <- c(15, 16, 7, 9, 17)

#study.phase <- c("Prior to SBV" = "black", "BVDX" = "#D73027", "Early MET" = "black", 
#                 "Late MET" = "black", "After MET" = "black")
#status.shapes <- c(22, 21, 22, 22, 22)

cust.theme <- theme_bw() + theme(text = element_text(size = 12),
                                 aspect.ratio = 3/8,
                                 legend.title = element_text(hjust = 0.5))

# Adjust study day so points for consecutive days don't overlap
all.days$SERIAL_adj <- with(all.days, SERIAL + as.numeric(factor(STATUS)) * 0.2 - 0.5)
all.days <- all.days %>%
                        relocate(SERIAL_adj, .after = SERIAL) 

Fig.1 <- ggplot(all.days, aes(y = SID, x = SERIAL_adj, fill = STATUS, shape = STATUS)) +
                geom_point(size = 5, color = "black", stroke = 1, alpha = 0.8) + 
                xlab("Study Day") +
                xlim(c(0, 70)) +
                scale_fill_manual(values = study.phase) +
                scale_shape_manual(values = status.shapes) +
                guides(fill = guide_legend(title = "Study phase")) +
                cust.theme +
                theme(legend.margin = margin(0, 0, 0, -6))

Fig.1
ggsave(Fig.1, filename = paste(outputDir, "Figure_1.pdf", sep = ""), 
       height = 8, width = 11)