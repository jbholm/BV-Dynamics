################################################################################
#Title: Spearman correlation of immune markers and bacterial abundances
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20240808
#See Readme file for details
################################################################################

# Clear workspace and restart R
#rm(list=ls())  
#.rs.restartR()

# Load libraries
library(devtools) 
library(magrittr)
require(readxl)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# Set seed
set.seed(54321)

# Set working 
workingDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/RDS_Files/"
setwd(workingDir)
outputDir <- "/Users/amandawilliams/IGS Dropbox/Amanda Williams/BV_Dynamics_2022/Manuscript_2024/Figures/"


# Load files
all.days <- read.table("BVDyn_Spearman_Corr_Input.RDS", header = T, row.names = 1, check.names = FALSE, sep = ',')


######################## For all participants by Status ########################

all.days <- all.days %>% rename(Ca_Lachnocurva_vaginae = UBA629_sp005465875)

species.list <- names(all.days)[19:739]
species.list <- data.frame(species.list)


# Create filter variable (Status)
stat <- unique(all.days$STATUS)

for (i in stat) {
  data.2 <- all.days %>% 
                      filter(STATUS == i) 

    ## Pull out immune data and format for analysis
    immune.4.analysis <- data.2[ , c(3, 749:760)]
    immune.4.analysis <- data.frame(immune.4.analysis)
    rownames(immune.4.analysis) <- immune.4.analysis$UID
    immune.4.analysis[is.na(immune.4.analysis)] <- 0
    immune.4.analysis <- immune.4.analysis[ , -1]
    
    # Create list of filtered variables (cytokines)
    cytokines <- colnames(immune.4.analysis)
    cytokines <- data.frame(cytokines) %>%
                                          rename(Cytokines = 1)
    
    
    ## Pull out amplicon data and format for analysis
    amp.4.analysis <- data.2[ , c(3, 19:748)]
    amp.4.analysis <- data.frame(amp.4.analysis)
    rownames(amp.4.analysis) <- amp.4.analysis$UID
    amp.4.analysis <- amp.4.analysis[ , -1]
    
    # Filter amplicon data for species with >X reads
    amp.4.analysis <- amp.4.analysis %>% drop_na(Lactobacillus_iners)
    amp.4.analysis[is.na(amp.4.analysis)] <- 0
    amp.4.analysis <- amp.4.analysis[as.logical(rowSums(amp.4.analysis != 0)), ]
    
    # Create list of filtered variables (species)
    species <- colnames(amp.4.analysis)
    species <- data.frame(species) %>%
                                      rename(Species = 1)
    
    ## Join data frames for analysis
    joined.data <- amp.4.analysis %>% merge(immune.4.analysis, by = 0)
    rownames(joined.data) <- joined.data$Row.names
    joined.data$Row.names <- NULL
    
    # Perform log10 transformation
    joined.data <- log2(joined.data + 1)
    
    
    ## Create data frame of every combination of metabolite and class with leading data frame name
    combinations <- crossing(species, cytokines)
    combinations.merge <- combinations %>% rename(Species = 1,
                                                  Cytokines = 2)
    combinations <- t(combinations)
    combinations <- data.frame(combinations)
    names(combinations) = str_sub(names(combinations), 2)
    
    ## Run the Spearman correlation test
    spearman.output <- lapply(combinations, function(x) {
                              (cor.test(joined.data[ ,x[1]], joined.data[ ,x[2]], 
                                        method = "spearman",
                                        exact = FALSE))
                            })
    
    ## Format spearman output list into data frame
    spearman.output.df <- data.frame(matrix(unlist(spearman.output), nrow = length(spearman.output), byrow = TRUE))
    
    spearman.output.df <- spearman.output.df[, c(1:5)]
    spearman.output.df <- spearman.output.df %>%
                                                rename(S.statistic = 1,
                                                       pavalue = 2,
                                                       rho = 3,
                                                       null.value = 4,
                                                       type = 5)
    
    combinations.merge <- data.frame(rowname = rownames(combinations.merge), combinations.merge, row.names = NULL)
    spearman.output.df <- data.frame(rowname = rownames(spearman.output.df), spearman.output.df, row.names = NULL)
    
    spearman.output.df <- spearman.output.df %>% full_join(combinations.merge, by = "rowname")
    spearman.output.df <- spearman.output.df[ , -1]
    spearman.output.df$rho <- as.numeric(spearman.output.df$rho)
    
   ## Unmelt ambient table and format for plotting
    corr.results <- spearman.output.df[ , c(6, 7, 3)]
    corr.results <- corr.results %>% drop_na(rho)
    corr.results <- corr.results[order(corr.results$Species, corr.results$Cytokines), ]
    
    corr.results <- dcast(corr.results, Species ~ Cytokines)
    corr.results <- corr.results %>% 
                                    arrange(Species)
    rownames(corr.results) <- corr.results$Species
    corr.results <- corr.results[ , -1]
    
    # Transform df into matrix
    corr.results.mat <- as.matrix(corr.results)
    mode(corr.results.mat)
    
    
  
    ## Subset species of interest (BV - related) and set order
    
    # Subset species of interest for analysis from data frame
    t.corr.results <- t(corr.results) %>%
                                          data.frame()
    lacto.gard <- t.corr.results[ , (names(t.corr.results) %in% order)] 
    lacto.gard.mat <- t(lacto.gard)
    
    # Correct taxa names
    rownames(lacto.gard.mat) <- gsub("_", " ", rownames(lacto.gard.mat))
    rownames(lacto.gard.mat) <- gsub("(^[A-Za-z])[a-z]+\\s([a-z]+)", "\\1. \\2", rownames(lacto.gard.mat))
  
}


## 08Nov2024 -- Johanna added

all.days <- readRDS("BVDyn_Spearman_Corr_Input.RDS")
all.days <- all.days %>% dplyr::rename(Ca_Lachnocurva_vaginae = UBA629_sp005465875)

species.list <- names(all.days)[19:739]
species.list <- data.frame(species.list)


data.2 <- all.days

## Pull out immune data and format for analysis
immune.4.analysis <- data.2[ , c(3, 749:760)]
immune.4.analysis <- data.frame(immune.4.analysis)
rownames(immune.4.analysis) <- immune.4.analysis$UID
immune.4.analysis[is.na(immune.4.analysis)] <- 0
immune.4.analysis <- immune.4.analysis[ , -1]

# Create list of filtered variables (cytokines)
cytokines <- colnames(immune.4.analysis)
cytokines <- data.frame(cytokines) %>%
  dplyr::rename(Cytokines = 1)


## Pull out amplicon data and format for analysis
amp.4.analysis <- data.2[ , c(3, 19:748)]
amp.4.analysis <- data.frame(amp.4.analysis)
rownames(amp.4.analysis) <- amp.4.analysis$UID
amp.4.analysis <- amp.4.analysis[ , -1]

# Filter amplicon data for species with >X reads
amp.4.analysis <- amp.4.analysis %>% drop_na(Lactobacillus_iners)
amp.4.analysis[is.na(amp.4.analysis)] <- 0
amp.4.analysis <- amp.4.analysis[as.logical(rowSums(amp.4.analysis != 0)), ]

# Create list of filtered variables (species)
species <- colnames(amp.4.analysis)
species <- data.frame(species) %>%
  dplyr::rename(Species = 1)

## Join data frames for analysis
joined.data <- amp.4.analysis %>% merge(immune.4.analysis, by = 0)
rownames(joined.data) <- joined.data$Row.names
joined.data$Row.names <- NULL

# Perform log10 transformation
joined.data <- log2(joined.data + 1)


## Create data frame of every combination of metabolite and class with leading data frame name
combinations <- crossing(species, cytokines)
combinations.merge <- combinations %>% dplyr::rename(Species = 1,
                                                     Cytokines = 2)
combinations <- t(combinations)
combinations <- data.frame(combinations)
names(combinations) = str_sub(names(combinations), 2)

## Run the Spearman correlation test
spearman.output <- lapply(combinations, function(x) {
  (cor.test(joined.data[ ,x[1]], joined.data[ ,x[2]], 
            method = "spearman",
            exact = FALSE))
})

## Format spearman output list into data frame
spearman.output.df <- data.frame(matrix(unlist(spearman.output), nrow = length(spearman.output), byrow = TRUE))

spearman.output.df <- spearman.output.df[, c(1:5)]
spearman.output.df <- spearman.output.df %>%
  dplyr::rename(S.statistic = 1,
                pavalue = 2,
                rho = 3,
                null.value = 4,
                type = 5)

combinations.merge <- data.frame(rowname = rownames(combinations.merge), combinations.merge, row.names = NULL)
spearman.output.df <- data.frame(rowname = rownames(spearman.output.df), spearman.output.df, row.names = NULL)

spearman.output.df <- spearman.output.df %>% full_join(combinations.merge, by = "rowname")
spearman.output.df <- spearman.output.df[ , -1]
spearman.output.df$rho <- as.numeric(spearman.output.df$rho)
spearman.output.df$pavalue <- as.numeric(spearman.output.df$pavalue)

all.days.m.spp <- reshape2::melt(all.days[,1:748], id.vars=names(all.days)[1:18], variable.name="Species", value.name = "abd")
all.days.m.spp <- all.days.m.spp %>% group_by(Species, STATUS) %>% summarise(mean_abd=mean(abd, na.rm=T))

corr.results.df.melt<-merge(spearman.output.df, all.days.m.spp, all.x=TRUE)
corr.results.df.melt$corr<-corr.results.df.melt$rho
corr.results.df.melt$corr_dir<-factor(ifelse(corr.results.df.melt$corr > 0, "Positive Correlations", "Negative Correlations"), levels=c("Positive Correlations", "Negative Correlations"))
corr.results.df.melt$STATUS<-factor(corr.results.df.melt$STATUS, levels=c("Prior to SBV", "BVDX", "Early MET", "Late MET", "After MET"))
corr.results.df.melt<-corr.results.df.melt[!is.na(corr.results.df.melt$Species), ]
corr.results.df.melt$sig<-ifelse(corr.results.df.melt$pavalue < 0.05, "*", ifelse(corr.results.df.melt$pavalue < 0.01, "**", ifelse(corr.results.df.melt$pavalue < 0.001, "***", "")))
#corr.results.df.melt$sig<-paste0(corr.results.df.melt$Species, "\n", "rho = ", round(corr.results.df.melt$rho, 2))
corr.results.df.melt<-corr.results.df.melt[str_count(corr.results.df.melt$Species, "_") > 0, ]

the.cytos<-c("sEcad", "IP10", "MMP9")
to.plot<-corr.results.df.melt[corr.results.df.melt$Cytokines %in% the.cytos & corr.results.df.melt$pavalue < 0.001 & abs(corr.results.df.melt$corr) > 0.6, ]
#to.plot<-corr.results.df.melt[corr.results.df.melt$pavalue < 0.01 & abs(corr.results.df.melt$rho) > 0.5, ]
to.plot$Cytokines<-factor(to.plot$Cytokines, levels=c("sEcad", "IP10", "MMP9"))
to.plot<-to.plot[complete.cases(to.plot), ]
to.plot$Species<-gsub("_", " ", to.plot$Species)
to.plot$STATUS<-gsub(" ", "\n", to.plot$STATUS)
to.plot$STATUS<-factor(to.plot$STATUS, levels=c("Prior\nto\nSBV", "BVDX", "Early\nMET", "Late\nMET", "After\nMET"))
to.plot$STATUS<-factor(to.plot$STATUS, levels=c( "After\nMET","Late\nMET", "Early\nMET", "BVDX", "Prior\nto\nSBV"))
to.plot$abd<-ifelse(to.plot$rho < 0, log2(to.plot$mean_abd+1)*-1, log2(to.plot$mean_abd+1))
#to.plot$Species<-paste0(to.plot$Species, ", ρ = ", round(to.plot$rho, 2))
ggplot(to.plot, aes(x=STATUS, y=reorder(Species, abd), fill=abd))+
  geom_tile()+
  #geom_text(size=2, aes(label = round(rho, 2), x = 0.7, y = reorder(Species, abd)))+
  facet_wrap(~Cytokines, scales="free_x")+
  scale_fill_gradientn(colours = c("skyblue4", "skyblue3","lightblue", "powderblue", "lightblue1", "white", "lavenderblush", "pink", "pink2", "red", "darkred"),
                       name="Log2\nMean\nSpecies\nAbundance")+
  theme( 
    text=element_text(size=6, family = "Arial", color="black"), 
    axis.text.y=element_text(size=8, color="black", face="bold"), 
    axis.text.x=element_text(size=6, color="black", angle=45, hjust=1, face="italic"), 
    legend.position="right", 
    legend.title.position = "top", 
    legend.title=element_text(hjust=0.5),
    line=element_line(linewidth=0.1), 
    strip.text = element_text(size=8, face="bold"), 
    strip.background = element_rect(linewidth = 0.2), 
    plot.background = element_rect(linewidth=0.2))+
  coord_flip()+
  #ggtitle("Species/Immune Marker Correlations (|ρ| > 0.6, p<0.05)\nRed Scale Abundances: Positive Correlation with Marker\nBlue Scale Abundances: Negative Correlation with Marker")+
  ylab("")+
  xlab("")
ggsave("all_corr.png", height=4, width=15, dpi = 600)
