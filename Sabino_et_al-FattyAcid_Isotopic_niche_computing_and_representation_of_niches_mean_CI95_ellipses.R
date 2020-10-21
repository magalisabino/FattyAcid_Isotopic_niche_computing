##--------------------------------------------------------------------------------------------------------
## SCRIPT : Characterisation of trophic overlap between three sympatric species of spiny lobsters in Seychelles
## Specific content : - Computation of fatty acid and isotopic niches
##                    - Calculation of niche metrics (i.e. niche size and probabilities of niche overlap)
##                    - Representation (plot) of niche mean ellipses and CI95% ellipses
## As part of :
##    Sabino et al. "Habitat degradation increases interspecific competition between three spiny lobsters in Seychelles"
##
## Author : Magali Sabino
## Last update : 2020-10-21
##
## For more information on the nicheROVER package, see:
##    - Swanson et al. 2015 "A new probabilistic method for quantifying n-dimensional ecological niches and niche overlap" Ecology, 96 (2), p. 318-324
##    - Lysy et al. 2017 "nicheROVER: (Niche) (R)egion and Niche (Over)lap Metrics for Multidimensional Ecological Niches" Web: https://cran.r-project.org/web/packages/nicheROVER/index.html
##    - An Ecologist's Guide to nicheROVER: Niche Region and Niche Overlap Metrics for Multidimensional Ecological Niches Web: https://cran.r-project.org/web/packages/nicheROVER/vignettes/ecol-vignette.html
##
####
##
## R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
## Copyright (C) 2018 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##-------------------------------------------------.-------------------------------------------------------

### 0 // Packages ##########################################################################################

lapply(c("tidyverse", "vegan", "nicheROVER"),
       library, character.only = TRUE)


### I // Data ##############################################################################################

## 1 / Dataset creation

# The dataset used in the paper associated to this script is available from the corresponding author on
# reasonable request. Thus, here we generate type datasets for each species containing fatty acid profiles
# (expressed as % total fatty acid) and stable isotope values using priors on distribution (i.e. mean and
# standard deviation). In this dataset, each line corresponds to one individual while fatty acid proportions
# and d13C and d15N values are in columns.

# We create a function to generate fatty acid proportions and to ensure that the sum of all fatty acids for
# one individual is equal to 100%.
Sample_prop_matrix <- 
  function(N_samp = 100, y_cols = 20, Mean_priors, SD_priors, Col_names, limit_inf = 0){
    Output_matrix <- matrix(nrow = N_samp, ncol = y_cols)
    colnames(Output_matrix) <- Col_names
    for (u in 1:N_samp){
      for (i in 1:y_cols){
        generated_val <- rnorm(n=1, mean=Mean_priors[i], sd=SD_priors[i])
        generated_val <- ifelse(generated_val < limit_inf, 0, generated_val)
        Output_matrix[u, i] <- generated_val
      }
    }
    rowSums_matrix <- rowSums(Output_matrix)
    Output_matrix <- Output_matrix/rowSums_matrix
    result <- Output_matrix * 100
  }

# We generate a dataset for each species to have different priors between species
SP1 <- as.data.frame(Sample_prop_matrix(N_samp = 40, y_cols = 20,
                               Mean_priors = c(0.3, 1.1, 4.4, 2.1, 1.3,
                                               1.5, 10.5, 12, 2.8, 1.2,
                                               1, 0.3, 8.7, 7.1, 0.2,
                                               0.5, 0.4, 15.2, 3.7, 2.5),
                               SD_priors = c(0.1, 0.3, 1.3, 0.3, 0.2,
                                             0.5, 2.5, 2, 0.9, 0.5,
                                             0.5, 0.1, 1.1, 2, 0.1,
                                             0.2, 0.1, 1.9, 1.3, 0.6),
                               Col_names = c("FA1", "FA2", "FA3", "FA4", "FA5",
                                             "FA6","FA7","FA8","FA9","FA10",
                                             "FA11","FA12","FA13","FA14","FA15",
                                             "FA16","FA17","FA18","FA19","FA20"),
                               limit_inf = 0)) %>% 
  mutate(species = "SP1",
         d13C = rnorm(n = 40, mean = -13.5, sd = 0.7),
         d15N = rnorm(n = 40, mean = 11.9, sd = 0.4))

SP2 <- as.data.frame(Sample_prop_matrix(N_samp = 40, y_cols = 20,
                          Mean_priors = c(0.4, 1.2, 6.2, 2, 1,
                                          1.6, 7.1, 8.9, 1.2, 1,
                                          0.9, 0.3, 8.5, 7.5, 0.3,
                                          0.6, 0.5, 15, 3.5, 2.3),
                          SD_priors = c(0.1, 0.2, 1.5, 0.4, 0.3,
                                        0.5, 1.4, 1, 0.2, 0.5,
                                        0.4, 0.1, 1.0, 1.9, 0.1,
                                        0.2, 0.1, 2, 1, 1),
                          Col_names = c("FA1", "FA2", "FA3", "FA4", "FA5",
                                        "FA6","FA7","FA8","FA9","FA10",
                                        "FA11","FA12","FA13","FA14","FA15",
                                        "FA16","FA17","FA18","FA19","FA20"),
                          limit_inf = 0)) %>% 
  mutate(species = "SP2",
         d13C = rnorm(n = 40, mean = -13, sd = 0.7),
         d15N = rnorm(n = 40, mean = 11, sd = 0.3))

SP3 <- as.data.frame(Sample_prop_matrix(N_samp = 40, y_cols = 20,
                          Mean_priors = c(0.3, 1.2, 6.5, 1.9, 1.5,
                                          4, 7.6, 8.5, 1, 1.1,
                                          1, 0.4, 11.5, 7.5, 0.2,
                                          0.6, 0.4, 15, 3.8, 2.7),
                          SD_priors = c(0.1, 0.3, 1.5, 0.7, 0.4,
                                        1.5, 2, 1.2, 0.2, 0.3,
                                        0.5, 0.1, 1.5, 2, 0.1,
                                        0.1, 0.1, 1.5, 1.5, 0.3),
                          Col_names = c("FA1", "FA2", "FA3", "FA4", "FA5",
                                        "FA6","FA7","FA8","FA9","FA10",
                                        "FA11","FA12","FA13","FA14","FA15",
                                        "FA16","FA17","FA18","FA19","FA20"),
                          limit_inf = 0)) %>% 
  mutate(species = "SP3",
         d13C = rnorm(n = 40, mean = -12, sd = 0.7),
         d15N = rnorm(n = 40, mean = 9.9, sd = 0.3))

# We merge all created datasets to have one big dataset with all generated data
trophic_data <- rbind(SP1,SP2,SP3)
rm(SP1,SP2,SP3)


## 2 / Fatty acid data cleasing - Get only proportions > 0.8% for data treatment

FAmean <- colMeans(trophic_data[,1:20], na.rm = TRUE) # Calculating mean for each fatty acid
others <- names(which(FAmean <0.8)) # Give names of fatty acids which are < 0.8%
trophic_data_clean <- trophic_data[,!names(trophic_data) %in% (others)] # Get all fatty acids > 0.8%
trophic_data_clean <- droplevels(trophic_data_clean)
trophic_data_clean$others <- rowSums(trophic_data[,others]) # Calculating sum of all fatty acids < 0.8%
rm(FAmean,others)


### II // Fatty acid niches ################################################################################

## 1 / Non-metric multidimensional scaling (nMDS)

# Fatty acid profiles contain a high number of dimensions (here, for fatty acids > 0.8%, n = 15), thus we
# used nMDS ordination with a Bray-Curtis dissimilarity matrix to constrain all 15 dimensions in only two.

names(trophic_data_clean)
MDS_FA <- metaMDS(trophic_data_clean[,1:15],distance = "bray",k = 2,try = 300)

data_FAniche <- as.data.frame(scores(MDS_FA)) # Using the scores function from vegan to extract the sitescores and convert to a data.frame
data_FAniche$species <- trophic_data_clean$species # Add species for future niche computing

## 2 / Niche ellipses computing

# The fish.par dataframe contains the posterior distribution of mu and sigma, i.e. the parameters of all
# computed ellipses

nsamples <- 1000
fish.par <- tapply(1:nrow(data_FAniche), data_FAniche$species,
                   function(ii) niw.post(nsamples = nsamples, X = data_FAniche[ii,1:2]))

## 3 / Niche metrics

# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})
# Niche size calculation: mean of all ellipse sizes and CI95%
rbind(est = colMeans(fish.size),
      quant = apply(fish.size, 2, quantile, prob = c(0.025,0.975)))

# Niche overlaps
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95,0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs = c(0.025, 0.975), na.rm=T) * 100

rm(nsamples, over.stat, over.mean, over.ICinf)

## 4 / Extraction of points coordinates of all calculated ellipses

# For this section, the niche.plot function of the nicheROVER package was modified
# in order to get access to all coordinates of all computed ellipses and for
# each species.

niche.par = fish.par
pfrac = 0.05
alpha=0.95

niso <- ncol(niche.par[[1]]$mu)
nspec <- length(niche.par)
npts <- 100
nell <- sapply(niche.par, function(x) nrow(x$mu))
species.names <- colnames(fish.size)
D <- combn(niso, 2)

Mu_all <- NULL
Sigma_all <- NULL
ell.coord_all <- NULL

for (ii in 1:nspec) {
  
  #ell.tmp <- array(NA, c(nell[ii], ncol(D), npts + 1, 2))
  for (jj in 1:nell[ii]) {
    
    for (kk in 1:ncol(D)) {
      
      Mu_prov <- as.data.frame(niche.par[[ii]]$mu[jj, D[, 
                                                        kk]])
      Mu_prov <- as.data.frame(t(Mu_prov))
      Mu_prov$SPP <- species.names[ii]
      Mu_prov$Iter <- jj
      
      Mu_all <- rbind(Mu_all, Mu_prov)
      
      
      
      Sigma_prov <- as.data.frame(niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                   kk], jj])
      
      Sigma_prov$SPP <- species.names[ii]
      Sigma_prov$Iter <- jj
      
      Sigma_all <- rbind(Sigma_all, Sigma_prov)
      
      ell.coord <- as.data.frame(ellipse(niche.par[[ii]]$mu[jj, D[, 
                                                                  kk]], V = niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                                                             kk], jj], alpha = alpha, n = npts))
      
      ell.coord_prov <- ell.coord
      ell.coord_prov$N_point <- c(1:101)
      ell.coord_prov$SPP <- species.names[ii]
      ell.coord_prov$Iter <- jj
      
      ell.coord_all <- rbind(ell.coord_all, ell.coord_prov)
    }
  }
}

rm(ii,jj,kk,alpha,niso)
rm(ell.coord)
rm(pfrac,nspec,npts,nell,species.names,D,Mu_all,niche.par,
   Sigma_all,Mu_prov,Sigma_prov,ell.coord_prov)


## 5 / Recovery of point coordinates for CI95% ellipses

ell.coord_inf_all <- NULL
ell.coord_sup_all <- NULL

for (i in 1:ncol(fish.size)){
  
  # Selection of ellipse sizes for only one species
  fish_selec <- fish.size[,i]
  fish_area <- data.frame(area = fish_selec, Iter = c(1:1000))
  
  
  # Selection of point coordinates of all ellipses for one species
  fish_name <- colnames(fish.size)[i]
  ell.coord_fish <- ell.coord_all[ell.coord_all$SPP == fish_name,]
  
  
  # Merging ellipse sizes with data related to ellipses
  ell.coord_fish <- ell.coord_fish %>%
    left_join(fish_area, by = "Iter")
  
  
  # Identification of 2.5% smallest and 2.5% largest ellipses to get point coordinates
  Middle_area <- fish_area[order(fish_area$area, decreasing = F),] # Ranking ellipse sizes from smallest to largest
  QUANT_inf_area <- Middle_area[0:250,] # Selection of size of 2.5% smallest ellipses
  QUANT_sup_area <- Middle_area[750:1000,] # Selection of size of 2.5% largest ellipses
  
  ell.coord_fish$CI_inf <- ifelse(ell.coord_fish$area %in% QUANT_inf_area[,1], 1, 0)
  ell.coord_fish$CI_sup <- ifelse(ell.coord_fish$area %in% QUANT_sup_area[,1], 1, 0)
  
  ell.coord_inf <- ell.coord_fish[ell.coord_fish$CI_inf == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  ell.coord_sup <- ell.coord_fish[ell.coord_fish$CI_sup == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  
  ell.coord_inf_all <- rbind(ell.coord_inf_all,ell.coord_inf)
  ell.coord_sup_all <- rbind(ell.coord_sup_all,ell.coord_sup)
}

rm(i)
rm(fish_selec, fish_area, fish_name, ell.coord_fish,
   Middle_area, QUANT_inf_area, QUANT_sup_area,
   ell.coord_inf, ell.coord_sup)
rm(fish.par, fish.size)

## 6 / Calculation of mean and CI95% ellipses for each species

# Mean ellipses: mean coordinates of all ellipses for each species
ell.coord_MEAN <- as.data.frame(
  ell.coord_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

# CI95% ellipses: mean coordinates of 2.5% smallest ellipses for each species
ell.coord_CIinf <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

# CI95% ellipses: mean coordinates of 2.5% largest ellipses for each species
ell.coord_CIsup <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)

## 7 / Recovery of nMDS data for further plotting

# Getting fatty acid nMDS coordinates to plot them on the niches plot
FA.scores <- as.data.frame(scores(MDS_FA, "species"))  # Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores$FA <- rownames(FA.scores) # To add species names = fatty acid names
head(FA.scores)

rm(MDS_FA)

## 8 / Plotting niches

data_FAniche %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = species), size = 2.5, shape = 21, color = "grey50") + # add all individual points
  geom_text(data = FA.scores, aes(x = NMDS1, y = NMDS2,label = FA), size = 3, fontface = "bold", color = "grey23") + # add the FA labels
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species")+
  theme(axis.title.x = element_text(face = "italic"), # remove x-axis labels
        axis.title.y = element_text(face = "italic"), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank())

rm(data_FAniche, ell.coord_MEAN, ell.coord_CIinf,
   ell.coord_CIsup, FA.scores)

### III // Isotopic niches ####

## 1 / Isotopic data

data_SIniche <- trophic_data %>% 
  select(species,d13C,d15N)

## 2 / Niche ellipses computing

# The fish.par dataframe contains the posterior distribution of mu and sigma, i.e. the parameters of all
# computed ellipses

nsamples <- 1000
fish.par <- tapply(1:nrow(data_SIniche), data_SIniche$species,
                   function(ii) niw.post(nsamples = nsamples, X = data_SIniche[ii,2:3]))

## 3 / Niche metrics

# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})
# Niche size calculation: mean of all ellipse sizes and CI95%
rbind(est = colMeans(fish.size),
      quant = apply(fish.size, 2, quantile, prob = c(0.025,0.975)))

# Niche overlaps
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95,0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

rm(nsamples, over.stat, over.mean, over.ICinf)

## 4 / Extraction of points coordinates of all calculated ellipses

# For this section, the niche.plot function of the nicheROVER package was modified
# in order to get access to all coordinates of all computed ellipses and for
# each species.

niche.par = fish.par
pfrac = 0.05
alpha=0.95

niso <- ncol(niche.par[[1]]$mu)
nspec <- length(niche.par)
npts <- 100
nell <- sapply(niche.par, function(x) nrow(x$mu))
species.names <- colnames(fish.size)
D <- combn(niso, 2)

Mu_all <- NULL
Sigma_all <- NULL
ell.coord_all <- NULL

for (ii in 1:nspec) {
  
  #ell.tmp <- array(NA, c(nell[ii], ncol(D), npts + 1, 2))
  for (jj in 1:nell[ii]) {
    
    for (kk in 1:ncol(D)) {
      
      Mu_prov <- as.data.frame(niche.par[[ii]]$mu[jj, D[, 
                                                        kk]])
      Mu_prov <- as.data.frame(t(Mu_prov))
      Mu_prov$SPP <- species.names[ii]
      Mu_prov$Iter <- jj
      
      Mu_all <- rbind(Mu_all, Mu_prov)
      
      
      
      Sigma_prov <- as.data.frame(niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                   kk], jj])
      
      Sigma_prov$SPP <- species.names[ii]
      Sigma_prov$Iter <- jj
      
      Sigma_all <- rbind(Sigma_all, Sigma_prov)
      
      ell.coord <- as.data.frame(ellipse(niche.par[[ii]]$mu[jj, D[, 
                                                                  kk]], V = niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                                                             kk], jj], alpha = alpha, n = npts))
      
      ell.coord_prov <- ell.coord
      ell.coord_prov$N_point <- c(1:101)
      ell.coord_prov$SPP <- species.names[ii]
      ell.coord_prov$Iter <- jj
      
      ell.coord_all <- rbind(ell.coord_all, ell.coord_prov)
    }
  }
}

rm(ii,jj,kk,alpha,niso)
rm(ell.coord)
rm(pfrac,nspec,npts,nell,species.names,D,Mu_all,niche.par,
   Sigma_all,Mu_prov,Sigma_prov,ell.coord_prov)


## 5 / Recovery of point coordinates for CI95% ellipses

ell.coord_inf_all <- NULL
ell.coord_sup_all <- NULL

for (i in 1:ncol(fish.size)){
  
  # Selection of ellipse sizes for only one species
  fish_selec <- fish.size[,i]
  fish_area <- data.frame(area = fish_selec, Iter = c(1:1000))
  
  
  # Selection of point coordinates of all ellipses for one species
  fish_name <- colnames(fish.size)[i]
  ell.coord_fish <- ell.coord_all[ell.coord_all$SPP == fish_name,]
  
  
  # Merging ellipse sizes with data related to ellipses
  ell.coord_fish <- ell.coord_fish %>%
    left_join(fish_area, by = "Iter")
  
  
  # Identification of 2.5% smallest and 2.5% largest ellipses to get point coordinates
  Middle_area <- fish_area[order(fish_area$area, decreasing = F),] # Ranking ellipse sizes from smallest to largest
  QUANT_inf_area <- Middle_area[0:250,] # Selection of size of 2.5% smallest ellipses
  QUANT_sup_area <- Middle_area[750:1000,] # Selection of size of 2.5% largest ellipses
  
  ell.coord_fish$CI_inf <- ifelse(ell.coord_fish$area %in% QUANT_inf_area[,1], 1, 0)
  ell.coord_fish$CI_sup <- ifelse(ell.coord_fish$area %in% QUANT_sup_area[,1], 1, 0)
  
  ell.coord_inf <- ell.coord_fish[ell.coord_fish$CI_inf == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  ell.coord_sup <- ell.coord_fish[ell.coord_fish$CI_sup == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  
  ell.coord_inf_all <- rbind(ell.coord_inf_all,ell.coord_inf)
  ell.coord_sup_all <- rbind(ell.coord_sup_all,ell.coord_sup)
}

rm(i)
rm(fish_selec, fish_area, fish_name, ell.coord_fish,
   Middle_area, QUANT_inf_area, QUANT_sup_area,
   ell.coord_inf, ell.coord_sup)
rm(fish.par, fish.size)

## 6 / Calculation of mean and CI95% ellipses for each species

# Mean ellipses: mean coordinates of all ellipses for each species
ell.coord_MEAN <- as.data.frame(
  ell.coord_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

# CI95% ellipses: mean coordinates of 2.5% smallest ellipses for each species
ell.coord_CIinf <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

# CI95% ellipses: mean coordinates of 2.5% largest ellipses for each species
ell.coord_CIsup <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_point) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
)

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)

## 6 / Plotting niches

data_SIniche %>% 
  group_by(species) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE)/sqrt(length(d13C)),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)/sqrt(length(d15N))) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = species), size = 2, shape = 21, color = "grey50") +
  geom_polygon(data = ell.coord_MEAN, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.5)+
  geom_polygon(data = ell.coord_CIinf, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.9)+
  geom_polygon(data = ell.coord_CIsup, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.9)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = species), shape = 23, color = 'grey50', size=5) +
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species")+
  theme(axis.title.x = element_text(face = "italic"), # remove x-axis labels
        axis.title.y = element_text(face = "italic"), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank())

rm(data_SIniche, ell.coord_MEAN, ell.coord_CIinf,
   ell.coord_CIsup)