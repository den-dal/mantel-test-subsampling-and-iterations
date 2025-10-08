# mantel-test-subsampling-and-iterations
````
install.packages("geosphere")
library(vegan)
library(ggplot2)
library(geosphere)
library(dplyr)
setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/rbcL")
list.files()
df = read.csv("OTU_lat_long.csv", header= TRUE)

# === Filter to the island of interest ===
df <- subset(df, Island == "MARCHENA")

# === List locations and number of samples for the island ===
location_counts <- as.data.frame(table(df$Location))
colnames(location_counts) <- c("Location", "N_samples")
print(location_counts)

min_sam <- 8

# === Keep only locations with >= # of samples ===
df <- df[df$Location %in% names(which(table(df$Location) >= min_sam)), ]
 
####Function: one subsampling run ===
  mantel_subsample <- function(df, n_sub = min_sam) {
    sub_df <- df %>%
      group_by(Location) %>%
      slice_sample(n = n_sub) %>%
      ungroup()
    
    abund = sub_df[,9:ncol(sub_df)]
    geo = data.frame(sub_df$Longitude, sub_df$Latitude)
    
    dist.abund = vegdist(abund, method = "bray")
    d.geo = distm(geo, fun = distHaversine)
    dist.geo = as.dist(d.geo)
    
    m <- mantel(dist.abund, dist.geo, method = "spearman",
                permutations = 9999, na.rm = TRUE)
    return(c(r = m$statistic, p = m$signif))
  }

# === Repeat subsampling ===
set.seed(123)
results <- replicate(100, mantel_subsample(df), simplify = TRUE)
results <- t(results)
results <- as.data.frame(results)

# Summary of results
summary(results$r)
mean(results$r)

sd(results$r)
mean(results$p < 0.05)   # proportion of significant runs
````
