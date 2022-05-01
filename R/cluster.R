library("cluster")
library("factoextra")
library("magrittr")

# Load  and prepare the data
data("USArrests")

my_data <- USArrests %>%
  na.omit() %>%          # Remove missing values (NA)
  scale()                # Scale variables

my_data = read.csv("/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/NMF_heatmap_Dong_data_dropna.csv",row.names = 1)
# View the firt 3 rows
head(my_data, n = 3)
# my_data[is.na(my_data)] = 0
my_data = t(my_data)

df = my_data
# Elbow method
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
## ---------------- Compute and visualize k-means clustering ----------
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
km.res <- kmeans(my_data, 3, nstart = 25)
# Visualize
library("factoextra")
fviz_cluster(km.res, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

km_cluster = km.res[["cluster"]]
write.csv(km_cluster,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/NMF_heatmap_cluster_group_km.csv")
# res.dist <- get_dist(USArrests, stand = TRUE, method = "pearson")

fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

## ------------ Partitioning clustering -------------------
fviz_nbclust(my_data, kmeans, method = "gap_stat")
# Compute and visualize k-means clustering
set.seed(123)
km.res <- kmeans(my_data, 3, nstart = 25)
# Visualize
library("factoextra")
fviz_cluster(km.res, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

# Compute PAM
library("cluster")
pam.res <- pam(my_data, 3)
# Visualize
fviz_cluster(pam.res)

pam_cluster = pam.res[["clustering"]]
write.csv(pam_cluster, "/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/NMF_heatmap_cluster_group_pam.csv")
## ------------ Hierarchical clustering -------------------

# Compute hierarchical clustering
res.hc <- my_data %>%          # USArrests
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering


# Visualize using factoextra
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

# Determining the optimal number of clusters
set.seed(123)

# Compute
library("NbClust")
res.nbclust <- my_data %>%
  scale() %>%
  NbClust(distance = "euclidean",
          min.nc = 2, max.nc = 10, 
          method = "complete", index ="all") 


# fviz_nbclust(my_data, kmeans, method = "gap_stat")
fviz_nbclust(res.nbclust, ggtheme = theme_minimal())

# Clustering validation statistics

set.seed(123)
# Enhanced hierarchical clustering, cut in 3 groups
res.hc <- my_data %>%
  scale() %>%
  eclust("hclust", k = 3, graph = T)

# Visualize with factoextra
fviz_dend(res.hc, palette = "jco",
          rect = TRUE, show_labels = T)

# Inspect the silhouette plot:
fviz_silhouette(res.hc)

# Which samples have negative silhouette? To what cluster are they closer?
# Silhouette width of observations
sil <- res.hc$silinfo$widths[, 1:3]
write.csv(sil,"/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/response_R2C5/NMF_heatmap_cluster_group.csv")
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]




