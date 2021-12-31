# Lakshmi Yechuri 
# 12110035

# importing required packages
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(dendextend)
library(dplyr)
library(corrplot)
library(dummies)
library(psych)
library(writexl)
library(factoextra)


# reading the dataset and inserting column names
winedf = read.csv("wine-data.csv", header = FALSE)
colnames(winedf) <- c("Type","Alcohol","Malic acid","Ash","Alcalinity of ash", 
                      "Magnesium", "Total phenols", "Flavanoids", 
                      "Nonflavanoid phenols", "Proanthocyanins", 
                      "Color intensity", "Hue", "OD280/OD315 of diluted wines", 
                      "Proline")
winedf$TypeC <- winedf$Type
winedf$TypeC <- factor(winedf$TypeC)

# creating dummy variables for wine types
winedf$Type <- factor(winedf$Type)
winedf <- dummy.data.frame(winedf, names = "Type",omit.constants=FALSE)
write.csv(summary(winedf), "winedf_summary.csv")


# creating a scatterplot to understand if any variables have strong correlation 
# if there is any need for PCA
pairs.panels(winedf[,4:16], method = "pearson", gap=0,
             bg=c("cyan","lavender","yellow")[winedf$TypeC],
             pch=21)

# PCA using prcomp function
wine.pca <- prcomp(winedf[,4:16],center = TRUE,scale. = TRUE)

summary(wine.pca)
pca.rotation <- data.frame(wine.pca$rotation)
write_xlsx(pca.rotation,"Wine data Loading.xlsx")

var.pca <- get_pca_var(wine.pca)
corrplot(var.pca$cos2, is.corr = FALSE, addrect = 2)

pc.scores <- wine.pca$x

winescores <- data.frame(cbind(winedf[,17],pc.scores))

df_grp_region <- winescores %>% group_by(V1)  %>%
  summarise('Health and Flavor Factor'  = mean(PC1),
            'Alcohol and color factor' = mean(PC2),
            'Alkalinity factor ' = mean(PC3),
            'Sweetness and Hue factor ' = mean(PC4),
            'Mineral and Aroma factor' = mean(PC5),
            .groups = 'drop')

write_xlsx(df_grp_region,"PCs scores summary .xlsx")

# Creating a biplot
g <- ggbiplot(wine.pca, obs.scale = 1,var.scale = 1, groups = winedf$TypeC,
              ellipse = TRUE, circle = TRUE, ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

# scale date for further analysis
standardized_data <- scale(winedf[,1:16])

# hierarchical clustering on all chemical measurements
d <- dist(standardized_data, method = "euclidean")
result.hc <- hclust(d, method = "ward.D2" )
result.hc <- as.dendrogram(result.hc)
myDendogram = color_branches(result.hc,k=3) 
plot(myDendogram, main = "H-clust with all chemical measurements")
finalclusters <- cutree(result.hc, k=3) 
aggregates.h.cm <- aggregate(winedf[,1:16],list(finalclusters),mean)
aggregates.h.cm.df = data.frame(Cluster=aggregates.h.cm[,1],
                                Freq=as.vector(table(finalclusters)),
                                aggregates.h.cm[,-1])
write_xlsx(aggregates.h.cm.df,"All Chemical H-clust aggregates.xlsx")

# hierarchical clustering on two most significant PC scores
dPC <- dist(pc.scores[,1:2], method = "euclidean")
result.hcPC <- hclust(dPC, method = "ward.D2" )
result.hcPC <- as.dendrogram(result.hcPC)
myDendogramPC = color_branches(result.hcPC,k=3) 
plot(myDendogramPC, main = "H-Clust with two most significant PC scores")
finalclustersPC <- cutree(result.hcPC, k=3) 
table(finalclustersPC)
aggregates.h.pc <- aggregate(pc.scores[,1:2],list(finalclusters),mean)
aggregates.h.pc.df = data.frame(Cluster=aggregates.h.pc[,1],Freq=as.vector(table(finalclustersPC)),aggregates.h.pc[,-1])
write_xlsx(aggregates.h.cm.df,"Top 2 PC Hirearchial cluster aggregates.xlsx")



# K means on all chemical measurements

# Determine number of clusters
clustercount <- matrix(nrow=10, ncol=1)
for (i in 1:10) clustercount[i] <- kmeans(standardized_data, 
                                          centers=i)$tot.withinss
plot(1:10, clustercount, type="b", main = "Elbow plot with all chemical measurements", 
     xlab="Number of clusters", ylab="Within groups sum of squares")

# k means fit
fit.all <- kmeans(standardized_data, centers=3, iter.max=10, nstart=10)

# computing cluster aggregates
aggregates.km <- aggregate(winedf[,1:16],list(fit.all$cluster),mean)
aggregates.km.df <- data.frame(Cluster=aggregates.km[,1],
                           Freq=as.vector(table(fit.all$cluster)),aggregates.km[,-1])
write_xlsx(aggregates.km.df,"All checmical K clust aggregates.xlsx")

# displaying clusters
fviz_cluster(fit.all, data = standardized_data,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",ellipse.type = "convex", 
             ggtheme = theme_bw(),
             main = "K-Clust with all chemical measurements ")



# k means using two most significant PC scores

## Determine number of clusters
pca.clustercount <- matrix(nrow=10, ncol=1)
for (i in 1:10) pca.clustercount[i] <- kmeans(pc.scores[,1:2], 
                                              centers=i)$tot.withinss
plot(1:10, pca.clustercount, type="b", main = "Elbow plot for PCA scores data",
     xlab="Number of clusters", ylab="Within groups sum of squares")

# k means fit
fit.pca <- kmeans(pc.scores[,1:2], centers=3, iter.max=10, nstart=10)
df.pca <- data.frame(fit.pca$centers)

# computing cluster aggregates
pca.aggregates <- aggregate(pc.scores[,1:2],list(fit.pca$cluster),mean)
pca.aggregates.df <- data.frame(Cluster=pca.aggregates[,1],Freq=as.vector(table(fit.pca$cluster)),pca.aggregates[,-1])
write_xlsx(pca.aggregates.df,"Top two PCA K clust aggregates.xlsx")

# displaying clusters

fviz_cluster(fit.pca, data = pc.scores[,1:2],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800" ), 
             geom = c("point"),
             ellipse.type = "convex", 
             ggtheme = theme_bw(),
             main = "K-Clust with two most significant PC scores ")

