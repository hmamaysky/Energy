setwd("C:/Users/j1kxb09/OneDrive - FR Banks/Documents/MNSC meetings/data")
topic <- read.csv("topicVecs (3).csv", stringsAsFactors = FALSE)
topic <- subset(topic, select = -X)
#install.packages("ggpubr")
#install.packages("factoextra")
#install.packages("NbClust")
#install.packages("gridExtra")
#install.packages("wordcloud2")
library(ggpubr)
library(factoextra)
library(NbClust)
library(tm)
library(dplyr)
library(wordcloud)
library(wordcloud2)
library(tidyr)
library(lubridate)
library(ggplot2)
library("cowplot")
library(gridExtra)
######## cluster plot #######
set.seed(12345)
res.km <- kmeans(scale(topic), 7, nstart = 25)
res.km$cluster

fviz_cluster(res.km, data = topic,
             palette = c("purple", "#2E9FDF", "#E7B800","#fc0439","navy", "green", "black"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
#1) extracting cluster assignments
cluster_assignments <- res.km$cluster

clustered_data <- split(topic, cluster_assignments)
cluster1 <- clustered_data[[1]]
cluster2 <- clustered_data[[2]]
cluster3 <- clustered_data[[3]]
cluster4 <- clustered_data[[4]]
cluster5 <- clustered_data[[5]]
cluster6 <- clustered_data[[6]]
cluster7 <- clustered_data[[7]]
###cluster 1###
cluster1_averages <- colMeans(cluster1)
cluster1_averages <- cluster1_averages[order(-cluster1_averages)]
cluster1_averages_df <- data.frame(Variable = names(cluster1_averages),
                                   Average = cluster1_averages)

top_30_words <- head(cluster1_averages_df[order(-
                                                  cluster1_averages_df$Average), ],30)

cl1 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5))

###cluster 2 ###

cluster2_averages <- colMeans(cluster2)
cluster2_averages <- cluster2_averages[order(-cluster2_averages)]
cluster2_averages_df <- data.frame(Variable = names(cluster2_averages),
                                   Average = cluster2_averages)

top_30_words <- head(cluster2_averages_df[order(-
                                                  cluster2_averages_df$Average), ],30)

cl2 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )

## cluster 3 ###
cluster3_averages <- colMeans(cluster3)
cluster3_averages <- cluster3_averages[order(-cluster3_averages)]
cluster3_averages_df <- data.frame(Variable = names(cluster3_averages),
                                   Average = cluster3_averages)

top_30_words <- head(cluster3_averages_df[order(-
                                                  cluster3_averages_df$Average), ],30)

cl3 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5))

####cluster4###

cluster4_averages <- colMeans(cluster4)
cluster4_averages <- cluster4_averages[order(-cluster4_averages)]
cluster4_averages_df <- data.frame(Variable = names(cluster4_averages),
                                   Average = cluster4_averages)

top_30_words <- head(cluster4_averages_df[order(-
                                                  cluster4_averages_df$Average), ],30)

cl4 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )


####cluster5###

cluster5_averages <- colMeans(cluster5)
cluster5_averages <- cluster5_averages[order(-cluster5_averages)]
cluster5_averages_df <- data.frame(Variable = names(cluster5_averages),
                                   Average = cluster5_averages)

top_30_words <- head(cluster5_averages_df[order(-
                                                  cluster5_averages_df$Average), ],30)

cl5 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )


####cluster6###

cluster6_averages <- colMeans(cluster6)
cluster6_averages <- cluster6_averages[order(-cluster6_averages)]
cluster6_averages_df <- data.frame(Variable = names(cluster6_averages),
                                   Average = cluster6_averages)

top_30_words <- head(cluster6_averages_df[order(-
                                                  cluster6_averages_df$Average), ],30)

cl6 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )

####cluster7###

cluster7_averages <- colMeans(cluster7)
cluster7_averages <- cluster7_averages[order(-cluster7_averages)]
cluster7_averages_df <- data.frame(Variable = names(cluster7_averages),
                                   Average = cluster7_averages)

top_30_words <- head(cluster7_averages_df[order(-
                                                  cluster7_averages_df$Average), ],30)

cl7 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )

###### 5 Clusters #######
set.seed(12345)
res.km <- kmeans(scale(topic), 5, nstart = 25)
res.km$cluster

fviz_cluster(res.km, data = topic,
             palette = c("#2E9FDF", "#E7B800","#fc0439","purple", "green"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#1) extracting cluster assignments
cluster_assignments <- res.km$cluster

clustered_data <- split(topic, cluster_assignments)
cluster1 <- clustered_data[[1]]
cluster2 <- clustered_data[[2]]
cluster3 <- clustered_data[[3]]
cluster4 <- clustered_data[[4]]
cluster5 <- clustered_data[[5]]

###cluster 1###
cluster1_averages <- colMeans(cluster1)
cluster1_averages <- cluster1_averages[order(-cluster1_averages)]
cluster1_averages_df <- data.frame(Variable = names(cluster1_averages),
                                   Average = cluster1_averages)

top_30_words <- head(cluster1_averages_df[order(-
                                                  cluster1_averages_df$Average), ],30)

cl1 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )

###cluster 2 ###

cluster2_averages <- colMeans(cluster2)
cluster2_averages <- cluster2_averages[order(-cluster2_averages)]
cluster2_averages_df <- data.frame(Variable = names(cluster2_averages),
                                   Average = cluster2_averages)

top_30_words <- head(cluster2_averages_df[order(-
                                                  cluster2_averages_df$Average), ],30)

cl2 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )

## cluster 3 ###
cluster3_averages <- colMeans(cluster3)
cluster3_averages <- cluster3_averages[order(-cluster3_averages)]
cluster3_averages_df <- data.frame(Variable = names(cluster3_averages),
                                   Average = cluster3_averages)

top_30_words <- head(cluster3_averages_df[order(-
                                                  cluster3_averages_df$Average), ],30)

cl3 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5))

####cluster4###

cluster4_averages <- colMeans(cluster4)
cluster4_averages <- cluster4_averages[order(-cluster4_averages)]
cluster4_averages_df <- data.frame(Variable = names(cluster4_averages),
                                   Average = cluster4_averages)

top_30_words <- head(cluster4_averages_df[order(-
                                                  cluster4_averages_df$Average), ],30)

cl4 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )


####cluster5###

cluster5_averages <- colMeans(cluster5)
cluster5_averages <- cluster5_averages[order(-cluster5_averages)]
cluster5_averages_df <- data.frame(Variable = names(cluster5_averages),
                                   Average = cluster5_averages)

top_30_words <- head(cluster5_averages_df[order(-
                                                  cluster5_averages_df$Average), ],30)

cl5 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5),
                 random.color = )





###### 6 Clusters #######
set.seed(12345)
res.km <- kmeans(scale(topic), 6, nstart = 25)
res.km$cluster

fviz_cluster(res.km, data = topic,
             palette = c("#2E9FDF", "#E7B800","#fc0439","purple", "green", "black"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           #1) extracting cluster assignments
cluster_assignments <- res.km$cluster

clustered_data <- split(topic, cluster_assignments)
cluster1 <- clustered_data[[1]]
cluster2 <- clustered_data[[2]]
cluster3 <- clustered_data[[3]]
cluster4 <- clustered_data[[4]]
cluster5 <- clustered_data[[5]]
cluster6 <- clustered_data[[6]]

par(mfrow = c(3,2))

###cluster 1###
cluster1_averages <- colMeans(cluster1)
cluster1_averages <- cluster1_averages[order(-cluster1_averages)]
cluster1_averages_df <- data.frame(Variable = names(cluster1_averages),
                                   Average = cluster1_averages)

top_30_words_cl1 <- head(cluster1_averages_df[order(-
                                                  cluster1_averages_df$Average), ],30)

top_30_words_cl1 <- top_30_words_cl1[-c(1,2),]
#cl1 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5)); title("Cluster 1: Pipelin/Explor", line = -3, adj=0.5)
cl1 <- wordcloud(words = top_30_words_cl1$Variable, freq = top_30_words_cl1$Average, scale = c(6,1)); title("Cluster 1: Pipelin/Explor", line = 0, adj=0.5, cex.main = 2)
###cluster 2 ###

cluster2_averages <- colMeans(cluster2)
cluster2_averages <- cluster2_averages[order(-cluster2_averages)]
cluster2_averages_df <- data.frame(Variable = names(cluster2_averages),
                                   Average = cluster2_averages)

top_30_words_cl2 <- head(cluster2_averages_df[order(-
                                                  cluster2_averages_df$Average), ],30)

top_30_words_cl2 <- top_30_words_cl2[-c(1,2),]
#cl2 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5)); title("Cluster 2: Gas/Energi", line = -3, adj = 0.5)
cl2 <- wordcloud(words = top_30_words_cl2$Variable, freq = top_30_words_cl2$Average, scale = c(6,1)); title("Cluster 2: Gas/Energi", line = 0, adj = 0.5, cex.main = 2)
## cluster 3 ###
cluster3_averages <- colMeans(cluster3)
cluster3_averages <- cluster3_averages[order(-cluster3_averages)]
cluster3_averages_df <- data.frame(Variable = names(cluster3_averages),
                                   Average = cluster3_averages)

top_30_words_cl3 <- head(cluster3_averages_df[order(-
                                                  cluster3_averages_df$Average), ],30)

top_30_words_cl3 <- top_30_words_cl3[-c(1,2),]

#custom_colors <- c("#ffc300","#c70039","#900c3f","#581845","#2180D3")

#cl3 <- wordcloud(words = top_30_words$Variable, 
#                 color = custom_colors,
#                 freq = top_30_words$Average, scale = c(3,0.5)); title("Cluster 3: Oil/Barrel", line = -3, adj = 0.5)

cl3 <- wordcloud(words = top_30_words_cl3$Variable, freq = top_30_words_cl3$Average, scale = c(6,1)); title("Cluster 3: Oil/Barrel", line = 0, adj = 0.5, cex.main = 2)
####cluster4###

cluster4_averages <- colMeans(cluster4)
cluster4_averages <- cluster4_averages[order(-cluster4_averages)]
cluster4_averages_df <- data.frame(Variable = names(cluster4_averages),
                                   Average = cluster4_averages)

top_30_words_cl4 <- head(cluster4_averages_df[order(-
                                                  cluster4_averages_df$Average), ],30)
top_30_words_cl4 <- top_30_words_cl4[-c(1,2),]
#cl4 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5)); title("Cluster 4: WTI/Barrel", line = - 3, adj = 0.5)
cl4 <- wordcloud(words = top_30_words_cl4$Variable, freq = top_30_words_cl4$Average, scale = c(6,1)); title("Cluster 4: WTI/Barrel", line = 0, adj = 0.5, cex.main = 2)

####cluster5###

cluster5_averages <- colMeans(cluster5)
cluster5_averages <- cluster5_averages[order(-cluster5_averages)]
cluster5_averages_df <- data.frame(Variable = names(cluster5_averages),
                                   Average = cluster5_averages)

top_30_words_cl5 <- head(cluster5_averages_df[order(-
                                                  cluster5_averages_df$Average), ],30)
top_30_words_cl5 <- top_30_words_cl5[-c(1,2),]
#cl5 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(3,0.5)); title("Cluster 5: Energi/Gas", line = - 3, adj = 0.5)
cl5 <- wordcloud(words = top_30_words_cl5$Variable, freq = top_30_words_cl5$Average, scale = c(6,1)); title("Cluster 5: Energi/Gas", line = 0, adj = 0.5, cex.main = 2)

####cluster6###

cluster6_averages <- colMeans(cluster6)
cluster6_averages <- cluster6_averages[order(-cluster6_averages)]
cluster6_averages_df <- data.frame(Variable = names(cluster6_averages),
                                   Average = cluster6_averages)

top_30_words_cl6 <- head(cluster6_averages_df[order(-
                                                  cluster6_averages_df$Average), ],30)
top_30_words_cl6 <- top_30_words_cl6[-c(1,2),]
cl6 <- wordcloud(words = top_30_words_cl6$Variable, freq = top_30_words_cl6$Average, scale = c(6,1)); title("Cluster 6: Oil/Fuel", line = 0, adj = 0.5, cex.main = 2)


par(mfrow = c(3,2))
cl1 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 1: Pipelin/Explor", line = 0, adj=0.5)
cl2 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 2: Gas/Energi", line = 0, adj=0.5)
cl3 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 3: Oil/Barrel", line = 0, adj=0.5)
cl4 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 4: WTI/Barrel", line = 0, adj=0.5)
cl5 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 5: Energi/Gas", line = 0, adj=0.5)
cl6 <- wordcloud(words = top_30_words$Variable, freq = top_30_words$Average, scale = c(9,1)); title("Cluster 6: Oil/Fuel", line = 0, adj=0.5)

##### heat map ####
setwd("C:/Users/j1kxb09/OneDrive - FR Banks/Documents/MNSC meetings/data")
topic <- read.csv("topicVecs (3).csv", stringsAsFactors = FALSE)
topic$X <- year(topic$X) 
set.seed(12345)
res.km <- kmeans(scale(topic[, -which(names(topic) == "X")]), 6, nstart = 25)
res.km$cluster
cluster_assignments <- res.km$cluster

clustered_data <- split(topic, cluster_assignments)
cluster1 <- clustered_data[[1]]
cluster2 <- clustered_data[[2]]
cluster3 <- clustered_data[[3]]
cluster4 <- clustered_data[[4]]
cluster5 <- clustered_data[[5]]
cluster6 <- clustered_data[[6]]

bind <- rbind(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6)
bind_sum <- bind %>%
  group_by(X) %>%
  summarise(row_count = n())

sum_cluster1 <- cluster1 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 1)


sum_cluster2 <- cluster2 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 2)


sum_cluster3 <- cluster3 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 3)


sum_cluster4 <- cluster4 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 4)


sum_cluster5 <- cluster5 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 5)

sum_cluster6 <- cluster6 %>%
  group_by(X) %>%
  summarise(row_count = n(), cluster = 6)

cluster_topic <- rbind(sum_cluster1, sum_cluster2, sum_cluster3, sum_cluster4, sum_cluster5, 
                       sum_cluster6)

cluster_topic_year <- cluster_topic %>%
  group_by(cluster) %>%
  summarise(sum = sum(row_count))

matrix_data <- cluster_topic %>%
  pivot_wider(names_from = X, values_from = row_count, values_fill = 0)

p <- ggplot(cluster_topic, aes(x = X, y = as.factor(cluster), fill = row_count)) +
  geom_tile() + ggtitle("Figure 1: Topic Visualization by cluster and year") + 
  scale_fill_gradient(low = "lightblue", high = "blue", breaks = c(10, 20, 30, 40, 50, 60), labels = c(10, 20, 30, 40, 50, 60)) + 
  labs(x = "Year", y = "Cluster", fill = "Number of Topic Allocations For a Given Cluster") + scale_y_discrete(breaks = unique(cluster_topic$cluster)) +
  scale_x_continuous(breaks = seq(min(cluster_topic$X), max(cluster_topic$X), by = 2)) + geom_text(aes(x = 2024, y = 1, label = "158", size = 4)) + 
                                                                                         geom_text(aes(x = 2024, y = 2, label = "111", size = 4)) +
                                                                                         geom_text(aes(x = 2024, y = 3, label = "173", size = 4)) +
                                                                                         geom_text(aes(x = 2024, y = 4, label = "583", size = 4)) +
                                                                                         geom_text(aes(x = 2024, y = 5, label = "157", size = 4)) +
                                                                                         geom_text(aes(x = 2024, y = 6, label = "397", size = 4)) + theme_minimal() + scale_size(guide = "none")
p 

bind_sum <- as.data.frame(bind_sum)
#bind_sum$constant_y <- 1
bind_sum <- ggplot(bind_sum, aes(x = X, y = as.factor(row_count), fill = row_count)) + geom_tile() +
  scale_fill_gradient(low = "#8cd3f5", high = "navy", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  labs(x = "Year", fill = "Number of Topic Allocations in a Given Year") + ggtitle("Figure 2: Topic Visualization by year") +
  scale_x_continuous(breaks = seq(min(cluster_topic$X), max(cluster_topic$X), by = 2)) +
  theme_minimal() + theme(axis.title.y = element_blank())
combined_plot <- grid.arrange(p, bind_sum, ncol = 1)
