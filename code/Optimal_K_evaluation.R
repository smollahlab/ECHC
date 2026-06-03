# Script to evaluate optimal k cluster determination
 
#Load input data
histon_data<-read.table("data/preprocessed/histone.csv", sep=",", header=T, row.names=1)
# Define a function to calculate total within-cluster sum of squares
wss <- function(k) {
  kmeans(histone_data, k, nstart = 10)$tot.withinss
}
 
# Apply this function for k = 1 to 10
k_values <- 1:10
wss_values <- sapply(k_values, wss)
 
# Plot the results
plot(k_values, wss_values, type="b", pch = 19, frame = FALSE, xlab="Number of clusters K", ylab="Total within-clusters sum of squares",
main = "Elbow Method")
