{\rtf1\ansi\ansicpg1252\cocoartf2868
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red10\green10\blue10;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c3922\c3922\c3922;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs37\fsmilli18667 \cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 # Script to evaluate optimal k cluster determination
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 \'a0
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 #Load input data
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 histon_data<-read.table("histone.csv", sep=",", header=T, row.names=1)
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 # Define a function to calculate total within-cluster sum of squares
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 wss <- function(k) \{
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 \'a0\'a0kmeans(histone_data, k, nstart =\'a010)$tot.withinss
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 \}
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 \'a0
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 # Apply this function for k = 1 to 10
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 k_values <- 1:10
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 wss_values <- sapply(k_values, wss)
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 \'a0
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 # Plot the results
\fs26\fsmilli13333 \cf0 \strokec3 \

\fs37\fsmilli18667 \cf2 \strokec2 plot(k_values, wss_values, type="b", pch = 19, frame = FALSE, xlab="Number of clusters K",\'a0ylab="Total within-clusters sum of squares",
\fs26\fsmilli13333 \cf0 \strokec3 \
\pard\pardeftab720\qj\partightenfactor0

\fs37\fsmilli18667 \cf2 \strokec2 main = "Elbow Method")
\fs28 \cf0 \strokec3 \
}