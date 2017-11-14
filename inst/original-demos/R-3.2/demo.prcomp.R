# Small script to illustrate the 'prcomp' function for performing a
# principal component analysis of a small data set.
library(ggplot2)

# Load the Iris data.
data(iris)

# Perform PCA using SVD.
out <- prcomp(iris[1:4])

# Generate text summary of PCA results.
print(out)

# Generate another text summary of PCA results.
print(summary(out))

# Generate visual summary of the PCA results.
print(ggplot(cbind(iris,out$x),aes(x = PC1,y = PC2,color = Species)) +
      geom_point() +
      theme_minimal())
