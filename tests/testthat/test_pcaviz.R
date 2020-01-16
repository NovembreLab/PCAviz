context("PCAviz")

# ----------------------------------------------------------------------
test_that("POPRES data set loads correctly",{
  data(popres)
  expect_equal(sd(popres$x$PC1),0.027,tolerance = 0.001,scale = 1)
})

# ----------------------------------------------------------------------
test_that("RegMap data set loads correctly",{
  data(regmap)
  expect_equal(sd(regmap$PC1),39.85,tolerance = 0.01,scale = 1)
})

# ----------------------------------------------------------------------
test_that("pcaviz(out.pca) invocation works with princomp() output",{
  data(iris)
  out1 <- pcaviz(princomp(iris[1:4]))
  out2 <- pcaviz(princomp(iris[1:4]),dat = iris)
  expect_output(print(summary(out1)),"principal components")
  expect_output(print(summary(out2)),"principal components")
})

# ----------------------------------------------------------------------
test_that("Alternative invocations of pcaviz() give same result",{
  data(iris)
  out.pca <- prcomp(iris[1:4])
  iris1 <- pcaviz(out.pca,dat = iris)
  iris2 <- with(out.pca,
                pcaviz(x = x,sdev = sdev,var = sum(sdev^2),
                       rotation = rotation,dat = iris))
  iris3 <- with(out.pca,
                pcaviz(dat = cbind(x,iris),sdev = sdev,var = sum(sdev^2),
                       rotation = rotation))
  iris4 <- with(out.pca,
                pcaviz(dat = cbind(x,iris),sdev = sdev,var = sum(sdev^2),
                       rotation = rotation,pc.cols = paste0("PC",1:4)))
  expect_output(print(summary(iris1)),"principal components")
  expect_identical(iris1,iris2)
  expect_identical(iris1,iris3)
  expect_identical(iris1,iris4)
})

# ----------------------------------------------------------------------
test_that("pcaviz.subset() correctly selects rows of POPRES data",{
  data(popres)
  popres1 <- pcaviz(dat = popres$x)

  # Remove the (11) Russian and Scottish samples.
  popres2 <- subset(popres1,!(abbrv == "RU" | abbrv == "Sct"))

  # Check that 11 samples were removed.
  expect_output(print(summary(popres1)),"1387")
  expect_output(print(summary(popres2)),"1376")
})

# ----------------------------------------------------------------------
test_that("pcaviz_transform2d() works correctly",{
  X <- rbind(c(0,0),
             c(1,0),
             c(0,1),
             c(1,1))
  Y <- rbind(c(-1,2),
             c(-1,7),
             c(-3,2),
             c(-3,7))
  colnames(X) <- c("PC1","PC2")
  colnames(Y) <- c("PC1","PC2")
  x <- pcaviz(x = X)
  y <- pcaviz_transform2d(x,
                          angle     = -90,
                          reflect.x = TRUE,
                          reflect.y = TRUE,
                          scale     = c(2,5),
                          a         = c(-1,2))
  expect_equal(y$data,as.data.frame(Y),tolerance = 1e-8,scale = 1)
})

# ----------------------------------------------------------------------
test_that("pcaviz_rescale_to_minimize_whitespace() works correctly",{
  data(iris)
  iris <- pcaviz(prcomp(iris[1:4]),dat = iris)
  out  <- pcaviz_reduce_whitespace(iris)
  expect_equal(diff(range(out$data$PC1)),
               diff(range(out$data$PC2)),
               tolerance = 1e-8)
})

# ----------------------------------------------------------------------
test_that("pcaviz_abbreviate_var() provides correct country abbreviations",{
  data(popres)
  popres <- pcaviz(dat = popres$x)
  popres <- subset(popres,
                   !(abbrv == "Sct" | country == "Serbia and Montenegro"))
  y <- pcaviz_abbreviate_var(popres,"country")$data$country.abbrv
  expect_identical(as.character(popres$data$abbrv),as.character(y))
})

# ----------------------------------------------------------------------
test_that("screeplot() does not produce an error on Iris data",{
  data(iris)
  iris <- pcaviz(prcomp(iris[1:4]),dat = iris)
  expect_s3_class(screeplot(iris),"ggplot")
})

# ----------------------------------------------------------------------
test_that("screeplot() generates error when s.d.'s are not provided",{
  data(iris)
  out.pca <- prcomp(iris[1:4])
  iris    <- with(out.pca,pcaviz(x = x,dat = iris))
  expect_error(screeplot(iris))
})

# ----------------------------------------------------------------------
test_that("pcaviz_violin() does not produce an error on Iris data",{
  data(iris)
  iris <- pcaviz(prcomp(iris[1:4]),dat = iris)
  expect_s3_class(pcaviz_violin(iris),"ggplot")
})

# ----------------------------------------------------------------------
test_that(paste("plot.pcaviz() generates a variety of plots from Iris",
                "data without error"),{
  data(iris)
  iris <- cbind(iris,data.frame(id = 1:150))
  iris <- transform(iris,id = as.character(id))
  iris <- pcaviz(prcomp(iris[1:4]),dat = iris)

  # Show "Species" as different colors and shapes.
  expect_s3_class(plot(iris,draw.points = TRUE),"ggplot")

  # Show a continuous variable (Petal Width) in different colors.
  expect_s3_class(plot(iris,draw.points = TRUE,color = "Petal.Width",
                       shape = "Species"),"ggplot")
  
  # Draw the rotated PC axes.
  iris.rot <- pcaviz_reduce_whitespace(iris)
  iris.rot <- pcaviz_rotate(iris.rot,15)
  expect_s3_class(plot(iris.rot,draw.points = TRUE,preserve.scale = TRUE),
                  "ggplot")

  # Plot PC1 vs. Petal Width, with lines showing the linear fit and
  # confidence interval.
  expect_s3_class(plot(iris.rot,coords = c("PC1","Petal.Width"),
                       draw.points = TRUE),"ggplot")

  # Plot all combinations of PCs 1-4.
  expect_s3_class(plot(iris.rot,coords = paste0("PC",1:4),draw.points = TRUE,
                       group = NULL,scale.pc.axes = 0.7),"ggplot")

  # Label all data points by their sample id.
  expect_s3_class(plot(iris,label = "id"),"ggplot")

  # Label all data points by "Species" assignment, using an
  # abbreviated label.
  expect_output(plot(iris),"Species.abbrv")
})

# ----------------------------------------------------------------------
test_that(paste("plot.pcaviz(plotly = TRUE) generates interactive plot",
                "from Iris data without error"),{
  data(iris)
  iris <- cbind(iris,data.frame(id = 1:150))
  iris <- transform(iris,id = as.character(id))
  iris <- pcaviz(prcomp(iris[1:4]),dat = iris)
  expect_s3_class(plot(iris,plotly = TRUE,tooltip = "id"),"plotly")
})
