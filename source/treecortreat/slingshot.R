BiocManager::install("slingshot")
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(TSCAN)
library(tidyverse)
means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
sce <- SingleCellExperiment(assays = List(counts = counts))

## bifurcating
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

crv1 <- getCurves(lin1)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
str(crv1)

dist_ind_curve1 = crv1@curves$curve1$dist_ind
dist_ind_curve2 = crv1@curves$curve2$dist_ind
dist_ind_curve = sapply(1:length(dist_ind_curve1), function(i) min(dist_ind_curve1[i], 
                                                                   dist_ind_curve2[i]))

dist_curve1 = crv1@curves$curve1$dist
orth_proj = crv1@curves$curve1$s

m = colMeans(rd)
tot_var = sapply(1:nrow(rd), function(r){
  sum((rd[r,]-m)^2)
}) %>% sum

tot_var = sum((sweep(rd, 2, colMeans(rd), '-'))^2) 
R2 = 1- dist_curve1/tot_var

sapply(1:ncol(rd), function(co){
  var(rd[,co]) * (nrow(rd)-1)
}) %>% sum
