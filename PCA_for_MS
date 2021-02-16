#Plotting:
#RBM15
#FTO
#SND1
#HNRPC
#YTHD2
#VennDiagram
#install.packages("VennDiagram")
library(VennDiagram)
#install.packages("yarr")
library(yarrr)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

#getwd()
##Set as working directory in order to read files
setwd("~/_KSHV_Project/RAP-MS")
#Get data from original files MS and FA
MS_A <- read.csv("MS_A.csv", TRUE, ",")
MS_B <- read.csv("MS_B.csv", TRUE, ",")
MS_C <- read.csv("MS_C.csv", TRUE, ",")

FA_A <- read.csv("FA_A.csv", TRUE, ",")
FA_B <- read.csv("FA_B.csv", TRUE, ",")
FA_C <- read.csv("FA_C.csv", TRUE, ",")

Scores_all <- read.table("ScoresMatrix.tex",
                         header = TRUE, na.strings = "NA")

class(Scores_all)
summary(Scores_all)
plot(Scores_all)

Scores.matrix <-data.matrix(Scores_all)
pca <- prcomp(t(Scores.matrix[,2:7]), scale=TRUE)
## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])
## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

#######
##
## NOTE: Everything that follow is just bonus stuff.
## It simply demonstrates how to get the same
## results using "svd()" (Singular Value Decomposition) or using "eigen()"
## (Eigen Decomposition).
##
#######

############################################
##
## Now let's do the same thing with svd()
##
## svd() returns three things
## v = the "rotation" that prcomp() returns, this is a matrix of eigenvectors
##     in other words, a matrix of loading scores
## u = this is similar to the "x" that prcomp() returns. In other words,
##     sum(the rotation * the original data), but compressed to the unit vector
##     You can spread it out by multiplying by "d"
## d = this is similar to the "sdev" value that prcomp() returns (and thus
##     related to the eigen values), but not
##     scaled by sample size in an unbiased way (ie. 1/(n-1)).
##     For prcomp(), sdev = sqrt(var) = sqrt(ss(fit)/(n-1))
##     For svd(), d = sqrt(ss(fit))
##
############################################

svd.stuff <- svd(scale(t(data.matrix), center=TRUE))

## calculate the PCs
svd.data <- data.frame(Sample=colnames(data.matrix),
                       X=(svd.stuff$u[,1] * svd.stuff$d[1]),
                       Y=(svd.stuff$u[,2] * svd.stuff$d[2]))
svd.data

## alternatively, we could compute the PCs with the eigen vectors and the
## original data
svd.pcs <- t(t(svd.stuff$v) %*% t(scale(t(data.matrix), center=TRUE)))
svd.pcs[,1:2] ## the first to principal components

svd.df <- ncol(data.matrix) - 1
svd.var <- svd.stuff$d^2 / svd.df
svd.var.per <- round(svd.var/sum(svd.var)*100, 1)

ggplot(data=svd.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", svd.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", svd.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("svd(scale(t(data.matrix), center=TRUE)")

############################################
##
## Now let's do the same thing with eigen()
##
## eigen() returns two things...
## vectors = eigen vectors (vectors of loading scores)
##           NOTE: pcs = sum(loading scores * values for sample)
## values = eigen values
##
############################################
cov.mat <- cov(scale(t(data.matrix), center=TRUE))
dim(cov.mat)

## since the covariance matrix is symmetric, we can tell eigen() to just
## work on the lower triangle with "symmetric=TRUE"
eigen.stuff <- eigen(cov.mat, symmetric=TRUE)
dim(eigen.stuff$vectors)
head(eigen.stuff$vectors[,1:2])

eigen.pcs <- t(t(eigen.stuff$vectors) %*% t(scale(t(data.matrix), center=TRUE)))
eigen.pcs[,1:2]

eigen.data <- data.frame(Sample=rownames(eigen.pcs),
                         X=(-1 * eigen.pcs[,1]), ## eigen() flips the X-axis in this case, so we flip it back
                         Y=eigen.pcs[,2]) ## X axis will be PC1, Y axis will be PC2
eigen.data

eigen.var.per <- round(eigen.stuff$values/sum(eigen.stuff$values)*100, 1)

ggplot(data=eigen.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", eigen.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", eigen.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("eigen on cov(t(data.matrix))")

class(MS_A)
barplot(MS_A$Score)
head(MS_A)
MS_A$
plot(MS_A$ProteinID, MS_A$Score)

barplot(MS_B$Score)

Norm_Score <- log(MS_A$Score,2)
MS_A.new <- cbind(MS_A,Norm_Score)
head(MS_A.new)

plot(MS_A.new$Coverage)

ID.A <- (MS_A$ProteinID)
ID.B <- (MS_B$ProteinID)
ID.C <- (MS_A$ProteinID)

#ID.A!=ID.B
ID.AB <- Reduce(intersect, list(MS_A$ProteinID,MS_B$ProteinID))
head(ID.AB)

#Show summary output
tail(MS_A)
tail(MS_B)
tail(MS_C)
#Shwo the header of all the MS samples
head(MS_A)
head(MS_B)
head(MS_C)

#To create a matrix
#Order raws
#MS_A.sorted <- sort(MS_A$ProteinID)
summary(MS_A)
#cor(MS_A[1:200,2], MS_B[1:200,2])
cor(MS_A[1:2000,2], MS_B[1:2000,2])

plot(MS_A[1:2384,5], FA_A[1:2384,5],
     col=2,
     pch=16,
     cex=1,
     #xlim = c(0,15000),
     #ylim = c(0,10000),
     xlab = "X",
     ylab = "S",
     main = "Correlations")
  points(MS_B[1:2384,5], FA_B[1:2384,5],
  pch = 16,
  cex= 1,
  col=3
  #col = "pink")
  #col = transparent("pink", trans.val = .8))
  #col = brewer.pal(8, "Pastel2"))
  )
  points(MS_C[1:2384,5], FA_C[1:2384,5],
         pch = 16,
         cex=1,
         #col = brewer.pal(9, "Pastel1")
         col=4)

#Venndiagram for MS
venn.diagram(
x= list(MS_A$ProteinID,MS_B$ProteinID, MS_C$ProteinID),
category.names = c("Domain I","Domain II","Domain III"),
filename = "test_vennMS.png",
output=TRUE
)

#Venndiagram for FA
venn.diagram(
  x= list(FA_A$ProteinID,FA_B$ProteinID, FA_C$ProteinID),
  category.names = c("Domain I","Domain II","Domain III"),
  filename = "test_vennFA.png",
  output=TRUE
)

#Venndiagram all
venn.diagram(
  x= list(FA_A$ProteinID,MS_A$ProteinID),
  category.names = c("Domain I-FA","Domain I-MS"),
  filename = "test_vennI.png",
  output=TRUE
)

#Venndiagram all
venn.diagram(
  x= list(FA_B$ProteinID,MS_B$ProteinID),
  category.names = c("Domain II-FA","Domain II-MS"),
  filename = "test_vennII.png",
  output=TRUE
)

#Venndiagram all
venn.diagram(
  x= list(FA_C$ProteinID,MS_C$ProteinID),
  category.names = c("Domain III-FA","Domain III-MS"),
  filename = "test_vennIII.png",
  output=TRUE
)

#Venndiagram for MS
venn.diagram(
  x= list(MS_A$ProteinID,MS_B$ProteinID, FA_A$ProteinID,FA_B$ProteinID),
  category.names = c("Domain I-MS","Domain II-MS","Domain I-FA","Domain II-FA"),
  filename = "test_vennAB.png",
  output=TRUE
)

venn.diagram(
  x= list(MS_C$ProteinID,MS_B$ProteinID, FA_C$ProteinID,FA_B$ProteinID),
  category.names = c("Domain III-MS","Domain II-MS","Domain III-FA","Domain II-FA"),
  filename = "test_vennBC.png",
  output=TRUE
)


pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")


plot(MS_A[730:2215,5], FA_A[730:2215,5],
     col=2,
     pch=16,
     cex=1,
     xlim = c(0,100),
     ylim = c(0,100),
     xlab = "MS",
     ylab = "FA",
     main = "Correlations")
points(MS_B[730:2215,5], FA_B[730:2215,5],
       pch = 16,
       cex= 1,
       col=3
       #col = "pink")
       #col = transparent("pink", trans.val = .8))
       #col = brewer.pal(8, "Pastel2"))
)
points(MS_C[730:2215,5], FA_C[730:2215,5],
       pch = 16,
       cex=1,
       #col = brewer.pal(9, "Pastel1")
       col=4)
head(Scores_all)
head(Scores.matrix)
barplot(Scores.matrix[731:991,2:7])
Scores.matrix[731:991,2:7]

Scores.norm2 <- log2(Scores.matrix[,2:7])
rownames(Scores.matrix) <- Scores_all$ProteinID

#Plot
plot(Scores.norm2,
     pch=16,
     cex=1,
     col="grey",
     #xlim = c(4,10),
     #ylim = c(4,10),
     xlab = "MS",
     ylab = "FA",
     main = "PCA1"
     #xaxt = "n",
     #yaxt = "n"
     )
points(Scores.norm2[77:79,],
       pch = 16,
       cex=2,
       col="purple")
points(Scores.norm2[11:13,],
       pch = 16,
       cex=2,
       #col = brewer.pal(9, "Pastel1")
       col="purple")
text(9,10, "SND1")
text(9.9,9.5, "HNRPC")
text(10.9,10.3, "YTHD2")
text(7.6,8 ,"RBM15")
text(8.5,7.9,"FTO")
