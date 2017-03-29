######################
# PCA of 12A and ETI #
######################

# PCA Analysis in R:
# Added 0.0001 to all values before running PCA
# Read in data and find principle components
# Make sure not to have spaces in the column titles and to leave the top left cell blank
x <- AllFPKMs2[c("Hemi.12A","ETI.120hr","ETI.240hr","ETI.360hr","ETI.576hr")]
colnames(x)<- c("WT 120hr","ETI 120hr","ETI 240hr", "ETI 360hr", "ETI 576hr")
x<-x[-(which(rownames(AllFPKMs2) %in% c("noe"))),]
xx <- prcomp(t(x))
# Visualize the results in 2D using biplot
biplot(xx, var.axes=FALSE, ylabs=NULL)
# Visualize the results in 3D using rgl library - may need to install rgl first
# May want to change the radius in the spheres3d command if too big or small
library(rgl)
plot3d(xx$x,xlab="PC1",ylab="PC2",zlab="PC3",type="s", col=rainbow(length(xx$x[,1])))
plot3d(xx$x,xlab="PC1",ylab="PC2",zlab="PC3",type="h", add = TRUE)
grid3d(side="z", at=list(z=0))
text3d(xx$x, text=rownames(xx$x), adj=1.5, cex=0.8)

sink("PCA ETI and 12A minus noe - percentages.txt")
summary(xx)
sink()

xx$rotation
write.table(xx$rotation, file="PCA ETI and 12A minus noe - genes contributions.txt", sep="\t")