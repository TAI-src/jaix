library(philentropy)
library(RColorBrewer)
library(fields)

plot_mat = function(mat, nbins=256, title=""){
  values = as.vector(mat)
  t_min = min(values[is.finite(values)], na.rm=TRUE)
  t_max = max(values[is.finite(values)], na.rm=TRUE)
  z_cap = max(abs(t_min), abs(t_max))

  palette <- colorRampPalette(c("red",'white','#0033BB'))(nbins)
  image.plot(t(mat[rev(1:nrow(mat)),]),col = palette, axes=FALSE, 
             lab.breaks=NULL, main=title, zlim=c(-z_cap,z_cap))
  
  
  # axis labels
  image(t(mat[rev(1:nrow(mat)),]), col = palette, axes=FALSE, add=TRUE,
        zlim=c(-z_cap,z_cap))
  axis(3, at=seq(0,1, length=ncol(mat)), labels=colnames(mat), lwd=0, pos=1, las=2)
  axis(2, at=seq(1,0, length=nrow(mat)), labels=rownames(mat), lwd=0, pos=-0.05, las=1)
  
  # add values
  e <- expand.grid(seq(0,1, length=ncol(mat)), seq(1,0, length=nrow(mat)))
  text(e, labels=t(format(mat, digits=2)))
  
}


process = function(data, field){
  results = matrix(NA, nrow=length(unique(data$id)), ncol=6)
  counter = 1
  for(fid in sort(unique(data$id))){
    vals = c()
    tmp = data[data$id==fid,]
    for(fold in unique(tmp$to)){
      self = tmp[tmp$from==fold & tmp$to==fold, field]
      others = tmp[tmp$from!=fold & tmp$to==fold, field]
      vals = c(vals, self-others)
    }
    results[counter,] = vals
    counter = counter +1
  }
  colnames(results) = paste("i", c("1.1", "1.2", "2.1", "2.2", "3.1", "3.2"), sep="_")
  rownames(results) = paste("d", 0:(nrow(results)-1), sep="_")  

  plot_mat(results, title = field)
}


data = read.csv("fvq1pkfe_results.csv")
data$set_matches[is.na(data$set_matches)]=0
pdf("fvq1pkfe_results.pdf")
for(field in c("bin_dist.mean", "bin_dist.med", "bin_dist.min", "bin_dist.max", "set_matches", "number")){
#for(field in c("bin_dist.max")){
  process(data,field)
}
dev.off()

