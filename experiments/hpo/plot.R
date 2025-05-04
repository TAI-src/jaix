library(philentropy)
library(RColorBrewer)
library(fields)
library(paletteer)
require(data.table)

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

# total training duration per id
data = read.csv("fvq1pkfe_env_steps.csv")
tmp = read.csv("fvq1pkfe_data.csv")
tmp = aggregate(tmp,by=list(tmp$f, tmp$fold), min)
data = merge(data,tmp,by.x=c("id","fold"), by.y=c("f", "fold"))

nColor=20
colors = paletteer_c("ggthemes::Red-Blue Diverging", n=nColor)
min_r = as.numeric(cut(data$r, nColor)) 

par(mar=c(5.1, 4.1, 4.1, 3.1), xpd=TRUE)
plot(log(data$env_step), data$id, col=colors[min_r], pch=19,
     main="Total Sweep Duration", xlab="Log of Total Training Time [s]",
     ylab="Dataset id")
for(d in unique(data$id)){
  lines(x=c(min(log(data$env_step)), max(log(data$env_step))), y=c(d,d), col="lightgrey")
}
lab = levels(cut(data$r, nColor))[c(seq(1,nColor,5),20)]
legend("topright", legend=lab, pch=19, col=colors[c(seq(1,nColor,5),20)],
       title="Best ensemble rank")

dev.off()

dec2bin <- function(x, n=9){
  res = paste(as.integer(intToBits(x)), collapse = "")
  res = as.integer(strsplit(substr(res, 1, n), "")[[1]])
  return(res)
}

dec2which = function(x, n=9){
  res = dec2bin(x,n)
  res = which(as.logical(res))
  return(res)
}

dec2which_str = function(x, n=9){
  return(paste(dec2which(x, n), collapse=";"))
}


bp = function(fmla, data, f){
  if(fmla[[2]]=="r"){
    ylab="Ensemble rank"
  }else if(fmla[[2]]=="t"){
    ylab="Training time [s]"
  }
  boxplot(fmla, data=data, main=paste("Dataset", f),
          col=cols, xlab="", ylab=ylab,
          xaxt="n")
  axis(side=1, at=c(4.5, 13.5, 22.5), labels=paste("fold", 0:2))
  legend("bottomright", inset=c(-0.2,0), legend=1:9,
         col=cols, pch=19, title=fmla[[3]][[2]])
}

mop = function(col_col, shape_col, data, f, highlight=NULL){
  shape_fac = factor(data[[shape_col]])
  if(length(levels(shape_fac))==3){
    scn = 2
    shps = shapes
  }else{
    scn = length(levels(shape_fac))-1
    shps = 0:scn
  }
  if (!is.null(highlight)){
    shps = rep(20, length(levels(shape_fac)))
    shps[levels(shape_fac)==highlight] = 8
  }
  tmp_shapes = shps[as.numeric(shape_fac)]
  tmp_cols = cols[as.numeric(factor(data[[col_col]]))]
  plot(data$r, data$t, pch=tmp_shapes, col=tmp_cols,
       xlab="Ensemble rank", ylab="Training time [s]",
       main=paste("Dataset", f))
  legend("topright", inset=c(-0.2,0), legend=0:scn,
         pch=shps, title=shape_col)
  legend("bottomright", inset=c(-0.2,0), legend=1:9,
         col=cols, pch=19, title=col_col)
}

vs_plot = function(tmp){
  plot(NA, xlim=c(0,max(tmp$r)), ylim=c(0, max(tmp$r)),
       xlab="Rank fold x", ylab="Rank fold y",
       main="Same ensemble - different folds")
  combs = c()
  for(f1 in 0:1){
    tmpa = tmp[tmp$fold==f1,]
    tmpa = tmpa[order(tmpa$x),]
    for(f2 in (f1+1):2){
      combs = c(combs, paste("f", f1, "vs", f2))
      tmpb = tmp[tmp$fold==f2,]
      tmpb = tmpb[order(tmpb$x),]
      points(tmpa$r, tmpb$r, pch=4, col=cols_3[length(combs)])
    }
  }
  legend("topright", inset=c(-0.25,0), legend=combs,title="folds",
         pch=19, col=cols_3)
}

data = read.csv("fvq1pkfe_data.csv")
shapes = c(1,4,5)
# nbapalettes::jazz_city
cols = c("#010101FF", "#5D2A2CFF", "#93282CFF",
         "#C8102EFF", "#DA291CFF", "#DC582AFF",
         "#E87722FF", "#FF9E1BFF", "#FFC72CFF")
cols_3 = c("#5B1414FF", "#AD722CFF", "#1A6384FF")
for(f in sort(unique(data$f))){
  pdf(paste("fvq1pkfe_", f, ".pdf", sep=""))
  par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
  
  tmp = data[data$f==f,]
  tmp$num_m = sapply(tmp$x, function (x) sum(dec2bin(x)))
  
  vs_plot(tmp)
  mop("num_m", "fold", tmp, f)
  for(fold in sort(unique(tmp$fold))){
    mop("num_m", "fold", tmp, paste(f,"highlight fold", fold), highlight=fold)
  }
  bp(r ~ num_m* fold, data=tmp, f)
  bp(t ~ num_m* fold, data=tmp, f)
  
  # decode x which
  which_x = data.frame(w = sapply(1:max(tmp$x), dec2which_str), x = 1:max(tmp$x))
  dt = data.table(which_x)
  dt = dt[, list(w=unlist(strsplit(w, ";"))), by=x]
  dt$which_m = as.numeric(dt$w)
  
  full = merge(tmp, dt, by="x")
  
  
  mop("which_m", "fold", full, f)
  for(fold in sort(unique(full$fold))){
    mop("num_m", "which_m", full[full$fold==fold,], paste(f,"fold", fold))
    for(w in sort(unique(full$which_m))){
      mop("num_m", "which_m", full[full$fold==fold,],
          paste(f,"fold", fold, "highlight wich_m", w-1), highlight = w)
    } 
  }
  bp(r ~ which_m* fold, data=full, f)
  bp(t ~ which_m* fold, data=full, f)
  
  for(fold in sort(unique(tmp$fold))){
    tmpf = tmp[tmp$fold==fold,]
    hist(tmpf$r, xlab="Ensemble rank",
         main=paste("Fold", fold))
  }
  singles = full[full$num_m ==1,]
  mop("which_m", "fold", singles, f)
  
  dev.off()
}


plot_across = function(data){
  plot(data$r, data$f, col=data$cols, pch=data$shapes,
       main=paste("Singles", paste(sort(unique(data$which)), collapse=" ")),
       xlab="Ensemble rank",
       ylab="Dataset id")
  for(d in unique(data$f)){
    lines(x=c(0, max(data$r)), y=c(d,d), col="lightgrey")
  }
  legend("topright", inset=c(-0.2,0), legend=0:2,
         pch=shapes, title="folds")
  legend("bottomright", inset=c(-0.2,0), legend=1:9,
         col=cols, pch=19, title="which_m")  
}

pdf("fvq1pkfe_singles.pdf")
par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
data = read.csv("fvq1pkfe_data.csv")
data$num_m = sapply(data$x, function (x) sum(dec2bin(x)))
data= data[data$num_m == 1,]
data$which = sapply(data$x, dec2which)
data$cols = cols[as.numeric(factor(data$which))]
data$shapes = shapes[as.numeric(factor(data$fold))]
plot_across(data)
for(w in sort(unique(data$which))){
  plot_across(data[data$which==w,])
}


dev.off()
