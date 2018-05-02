
fm = list.files(".",pattern = "density")
fm = fm[grep("*.txt",fm)]
nrow = ceiling(length(fm)/2)

if(length(fm)==1) pdf("DensityPlots_EBprotV2.pdf",height=6,width =6,useDingbats = F)
if(nrow>1) {
  pdf("DensityPlots_EBprotV2.pdf",height=12,width =12)
  par(mfrow=c(2,2))
}
if(nrow==1 & length(fm)>1) {
  pdf("DensityPlots_EBprotV2.pdf",height=6,width =12)
  par(mfrow=c(1,2))
}
for(i in 1:length(fm)){
  denfile = read.delim(fm[i],as.is=T)
  
  label = strsplit(fm[i],".txt")[[1]][1]
  label = strsplit(label,"densityfile_")[[1]][2]
  
  plot(denfile$Ratio, denfile$Density,cex=0.4,xlab="Ratio",ylab="Density",type="l", main=label,lwd=2)
  points(denfile$Ratio, denfile$Null_density,cex=0.3,col="darkgrey",type="l",lwd=1.5)
  points(denfile$Ratio, denfile$Pos_density,cex=0.3,col=2,type="l",lwd=1.5)
  points(denfile$Ratio, denfile$Neg_density,cex=0.3,col=4,type="l",lwd=1.5) 
  
  legend("topleft",c("overall density","null-density","up-density","down-density"),col=c("black","darkgrey","red","blue"),
         lty=1,cex=0.7,lwd=2)
}
dev.off()


results = read.delim("EBprot_results.txt",as.is=T)
nums = (length(names(results))-1)/6

if(nums==1) pdf("EBprotV2_PPscoreplot.pdf",height=6,width=6,useDingbats = F)
if(nums>2) {
  pdf("EBprotV2_PPscoreplot.pdf",height=12,width=12,useDingbats = F)
  par(mfrow=c(2,2))
}
if(nums==2) {
  pdf("EBprotV2_PPscoreplot.pdf",height=6,width=12,useDingbats = F)
  par(mfrow=c(1,2))
}

for(i in 1:nums){
  tmp_l = strsplit(names(results)[((i-1)*6)+2],"_")[[1]][2]
  index_r = ((i-1)*6)+2
  index_p = ((i-1)*6)+5
  index_b = (i*6)+1
  upperL = max(results[,index_r], na.rm=T)+0.05
  lowerL = min(results[,index_r],na.rm=T)-0.05
  plot(results[,index_r],results[,index_p],cex=0.5,main =tmp_l ,pch=19,xlab ="MedianLog2ratio",ylab="PPscore", xlim=c(lowerL,upperL), ylim = c(-1,1))
  cidup = which(results[,index_b]<0.05 & results[,index_r]>0)
  ciddown = which(results[,index_b]<0.05 & results[,index_r]<0)
  cutoffup = min(abs(results[cidup,index_p]))
  cutoffdown =min(abs(results[ciddown,index_p]))
  points(results[c(cidup,ciddown),index_r],results[c(cidup,ciddown),index_p],col=2,pch=19,cex=0.5)
  abline(h=cutoffup,lty=2,col=3)
  abline(h=-cutoffdown,lty=2,col=3)  
  legend("topleft",c("sig at 5% FDR"),cex=0.8,col="red",pch=19,bty="n")
}
dev.off()
