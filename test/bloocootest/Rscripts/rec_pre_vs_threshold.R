#!/usr/bin/Rscript

args=commandArgs(TRUE)

tab_filename = args[1]
output_filename = args[2]

## h=T si nom des colonnes = première ligne (sinon h=F)
tab = read.table(tab_filename, h=T)

bitmap(output_filename,type="png256",width=12,height=7,res=400)
plot(range(tab$cov_t),c(0,100),type="n",main="recall/precision vs bloom threshold (fixed error rate = 1%)",xlab="bloom threshold",ylab="recall and precision (%)")
lines(tab$cov_t,tab$rec,type="b",pch=1,col=1)
lines(tab$cov_t,tab$pre,col=2,pch=3,type="b")
legend("bottomright",legend=c("recall","precision"),col=c(1,2),pch=c(1,3))
d=dev.off() ## pour fermer le device png
