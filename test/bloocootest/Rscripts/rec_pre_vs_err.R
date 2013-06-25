#!/usr/bin/Rscript

args=commandArgs(TRUE)

tab_filename = args[1]
output_filename = args[2]

## h=T si nom des colonnes = première ligne (sinon h=F)
tab = read.table(tab_filename, h=T)

bitmap(output_filename,type="png256",width=12,height=7,res=400)
plot(range(tab$err),c(0,100),type="n",main="recall/precision vs error rate (fixed bloom threshold = 4)",xlab="error rate",ylab="recall and precision (%)")
lines(tab$err,tab$rec,type="b",pch=1,col=1)
lines(tab$err,tab$pre,col=2,pch=3,type="b")
legend("bottomright",legend=c("recall","precision"),col=c(1,2),pch=c(1,3))
d=dev.off() ## pour fermer le device png
