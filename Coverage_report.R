#!/usr/bin/env Rscript

# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =======
#   Attributions
# This code is partially copied from https://gist.github.com/stephenturner/9396409#file-exome-coverage-multi-r 
# by Stephen Turner


bamname=commandArgs(trailingOnly = TRUE)[1]
regionsname=commandArgs(trailingOnly = TRUE)[2]
outprefix=commandArgs(trailingOnly = TRUE)[3]
options(echo=TRUE)
print(bamname)

cov=read.table(text=system(paste("bedtools coverage -abam ",bamname," -b ",regionsname," -hist | grep ^all",sep=""),intern=T))

cov_cumul <- 1-cumsum(cov[,5])

library(RColorBrewer)
cols <- brewer.pal(3, "Dark2")
col=cols[1]

# Save the graph to a file
png(paste(outprefix,"_exome-coverage-ecdf.png",sep=""), h=600, w=900, pointsize=20)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[2:nrow(cov), 2], cov_cumul[1:(length(cov_cumul)-1)], type='l', ,lwd=3,col=col,xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage",xlim=c(0,60000))

dev.off()

png(paste(outprefix,"_exome-coverage-hist.png",sep=""),h=600,w=900,pointsize=20)
seqs=seq(-1,60000,1000)
indices=.bincode(cov[,2],seqs)
covhist=rep(0,length(seqs))
names(covhist)=seqs
covhist[unique(indices)]=tapply(cov[,5],indices,sum)
b=barplot(covhist,axes = F,ylab = "Percent of bases",names.arg = "",xlab="Bases read depth",main="Depth by base")
axis(1,at =b[seq(1,length(b),5)],labels=(seqs+1)[seq(1,length(b),5)])
axis(2)
dev.off()
