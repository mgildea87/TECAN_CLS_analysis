setwd("Z:/mag456/TECAN")
sample = "2018_10_2 Tecan data.csv"
data<-read.csv(sample, skip = 57)
data = as.matrix(data)
rownames(data) <- data[,1]
data <- data[,-1]
class(data) <- "numeric"
tp = ncol(data)


#consolodating mean OD values
samplet.list <- grep("([A-Z][1-9])",row.names(data))
con.data = matrix(ncol = ncol(data), nrow = length(samplet.list))
row.names(con.data) = row.names(data)[samplet.list]
for (i in 1:length(samplet.list)){
  con.data[i,] = data[samplet.list[i]+3,]
}

#consolidate rows
if (nrow(con.data) > 96){
    con.data2 = matrix(ncol = ncol(con.data)*2, nrow=nrow(con.data)/2)
    for (i in 1:96){
      con.data2[i,] = c(con.data[i,], con.data[i+96,])
    }
    colnames(con.data2) <- c(1:ncol(con.data2))
    row.names(con.data2) <- row.names(con.data)[1:96]
    nas = which(is.na(con.data2[1,]))[1]
    con.data2 = con.data2[,1:nas]
} else {
  con.data2 = con.data
}

# time points (seconds)
times = data[2,]
times2 = c(times+times[tp]+689.3)
times = c(times, times2)

window = 10
first.deriv <- matrix(ncol = ncol(con.data2), nrow = nrow(con.data2))
row.names(first.deriv) = row.names(con.data2)
for (j in 1:nrow(con.data2)){  
  for (i in 1:(ncol(con.data2)-window)){
    first.deriv[j,i] = summary(lm(con.data2[j,i:(i+window)] ~ times[i:(i+window)]))$coefficients[2,1]
  }
}

#tMAX
Tmax = vector(length = nrow(first.deriv), mode = "numeric")
for (i in 1:nrow(first.deriv)){
  Tmax[i] = times[which(first.deriv[i,] == max(first.deriv[i,], na.rm=TRUE))]
}

#Background window
t2 = vector(length = nrow(first.deriv), mode = "numeric")
for (i in 1:nrow(first.deriv)){
  t = which(times == Tmax[i])
  t2[i] = which(first.deriv[i,5:t] > 0.000003)[1]
} 


#doubling time
doubling.time = vector(length = nrow(first.deriv), mode = "numeric")
r.squared = vector(length = nrow(first.deriv), mode = "numeric")
Log2.first.deriv = log2(first.deriv)
for (j in 1:nrow(Log2.first.deriv)){
    if (is.na(t2[j])){
      next
    }
    back.window = t2[j]+5
    r = mean(con.data2[j,5:back.window])
    v = vector(length = ncol(Log2.first.deriv), mode = "numeric")
    for (i in 1:length(v)){
      v[i] = con.data2[j,i]-r
    }
    v = log2(v)
    l = which(times == Tmax[j])
    doubling.time[j] = summary(lm(v[(l-25):(l)]~times[(l-25):(l)]))$coefficients[2,1]
    r.squared[j] = summary(lm(v[(l-25):(l)]~times[(l-25):(l)]))$r.squared
  }
for (i in 1:length(doubling.time)){
  doubling.time[i] = 1/doubling.time[i]/60
}



#make first deriv plots
num.plots <- 96
my.plots <- vector(num.plots, mode='list')

for (i in 1:num.plots) {
  plot(times[1:tp], first.deriv[i,1:tp],type = "o",main = row.names(con.data2)[i],col=ifelse(times==Tmax[i], "red", "dark green"),
       pch=ifelse(times==Tmax[i], 19, 19), cex=ifelse(times==Tmax[i], 1.5, 0.75), ylab = "Window Slope", xlab = "Time[s]")
  my.plots[[i]] <- recordPlot()
}
graphics.off()

pdf(paste(sample,'first_deriv.pdf', sep = "_"), onefile=TRUE)
for (my.plot in my.plots) {
  replayPlot(my.plot)
}
graphics.off()

#make growth deriv plots
num.plots <- 96
my.plots <- vector(num.plots, mode='list')

for (i in 1:num.plots) {
  plot(times[1:ncol(con.data2)], con.data2[i,],type = "o",main = row.names(con.data2)[i],col=ifelse(times==Tmax[i], "red", "dark green"),
       pch=ifelse(times==Tmax[i], 19, 19), cex=ifelse(times==Tmax[i], 1.5, 0.75),ylab = "O.D", xlab = "Time[s]")
  my.plots[[i]] <- recordPlot()
}
graphics.off()

pdf(paste(sample,'growth_curve.pdf', sep = "_"), onefile=TRUE)
for (my.plot in my.plots) {
  replayPlot(my.plot)
}
graphics.off()

write.csv(Tmax, paste(sample, "Tmax.csv", sep = "_"))
write.csv(doubling.time, paste(sample, "doubling_time.csv", sep = "_"))
