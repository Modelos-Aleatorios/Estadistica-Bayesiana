
#Ejercicio 2.3
library(readr)
normaldata <- read_delim("normaldata", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
larcenies<-c(normaldata$X1,normaldata$X2,normaldata$X3,normaldata$X4,normaldata$X5)
hist(larcenies,breaks = 21)


#Ejercicio 2.4
library(readr)
CMBdata <- read_csv("CMBdata", col_names = FALSE)
g<-CMBdata$X1

image_matrix<-matrix(g,nrow = 800,ncol = 800,byrow = TRUE)
image(image_matrix, axes = FALSE, col = grey(seq(0, 1, length = 256)))

gen_cords_pairs<-function() {
  coords<-c(trunc(runif(4)*800))
  if(coords[1]>coords[3]) {
    temp<-coords[1]
    coords[1]<-coords[3]
    coords[3]<-temp
  }
  if(coords[2]>coords[4]) {
    temp<-coords[2]
    coords[2]<-coords[4]
    coords[4]<-temp
  }
  coords
}

coords<-gen_cords_pairs()

r_matrix<-image_matrix[coords[1]:coords[3],coords[2]:coords[4]]
image(r_matrix,axes = FALSE, col = grey(seq(0, 1, length = 256)))

image(r_matrix,axes = FALSE, col = grey(1:1000/1000))

hist(r_matrix,breaks=100,freq = FALSE)

normal_fit<-function(data) {
  xfit<-seq(min(data),max(data),length=100)
  yfit<-dnorm(xfit,mean=mean(data),sd=sd(data))
  lines(xfit, yfit, col="black", lwd=2)
}

normal_fit(r_matrix)
normal_fit(CMBdata$X1)


#Normal curve fit
hist(CMBdata$X1,breaks = 100, freq = FALSE)
normal_fit(CMBdata$X1)


