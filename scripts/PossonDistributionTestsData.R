read.table ("1Mbfiltered.ldepth.mean", header=T)
#Open up the Depth file and put it into a R compatible table
depth_table <-read.table ("1Mbfiltered.ldepth.mean", header=T)
#Calculate mean,var and density
depth <- density (depth_table$MEAN_DEPTH)
mean <- mean (depth_table$MEAN_DEPTH)
var <- var (depth_table$MEAN_DEPTH)
qpois (c(0.0025, 0.9975),mean)

#Plot the Poisson Distribution
plot(depth,frame =FALSE, col ="red",
     main = "Density plot of mean site depth of 1Mb data",
     xlim = c(10,80),
     xlab = "Mean site depth" )
lines (dpois (seq(0, max(depth_table$MEAN_DEPTH),1),mean),
       col = "green", lty =2 )
#Changed from original code from 0.005 to 0.0025 and from 0.995 to 0.9975 
abline (v = qpois(c(0.0025,0.9975), mean), col=c("blue", "blue"),
        lty =3:3)
legend (50, 0.06, legend=c("Actual Data", "Poisson Curve", "qpois 99.5%"),
        col=c("red", "green", "blue"), lty=1:2, cex=0.8)
#legend (50,0.06, legend=c("ATsweeps","qpois =0.005", "qpois=0.995"), col=c("red","green","blue"), lty=1:1, cex=0.8)
plot(density(dpois(depth_table$MEAN_DEPTH,
                   mean(depth_table$MEAN_DEPTH))),
     xlim=c(-0.002,0.002))
plot(dpois(depth_table$MEAN_DEPTH, mean(depth_table$MEAN_DEPTH)))
plot(dpois(depth,mean))


