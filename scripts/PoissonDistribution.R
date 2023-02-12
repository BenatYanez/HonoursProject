#Modified from Matthew Fellbaum
#Last modified 10/02/2023
#Put the name of the file on "file name"
read.table ("file name.ldepth.mean", header=T)
#Open up the Depth file and put it into a R compatible table, put the name of the file on "file name"
depth_table <-read.table ("file name.ldepth.mean", header=T)
#Calculate mean,var and density
depth <- density (depth_table$MEAN_DEPTH)
mean <- mean (depth_table$MEAN_DEPTH)
var <- var (depth_table$MEAN_DEPTH)
qpois (c(0.0025, 0.9975),mean)

#Plot the mean site depth data
plot(depth,frame =FALSE, col ="red",
     main = "Density plot of mean site depth of the entire genome",
     xlim = c(10,80),
     xlab = "Mean site depth" )
#Plot a line showing the Poisson curve
lines (dpois (seq(0, max(depth_table$MEAN_DEPTH),1),mean),
       col = "green", lty =2 )
#Add vertical lines representing the 99.5% confidence intervals and add a legend for the 3 lines
abline (v = qpois(c(0.0025,0.9975), mean), col=c("blue", "blue"),
        lty =3:3)
legend (50, 0.06, legend=c("Actual Data", "Poisson Curve", "qpois 99.5%"),
        col=c("red", "green", "blue"), lty=1:2, cex=0.8)
#Create a density plot
plot(density(dpois(depth_table$MEAN_DEPTH,
                   mean(depth_table$MEAN_DEPTH))),
     xlim=c(-0.002,0.002))
plot(dpois(depth_table$MEAN_DEPTH, mean(depth_table$MEAN_DEPTH)))
plot(dpois(depth,mean))


