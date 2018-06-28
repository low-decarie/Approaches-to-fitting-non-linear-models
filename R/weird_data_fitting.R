#Set random number seed so that fits have same solution ###
set.seed(18)

#Create data  ###
x <- 1:100
weird_data <- data.frame(x=x,
                         y=(x*(sin(0.2*x)+1))^2)
plot(weird_data, type="l")


#Save weird data
save(weird_data, file="./Data/weird_data.RData")
write.csv(weird_data, file="./Data/weird_data.csv")

#Step 1 transform  ###
plot(sqrt(y)~x, data=weird_data,type="l")

#Step 2 fit with quadratic  ###
weird_lm_fit <- lm(sqrt(y)~x+I(x^2), data=weird_data)
plot(weird_data)
lines(weird_data$x,predict(weird_lm_fit)^2)

#Step 3 segmentize the quadratic  ###
weird_lm_fit_seg <- segmented(weird_lm_fit,
                              psi=c(30,45,70,80,90),
                              control=seg.control(it.max=100,
                                                  n.boot=100))
plot(y~x, data=weird_data)
lines(weird_data$x,predict(weird_lm_fit_seg)^2)
summary(weird_lm_fit_seg)

#Run step 3 to see how solutions change due to 