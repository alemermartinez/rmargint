library(mgcv)
salida.gam <- gam(Volume~s(Height)+s(Girth,k=20), family=gaussian(),data=trees)
summary(salida.gam)
deviance(salida.gam)
plot(salida.gam)
sum(residuals(salida.gam)^2)


set.seed(123)
Volume <- trees$Volume
Height <- trees$Height
Girth <- trees$Girth
X <- cbind(Height,Girth)

library(rmargint)
salida.margint <- margint.cl(X,Volume,windows=c(sd(Height),sd(Girth)))
summary(salida.margint)
plot(salida.margint)

residuals(salida.gam)
residuals(salida.margint)
