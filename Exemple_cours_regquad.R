Taille <- rnorm(30)*10+160
Taille

Masse <- (Taille + rnorm(30,0,2))/3

plot(Taille,Masse)

cor(Taille,Masse)

lm(Masse~Taille)
abline(lm(Masse~Taille),col="red")

Masse2 <- Taille+(Taille-mean(Taille))^2+rnorm(30,0,10)
plot(Taille,Masse2)
cor(Taille,Masse2)

plot(Taille,Masse2)
modelquad <- lm(Masse2~Taille+I(Taille^2))
modelquad
predict(modelquad)
points(Taille,predict(modelquad),col="blue",pch="+")
abline(lm(Masse2~Taille),col="red")
lines((1300:1800)/10,predict(modelquad,newdata = data.frame(Taille=(1300:1800)/10)),col="blue")

modeldeg1 <- lm(Masse2~Taille)
plot(Taille,residuals(modeldeg1))

plot(Taille,residuals(modelquad))

Xmatrix=cbind(Intercept=1,Taille,Taille2=Taille*Taille)
t(Xmatrix)%*%Xmatrix
solve(t(Xmatrix)%*%Xmatrix)%*%t(Xmatrix)%*%Masse2
modelquad
modelquad_withX=lm(Masse2~Taille+I(Taille^2),x = TRUE)
all(modelquad_withX$x==Xmatrix)

summary(modelquad)
#Residual standard error: 8.615 on 27 degrees of freedom

Xmatrix%*%solve(t(Xmatrix)%*%Xmatrix)%*%t(Xmatrix)%*%Masse2
#identique à 
predict(modelquad)

#(sigma^2)^chapeau
sigma2hat<-sum((Masse2-Xmatrix%*%solve(t(Xmatrix)%*%Xmatrix)%*%t(Xmatrix)%*%Masse2)^2)/(30-(2+1))

#sigma^chapeau=sqrt((sigma^2)^chapeau)
sigmahat<-sqrt(sigma2hat)
#identique à 
#Residual standard error: 8.615 on 27 degrees of freedom

summary(modelquad)
sigma2hat*solve(t(Xmatrix)%*%Xmatrix)

summary(modelquad)
sqrt(diag(sigma2hat*solve(t(Xmatrix)%*%Xmatrix)))
