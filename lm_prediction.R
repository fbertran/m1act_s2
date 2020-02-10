install.packages("faraway")
library(faraway)
help(package="faraway")

load('~/Documents/Enseignements/Stat_L3_Master/2018/Modele_Regression/vulnerability.RData')
str(vul)

fit_univ = lm(ln_death_risk~ln_events, data=vul,x = TRUE,y=TRUE)
predict(fit_univ)

newdata=data.frame(ln_events=3.4)
pred=predict(fit_univ,newdata,interval="predict")
print(pred)

ic=predict(fit_univ,interval="confidence")

X = fit_univ$x
Y = fit_univ$y
Xplus = matrix(c(1,3.4),nrow=1)

BetaChapeau = solve(t(X)%*%X)%*%t(X)%*%Y

YplusPred = Xplus%*%BetaChapeau

qt(0.975,nrow(X)-1-1) ; #qt(0.025,nrow(X)-1-1)
CMres = (anova(fit_univ)$`Mean Sq`)[2]

c(YplusPred,
YplusPred-sqrt(CMres*(Xplus%*%solve(t(X)%*%X)%*%t(Xplus)+1))*qt(0.975,nrow(X)-1-1),
YplusPred+sqrt(CMres*(Xplus%*%solve(t(X)%*%X)%*%t(Xplus)+1))*qt(0.975,nrow(X)-1-1))

pred


conf=predict(fit_univ,newdata,interval="confidence")
c(YplusPred,
  YplusPred-sqrt(CMres*(Xplus%*%solve(t(X)%*%X)%*%t(Xplus)))*qt(0.975,nrow(X)-1-1),
  YplusPred+sqrt(CMres*(Xplus%*%solve(t(X)%*%X)%*%t(Xplus)))*qt(0.975,nrow(X)-1-1))

install.packages('HH')
library(HH)
ci.plot(fit_univ)


X = vul[,c(3:6)]
cor_mat = cor(X)

#Calcul des valeurs propres et vecteurs propres
propres = eigen(cor_mat)
propres$values[1] / propres$values



