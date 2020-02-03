mod1 <- lm(y1~x1,data=anscombe)
mod2 <- lm(y2~x2,data=anscombe)
mod3 <- lm(y3~x3,data=anscombe)
mod4 <- lm(y4~x4,data=anscombe)

for(i in 1:4){
  y=getElement(anscombe,paste("y",i,sep=""))
  x=getElement(anscombe,paste("x",i,sep=""))
  assign(paste("mod",i,sep="_"),lm(y~x))
  rm(y,x)
}

ls()
coeff_matrix=NULL
for(i in 1:4){
  coeff_matrix <- rbind(coeff_matrix,coef(get(paste("mod",i,sep=""))))
}
coeff_matrix

R2_vect = NULL
for(i in 1:4){
  R2_vect <- rbind(R2_vect,summary(get(paste("mod",i,sep="")))$r.squared)
}

r_vect = NULL
for(i in 1:4){
  r_vect <- rbind(r_vect,cor(getElement(anscombe,paste("x",i,sep="")),
                             getElement(anscombe,paste("y",i,sep=""))))
}

r_vect^2==R2_vect
all.equal(r_vect^2,R2_vect)
#numeric ≥ 0. Differences smaller than tolerance are not reported. 
#The default value is close to 1.5e-8.

layout(matrix(1:4,nrow=2,ncol=2))
for(i in 1:4){
  plot(getElement(anscombe,paste("x",i,sep="")),
       getElement(anscombe,paste("y",i,sep="")),
       ylab=paste("y",i,sep=""), xlab=paste("x",i,sep=""))
  title(main=paste("RLS de",paste("y",i,sep=""),"p/r à",paste("x",i,sep="")),sub=paste("mod",i,sep=""))
  abline(get(paste("mod",i,sep="")),col="red",lwd=2,lty=2)
}



