setwd("D:/geoestadistica")
source("programas estadistica en R.txt")
require(ape)
require(vegan)
require(MASS)
require(akima)
require(gstat)
require(geoR)
require(lattice)
require(maptools)
require(rgdal)
require(GeoXp)
require(car)
require(Rcmdr)
trellis.par.set(sp.theme())
dir()

poligono=readShapePoly('lote.shp')
datos=read.table("datos.txt",sep="\t",dec=".",header=T)
xy=SpatialPoints(datos[c("x","y")])
plot(poligonos,axes=T)
text(datos$x,datos$y,datos$z,cex=0.6,col="blue")

#####Varianza constante
with(datos, tapply(z, Region, var, na.rm=TRUE))
leveneTest(z ~ Region, data=datos, center="median")

#####Media constante
with(datos, plotMeans(z, Region, error.bars="conf.int", level=0.95))
with(datos, tapply(z, Region, median, na.rm=TRUE))
kruskal.test(z ~ Region, data=datos)


#####Prueba de Mantel
w=as.matrix(dist(datos[,-3]))
u=as.matrix((dist(datos[,3]))^2)
mantel.test(w,u,graph=T)
mantel(w,u)


datos1=datos
coordinates(datos1)=c("x","y")
windows()
trellis.par.set(sp.theme())
driftmap(datos1,'z')

trellis.par.set(sp.theme())
spplot(datos1,"z")
intz=interp(x=datos$x,y=datos$y,z=datos$z)
contour(intz$z)
levelplot(intz$z)

muestra=as.data.frame(spsample(poligono,type='regular',n=10000))
names(muestra)=c("x","y")
##Trend surface
mod1reg=lm(z~x+y+I(x*y)+I(x^2)+I(y^2)+I(x*(y^2))+I((x^2)*y)+I((x^2)*(y^2))+I(x^3)+I(y^3),data=datos)
summary(mod1reg)
mod2reg=stepAIC(mod1reg,scope=list(upper=mod1reg$formula,lower=~1),direction="both")
summary(mod2reg)
mod3reg=lm(z~x+y+I(x*y)+I(x^2)+I(y^2),data=datos)
summary(mod3reg)
mod4reg=stepAIC(mod3reg,scope=list(upper=mod3reg$formula,lower=~1),direction="both")
summary(mod4reg)
ts3z=predict(mod2reg,muestra)
ts2z=predict(mod4reg,muestra)

#####Distancia inversa ponderada
idp=seq(0,20,0.5)
emc=numeric(0)
for(i in 1:length(idp)){
  zidw.cv <- gstat.cv(gstat(formula = z ~ 1, data = datos1, set = list(idp = idp[i])))
  emc[i]=mean((zidw.cv$residual^2))
}
cbind(idp,emc)
(idpmin=idp[which.min(emc)])
gridded(muestra)=c("x","y")
x11()
trellis.par.set(sp.theme())
idwz=idw(z~1,datos1,muestra,idp=idpmin)
predz=idwz[,-2]
names(predz)=c("IDW")
predz$TrendSurface2ndOrder=ts2z
predz$TrendSurface3rdOrder=ts3z
li=list("sp.polygons",poligono)
pts=list("sp.points",xy,pch=3,col="black",cex=0.2)
spplot(predz,c('TrendSurface2ndOrder','TrendSurface3rdOrder','IDW'),as.table=F,main="",sp.layout=list(li,pts),scales=list(draw=T),contour=T,labels=F,pretty=TRUE,col="black",key.space="right",col.regions=topo.colors(100))

#####Ajuste de medianas
zmat=tapply(datos$z,list(factor(datos$y),factor(datos$x)),function(x)x)
zmat
zmp=medpolish(zmat,na.rm=T)
zmp
roweff=matrix(rep(zmp$row,length(zmp$col)),ncol=length(zmp$col),byrow=F)
roweff
coleff=matrix(rep(zmp$col,length(zmp$row)),ncol=length(zmp$col),byrow=T)
coleff
modmp=zmp$overall+roweff+coleff
modmp
aux=cbind(as.numeric(zmp$residuals),as.numeric(modmp),as.numeric(roweff),as.numeric(coleff))
aux=aux[!is.na(aux)[,1],]
zmedpol=data.frame(datos,resMP=aux[,1],trend=aux[,2],roweff=aux[,3],coleff=aux[,4])
x11()
ztrend=tapply(zmedpol$trend,list(factor(zmedpol$y),factor(zmedpol$x)),function(x)x)
zmin=min(min(datos$z),min(zmedpol$trend))
zmax=max(max(datos$z),max(zmedpol$trend))
op=par(mfrow=c(1,2))
image(t(zmat),zlim=c(zmin,zmax),main="z")
image(t(ztrend),zlim=c(zmin,zmax),main="Ajustados por suavizamiento de medianas")
par(op)

###Variograma
geo=as.geodata(datos,coords.col=1:2,data.col=3)
x11()
trellis.par.set(sp.theme())
plot(geo,scatter3d=T)
trend=trend.spatial(formula(mod2reg),datos)
var1=variog(geo,trend=trend)
var2=variog(geo,trend=trend,max.dist=4)
op=par(mfrow=c(1,2))
plot(var1)
plot(var2)
par(op)
plot(variog4(geo,trend=trend,max.dist=4))
ev=eyefit(var2)
ev
####PAR¡METROS SELECCIONADOS########### 
> ev
cov.model sigmasq  phi tausq kappa kappa2   practicalRange
1         exponential    1.82 0.43  2.74  <NA>   <NA> 1.28817240939149
2 powered.exponential    2.08 0.78  2.74   0.5   <NA> 7.00005406685125
3           spherical    1.95 1.81  2.74  <NA>   <NA>             1.81
4            gaussian    0.91 1.21  3.78  <NA>   <NA> 2.09429000884074
#######################################

mod1exp=variofit(var2,ini=ev[[1]]$cov.pars,nugget=ev[[1]]$nugget,fix.nugget=F,cov.model=ev[[1]]$cov.model,weights="equal")
mod2exp=variofit(var2,ini=ev[[1]]$cov.pars,nugget=ev[[1]]$nugget,fix.nugget=F,cov.model=ev[[1]]$cov.model,weights="npairs")
mod3exp=variofit(var2,ini=ev[[1]]$cov.pars,nugget=ev[[1]]$nugget,fix.nugget=F,cov.model=ev[[1]]$cov.model,weights="cressie")
mod4exp=likfit(geo,trend=trend,ini=ev[[1]]$cov.pars,nugget=ev[[1]]$nugget,fix.nugget=F,cov.model=ev[[1]]$cov.model,lik.method = "ML")
mod5exp=likfit(geo,trend=trend,ini=ev[[1]]$cov.pars,nugget=ev[[1]]$nugget,fix.nugget=F,cov.model=ev[[1]]$cov.model,lik.method = "REML")

mod1powexp=variofit(var2,ini=ev[[2]]$cov.pars,nugget=ev[[2]]$nugget,fix.nugget=F,kappa=ev[[2]]$kappa,fix.kappa=F,cov.model=ev[[2]]$cov.model,weights="equal")
mod2powexp=variofit(var2,ini=ev[[2]]$cov.pars,nugget=ev[[2]]$nugget,fix.nugget=F,kappa=ev[[2]]$kappa,fix.kappa=F,cov.model=ev[[2]]$cov.model,weights="npairs")
mod3powexp=variofit(var2,ini=ev[[2]]$cov.pars,nugget=ev[[2]]$nugget,fix.nugget=F,kappa=ev[[2]]$kappa,fix.kappa=F,cov.model=ev[[2]]$cov.model,weights="cressie")
mod4powexp=likfit(geo,trend=trend,ini=ev[[2]]$cov.pars,nugget=ev[[2]]$nugget,fix.nugget=F,kappa=ev[[2]]$kappa,fix.kappa=F,cov.model=ev[[2]]$cov.model,lik.method = "ML")
mod5powexp=likfit(geo,trend=trend,ini=ev[[2]]$cov.pars,nugget=ev[[2]]$nugget,fix.nugget=F,kappa=ev[[2]]$kappa,fix.kappa=F,cov.model=ev[[2]]$cov.model,lik.method = "REML")

mod1sph=variofit(var2,ini=ev[[3]]$cov.pars,nugget=ev[[3]]$nugget,fix.nugget=F,cov.model=ev[[3]]$cov.model,weights="equal")
mod2sph=variofit(var2,ini=ev[[3]]$cov.pars,nugget=ev[[3]]$nugget,fix.nugget=F,cov.model=ev[[3]]$cov.model,weights="npairs")
mod3sph=variofit(var2,ini=ev[[3]]$cov.pars,nugget=ev[[3]]$nugget,fix.nugget=F,cov.model=ev[[3]]$cov.model,weights="cressie")
mod4sph=likfit(geo,trend=trend,ini=ev[[3]]$cov.pars,nugget=ev[[3]]$nugget,fix.nugget=F,cov.model=ev[[3]]$cov.model,lik.method = "ML")
mod5sph=likfit(geo,trend=trend,ini=ev[[3]]$cov.pars,nugget=ev[[3]]$nugget,fix.nugget=F,cov.model=ev[[3]]$cov.model,lik.method = "REML")

mod1gau=variofit(var2,ini=ev[[4]]$cov.pars,nugget=ev[[4]]$nugget,fix.nugget=F,cov.model=ev[[4]]$cov.model,weights="equal")
mod2gau=variofit(var2,ini=ev[[4]]$cov.pars,nugget=ev[[4]]$nugget,fix.nugget=F,cov.model=ev[[4]]$cov.model,weights="npairs")
mod3gau=variofit(var2,ini=ev[[4]]$cov.pars,nugget=ev[[4]]$nugget,fix.nugget=F,cov.model=ev[[4]]$cov.model,weights="cressie")
mod4gau=likfit(geo,trend=trend,ini=ev[[4]]$cov.pars,nugget=ev[[4]]$nugget,fix.nugget=F,cov.model=ev[[4]]$cov.model,lik.method = "ML")
mod5gau=likfit(geo,trend=trend,ini=ev[[4]]$cov.pars,nugget=ev[[4]]$nugget,fix.nugget=F,cov.model=ev[[4]]$cov.model,lik.method = "REML")


op=par(mfrow=c(2,2))
trellis.par.set(sp.theme())
plot(var2,main="Powered Exponential",ylim=c(0,5))
lines(mod1powexp,col=1)
lines(mod2powexp,col=2)
lines(mod3powexp,col=3)
lines(mod4powexp,col=4)
legend(1.2,3, legend = c("OLS", "WLS - npairs", "WLS - cressie", "ML", "REML"), col = 1:5, lty = c(1, 1, 2, 2), lwd = c(1, 2, 1, 2), cex = 0.9)

plot(var2,main="Spherical",ylim=c(0,5))
lines(mod1sph,col=1)
lines(mod2sph,col=2)
lines(mod3sph,col=3)
lines(mod4sph,col=4)
lines(mod5sph,col=5)
legend(1.2,3, legend = c("OLS", "WLS - npairs", "WLS - cressie", "ML", "REML"), col = 1:5, lty = c(1, 1, 2, 2), lwd = c(1, 2, 1, 2), cex = 0.9)

plot(var2,main="Gaussian",ylim=c(0,5))
lines(mod1gau,col=1)
lines(mod2gau,col=2)
lines(mod3gau,col=3)
lines(mod4gau,col=4)
lines(mod5gau,col=5)
legend(1.2,3, legend = c("OLS", "WLS - npairs", "WLS - cressie", "ML", "REML"), col = 1:5, lty = c(1, 1, 2, 2), lwd = c(1, 2, 1, 2), cex = 0.9)

plot(var2,main="Exponential",ylim=c(0,5))
lines(mod1exp,col=1)
lines(mod2exp,col=2)
lines(mod3exp,col=3)
lines(mod4exp,col=4)
lines(mod5exp,col=5)
legend(1.2,3, legend = c("OLS", "WLS - npairs", "WLS - cressie", "ML", "REML"), col = 1:5, lty = c(1, 1, 2, 2), lwd = c(1, 2, 1, 2), cex = 0.9)
par(op)

###OPCI”N geoR
cvolspowexp=xvalid(geo,model=mod1powexp,trend=trend,max.dist=4)
cvwlsnpairspowexp=xvalid(geo,model=mod2powexp,trend=trend,max.dist=4)
cvwlscressiepowexp=xvalid(geo,model=mod3powexp,trend=trend,max.dist=4)
cvmlpowexp=xvalid(geo,model=mod4powexp,trend=trend,max.dist=4)

cvolssph=xvalid(geo,model=mod1sph,trend=trend,max.dist=4)
cvwlsnpairssph=xvalid(geo,model=mod2sph,trend=trend,max.dist=4)
cvwlscressiesph=xvalid(geo,model=mod3sph,trend=trend,max.dist=4)
cvmlsph=xvalid(geo,model=mod4sph,trend=trend,max.dist=4)
cvremlsph=xvalid(geo,model=mod5sph,trend=trend,max.dist=4)

cvolsgau=xvalid(geo,model=mod1gau,trend=trend,max.dist=4)
cvwlsnpairsgau=xvalid(geo,model=mod2gau,trend=trend,max.dist=4)
cvwlscressiegau=xvalid(geo,model=mod3gau,trend=trend,max.dist=4)
cvmlgau=xvalid(geo,model=mod4gau,trend=trend,max.dist=4)
cvremlgau=xvalid(geo,model=mod5gau,trend=trend,max.dist=4)

cvolsexp=xvalid(geo,model=mod1exp,trend=trend,max.dist=4)
cvwlsnpairsexp=xvalid(geo,model=mod2exp,trend=trend,max.dist=4)
cvwlscressieexp=xvalid(geo,model=mod3exp,trend=trend,max.dist=4)
cvmlexp=xvalid(geo,model=mod4exp,trend=trend,max.dist=4)
cvremlexp=xvalid(geo,model=mod5exp,trend=trend,max.dist=4)

sqrt((mean(cvolspowexp$error)^2)+(mean(cvolspowexp$std.error)^2)+((var(cvolspowexp$std.error)-1)^2))
sqrt((mean(cvwlsnpairspowexp$error)^2)+(mean(cvwlsnpairspowexp$std.error)^2)+((var(cvwlsnpairspowexp$std.error)-1)^2))
sqrt((mean(cvwlscressiepowexp$error)^2)+(mean(cvwlscressiepowexp$std.error)^2)+((var(cvwlscressiepowexp$std.error)-1)^2))
sqrt((mean(cvmlpowexp$error)^2)+(mean(cvmlpowexp$std.error)^2)+((var(cvmlpowexp$std.error)-1)^2))

sqrt((mean(cvolsmat$error)^2)+(mean(cvolsmat$std.error)^2)+((var(cvolsmat$std.error)-1)^2))
sqrt((mean(cvwlsnpairsmat$error)^2)+(mean(cvwlsnpairsmat$std.error)^2)+((var(cvwlsnpairsmat$std.error)-1)^2))
sqrt((mean(cvwlscressiemat$error)^2)+(mean(cvwlscressiemat$std.error)^2)+((var(cvwlscressiemat$std.error)-1)^2))
sqrt((mean(cvmlmat$error)^2)+(mean(cvmlmat$std.error)^2)+((var(cvmlmat$std.error)-1)^2))
sqrt((mean(cvremlmat$error)^2)+(mean(cvremlmat$std.error)^2)+((var(cvremlmat$std.error)-1)^2))

sqrt((mean(cvolssph$error)^2)+(mean(cvolssph$std.error)^2)+((var(cvolssph$std.error)-1)^2))
sqrt((mean(cvwlsnpairssph$error)^2)+(mean(cvwlsnpairssph$std.error)^2)+((var(cvwlsnpairssph$std.error)-1)^2))
sqrt((mean(cvwlscressiesph$error)^2)+(mean(cvwlscressiesph$std.error)^2)+((var(cvwlscressiesph$std.error)-1)^2))
sqrt((mean(cvmlsph$error)^2)+(mean(cvmlsph$std.error)^2)+((var(cvmlsph$std.error)-1)^2))
sqrt((mean(cvremlsph$error)^2)+(mean(cvremlsph$std.error)^2)+((var(cvremlsph$std.error)-1)^2))

sqrt((mean(cvolsgau$error)^2)+(mean(cvolsgau$std.error)^2)+((var(cvolsgau$std.error)-1)^2))
sqrt((mean(cvwlsnpairsgau$error)^2)+(mean(cvwlsnpairsgau$std.error)^2)+((var(cvwlsnpairsgau$std.error)-1)^2))
sqrt((mean(cvwlscressiegau$error)^2)+(mean(cvwlscressiegau$std.error)^2)+((var(cvwlscressiegau$std.error)-1)^2))
sqrt((mean(cvmlgau$error)^2)+(mean(cvmlgau$std.error)^2)+((var(cvmlgau$std.error)-1)^2))
sqrt((mean(cvremlgau$error)^2)+(mean(cvremlgau$std.error)^2)+((var(cvremlgau$std.error)-1)^2))

sqrt((mean(cvolsexp$error)^2)+(mean(cvolsexp$std.error)^2)+((var(cvolsexp$std.error)-1)^2))
sqrt((mean(cvwlsnpairsexp$error)^2)+(mean(cvwlsnpairsexp$std.error)^2)+((var(cvwlsnpairsexp$std.error)-1)^2))
sqrt((mean(cvwlscressieexp$error)^2)+(mean(cvwlscressieexp$std.error)^2)+((var(cvwlscressieexp$std.error)-1)^2))
sqrt((mean(cvmlexp$error)^2)+(mean(cvmlexp$std.error)^2)+((var(cvmlexp$std.error)-1)^2))
sqrt((mean(cvremlexp$error)^2)+(mean(cvremlexp$std.error)^2)+((var(cvremlexp$std.error)-1)^2))

muestra=spsample(poligono,n=10000,"regular")
muestra1=as.data.frame(muestra)
names(muestra1)=c("x","y")
predUK=krige.conv(geo,loc=muestra1,krige=krige.control(type.krige = "ok",obj.model = mod5gau,trend.l=trend.spatial(~y+I(x*y)+I(x^2)+I(y^2)+I((x^2)*y)+I(x^3)+I(y^3),muestra1),trend.d=trend))
z.uk=data.frame(muestra1,var1.pred=predUK$predict,var1.var=predUK$krige.var)
gridded(z.uk)=c("x","y")
summary(z.uk)
li=list("sp.polygons",poligono)
pts=list("sp.points",xy,pch=3,col="black",cex=0.2)
spplot(z.uk,c("var1.pred"),as.table=TRUE,main="",sp.layout=list(li,pts),contour=TRUE,labels=FALSE,pretty=TRUE,col="black",col.regions=topo.colors(100))
x11()
spplot(z.uk,c("var1.var"),as.table=TRUE,main="",sp.layout=list(li,pts),contour=TRUE,labels=FALSE,pretty=TRUE,col="black",col.regions=topo.colors(100))

###OpciÛn gstat
(ve.fit1powexp <- vgm(mod1powexp$cov.pars[1],"Exc",mod1powexp$cov.pars[2],mod1powexp$nugget,kappa=mod1powexp$kappa))
(ve.fit2powexp <- vgm(mod2powexp$cov.pars[1],"Exc",mod2powexp$cov.pars[2],mod2powexp$nugget,kappa=mod2powexp$kappa))
(ve.fit3powexp <- vgm(mod3powexp$cov.pars[1],"Exc",mod3powexp$cov.pars[2],mod3powexp$nugget,kappa=mod3powexp$kappa))
(ve.fit4powexp <- vgm(mod4powexp$cov.pars[1],"Exc",mod4powexp$cov.pars[2],mod4powexp$nugget,kappa=mod4powexp$kappa))
(ve.fit5powexp <- vgm(mod5powexp$cov.pars[1],"Exc",mod5powexp$cov.pars[2],mod5powexp$nugget,kappa=mod5powexp$kappa))


ve.fit1mat <- as.vgm.variomodel(mod1mat)
ve.fit2mat <- as.vgm.variomodel(mod2mat)
ve.fit3mat <- as.vgm.variomodel(mod3mat)
ve.fit4mat <- as.vgm.variomodel(mod4mat)
ve.fit5mat <- as.vgm.variomodel(mod5mat)

ve.fit1sph <- as.vgm.variomodel(mod1sph)
ve.fit2sph <- as.vgm.variomodel(mod2sph)
ve.fit3sph <- as.vgm.variomodel(mod3sph)
ve.fit4sph <- as.vgm.variomodel(mod4sph)
ve.fit5sph <- as.vgm.variomodel(mod5sph)

ve.fit1gau <- as.vgm.variomodel(mod1gau)
ve.fit2gau <- as.vgm.variomodel(mod2gau)
ve.fit3gau <- as.vgm.variomodel(mod3gau)
ve.fit4gau <- as.vgm.variomodel(mod4gau)
ve.fit5gau <- as.vgm.variomodel(mod5gau)

ve.fit1exp <- as.vgm.variomodel(mod1exp)
ve.fit2exp <- as.vgm.variomodel(mod2exp)
ve.fit3exp <- as.vgm.variomodel(mod3exp)
ve.fit4exp <- as.vgm.variomodel(mod4exp)
ve.fit5exp <- as.vgm.variomodel(mod5exp)

cvolspowexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit1powexp, maxdist = 140)
cvwlsnpairspowexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit2powexp, maxdist = 140)
cvwlsncressiepowexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit3powexp, maxdist = 140)
cvmlpowexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit4powexp, maxdist = 140)
cvremlpowexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit5powexp, maxdist = 140)

cvolsmat <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit1mat, maxdist = 140)
cvwlsnpairsmat <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit2mat, maxdist = 140)
cvwlsncressiemat <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit3mat, maxdist = 140)
cvmlmat <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit4mat, maxdist = 140)
cvremlmat <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit5mat, maxdist = 140)

cvolssph <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit1sph, maxdist = 140)
cvwlsnpairssph <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit2sph, maxdist = 140)
cvwlsncressiesph <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit3sph, maxdist = 140)
cvmlsph <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit4sph, maxdist = 140)
cvremlsph <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit5sph, maxdist = 140)

cvolsgau <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit1gau, maxdist = 140)
cvwlsnpairsgau <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit2gau, maxdist = 140)
cvwlsncressiegau <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit3gau, maxdist = 140)
cvmlgau <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit4gau, maxdist = 140)
cvremlgau <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit5gau, maxdist = 140)

cvolsexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit1exp, maxdist = 140)
cvwlsnpairsexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit2exp, maxdist = 140)
cvwlsncressieexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit3exp, maxdist = 140)
cvmlexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit4exp, maxdist = 140)
cvremlexp <- krige.cv(head ~ x+I(x*y)+I(y^2), datos1, ve.fit5exp, maxdist = 140)

sqrt((mean(cvolspowexp$residual))^2+(mean(cvolspowexp$zscore))^2+(var(cvolspowexp$zscore)-1)^2)
sqrt((mean(cvwlsnpairspowexp$residual))^2+(mean(cvwlsnpairspowexp$zscore))^2+(var(cvwlsnpairspowexp$zscore)-1)^2)
sqrt((mean(cvwlsncressiepowexp$residual))^2+(mean(cvwlsncressiepowexp$zscore))^2+(var(cvwlsncressiepowexp$zscore)-1)^2)
sqrt((mean(cvmlpowexp$residual))^2+(mean(cvmlpowexp$zscore))^2+(var(cvmlpowexp$zscore)-1)^2)
sqrt((mean(cvremlpowexp$residual))^2+(mean(cvremlpowexp$zscore))^2+(var(cvremlpowexp$zscore)-1)^2)

sqrt((mean(cvolsmat$residual))^2+(mean(cvolsmat$zscore))^2+(var(cvolsmat$zscore)-1)^2)
sqrt((mean(cvwlsnpairsmat$residual))^2+(mean(cvwlsnpairsmat$zscore))^2+(var(cvwlsnpairsmat$zscore)-1)^2)
sqrt((mean(cvwlsncressiemat$residual))^2+(mean(cvwlsncressiemat$zscore))^2+(var(cvwlsncressiemat$zscore)-1)^2)
sqrt((mean(cvmlmat$residual))^2+(mean(cvmlmat$zscore))^2+(var(cvmlmat$zscore)-1)^2)
sqrt((mean(cvremlmat$residual))^2+(mean(cvremlmat$zscore))^2+(var(cvremlmat$zscore)-1)^2)

sqrt((mean(cvolssph$residual))^2+(mean(cvolssph$zscore))^2+(var(cvolssph$zscore)-1)^2)
sqrt((mean(cvwlsnpairssph$residual))^2+(mean(cvwlsnpairssph$zscore))^2+(var(cvwlsnpairssph$zscore)-1)^2)
sqrt((mean(cvwlsncressiesph$residual))^2+(mean(cvwlsncressiesph$zscore))^2+(var(cvwlsncressiesph$zscore)-1)^2)
sqrt((mean(cvmlsph$residual))^2+(mean(cvmlsph$zscore))^2+(var(cvmlsph$zscore)-1)^2)
sqrt((mean(cvremlsph$residual))^2+(mean(cvremlsph$zscore))^2+(var(cvremlsph$zscore)-1)^2)

sqrt((mean(cvolsgau$residual))^2+(mean(cvolsgau$zscore))^2+(var(cvolsgau$zscore)-1)^2)
sqrt((mean(cvwlsnpairsgau$residual))^2+(mean(cvwlsnpairsgau$zscore))^2+(var(cvwlsnpairsgau$zscore)-1)^2)
sqrt((mean(cvwlsncressiegau$residual))^2+(mean(cvwlsncressiegau$zscore))^2+(var(cvwlsncressiegau$zscore)-1)^2)
sqrt((mean(cvmlgau$residual))^2+(mean(cvmlgau$zscore))^2+(var(cvmlgau$zscore)-1)^2)
sqrt((mean(cvremlgau$residual))^2+(mean(cvremlgau$zscore))^2+(var(cvremlgau$zscore)-1)^2)

sqrt((mean(cvolsexp$residual))^2+(mean(cvolsexp$zscore))^2+(var(cvolsexp$zscore)-1)^2)
sqrt((mean(cvwlsnpairsexp$residual))^2+(mean(cvwlsnpairsexp$zscore))^2+(var(cvwlsnpairsexp$zscore)-1)^2)
sqrt((mean(cvwlsncressieexp$residual))^2+(mean(cvwlsncressieexp$zscore))^2+(var(cvwlsncressieexp$zscore)-1)^2)
sqrt((mean(cvmlexp$residual))^2+(mean(cvmlexp$zscore))^2+(var(cvmlexp$zscore)-1)^2)
sqrt((mean(cvremlexp$residual))^2+(mean(cvremlexp$zscore))^2+(var(cvremlexp$zscore)-1)^2)

muestra=spsample(poligono,n=10000,"regular")
muestra1=as.data.frame(muestra)
names(muestra1)=c("x","y")
gridded(muestra1)=c("x","y")
head.uk=krige(head~x+I(x*y)+I(y^2),datos1,muestra1,ve.fit3gau,maxdist=140)
summary(head.uk)
li=list("sp.polygons",poligono)
pts=list("sp.points",xy,pch=3,col="black",cex=0.2)
spplot(head.uk,c("var1.pred"),as.table=TRUE,main="",sp.layout=list(li,pts),contour=TRUE,labels=FALSE,pretty=TRUE,col="black",col.regions=topo.colors(100))
x11()
spplot(head.uk,c("var1.var"),as.table=TRUE,main="",sp.layout=list(li,pts),contour=TRUE,labels=FALSE,pretty=TRUE,col="black",col.regions=topo.colors(100))


