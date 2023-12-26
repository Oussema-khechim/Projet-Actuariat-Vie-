library(demography)
library(rainbow)
library(graphics)
library(StMoMo)
library(survival)
library(fitdistrplus)




#Question1
donnee <- read.table(file = "C:/Users/oussama/Downloads/Human-mortality-dashboard-main/Human-mortality-dashboard-main/mx.txt", header = TRUE, fill = TRUE)
x <- data.matrix(donnee[, 1])
y <- data.matrix(donnee[, 2])
z <- data.matrix(donnee[, 4])
dffff <- data.frame(Year = x, Age = y, Male = z)
summary(dffff)

France1 <- read.demogdata(file="C:/Users/oussama/Downloads/Human-mortality-dashboard-main/Human-mortality-dashboard-main/mx.txt",popfile="C:/Users/oussama/Downloads/Human-mortality-dashboard-main/Human-mortality-dashboard-main/Exposures_1x1.txt",type="mortality", label="France")

#par(mar = rep(2, 4))
Fr_years <- 1950:2000
plot(France1, series = "male", years = Fr_years, datatype = "rate", main = "France, Homme, 1950 - 2018", xlab = "Age", ylab = "Taux (log)")
colfunc <- colorRampPalette(c("red", "blue"))

cohort <- function(year, rates, log=FALSE)
{
  xg <- diag(rates[, colnames(rates) >= year])
  names(xg) <- paste(year, rownames(rates)[1:length(xg)], sep="-x=")
  if(log)
    xg <- log(xg)
  xg
}
#QUESTION 2

plot(cohort(1950, France1$rate$male, log=T),
     col=colfunc(length(1950)), 
     type="l",
     ylim=c(-11,5),
     main="France: Cohortes",
     xlab = "age", 
     ylab = "Taux de Mortalite")


#question 5 

cohort1950_m <- cohort(1950, France1$rate$male)

plot(France1$age, log(France1$rate$male[,"1950"]), main ='log mortality rates (FR_male, 1950)',
     xlab = "Ages x", ylab = "log mortality rates", type = "l")

lines(0:(length(cohort1950_m)-1), log(cohort1950_m), main ='male log mortality rates (FR, 1950)',
      xlab = "Ages x", ylab = "log mortality rates", type = "l",col='red')

legend(-4, -0.5,legend = c("lecture longitudinale", "lecture cohorte"),
       col=c("black","red"),lty = 1, cex=0.7,
       box.lty = 0
)


#Intervalle de confiance

fit.norm<-fitdist(cohort(1950, France1$rate$total, log=T), "norm" )
fit.norm$estimate

#NB : Niveau de con???ance de 99% ===> z(alpha/2)=2.576
ect = fit.norm$estimate["sd"]
moy_emp = fit.norm$estimate["mean"]
IC_inf = moy_emp-2.576*ect/sqrt(2)
IC_sup = moy_emp+2.576*ect/sqrt(2)



  
  plot(cohort(1950, France1$rate$male, log=T),
       col=colfunc(length(1950)), 
       type="l",
       ylim=c(-11,5),
       main="FR: Cohorte ",
       xlab = "age", 
       ylab = "Taux de Mortalite")
  
  abline(h=moy_emp,col="blue", lwd=3, lty=2)
  abline(h=IC_inf,col="green", lwd=3, lty=2)
  abline(h=IC_sup,col="green", lwd=3, lty=2)
  
  
#Question 3 --Modèle Lee-Carter
 
  par(mfrow=c(1,1))
  FR_ages = c(0,10,20,30,40,50,60,70,80,90,100)
  
  plot(France1,
       series="total",
       datatype="rate", 
       plot.type="time",
       age = FR_ages,
       main="total male death rates (1950 - 2018) ",axes = F)
  # on fixe les axes comme suit :
  axis(side = 1, at=1816:2018)
  axis(side = 2, at=-10:0)
  legend(x="bottomright", legend = FR_ages,
         col = rainbow(length(FR_ages)*1.25), lty = 1, cex=0.6,
         box.lwd = 0.3)
  
  #choix du plage
  
  #plage d ages
  ages.fit = 0:100
  #periode de calibration
  years.fit = 1950:2018
 
  # fitting Lee Carter model :
  
  lca.total <- lca(France1_ls_m, series="total", adjust="dt",years =years.fit ,ages = ages.fit)
  lca.male <- lca(France1_ls_m, series="male", adjust="dt",years =years.fit ,ages = ages.fit)
  
  
    plot(lca.male$ax, main="Coef. ax sur donnees francaise", xlab="Age", ylab="ax", type="l")
    lines(x=lca.male$age, y=lca.male$ax, main="ax", lty=2)
    legend("bottomright","Male", cex=0.8,  lty=1:2)
  
    
## Analyse des paramètres :
    
    # - αx :  la valeur moyenne des logs de la mortalité instantanné ( ln µ( x t, ) au cours du temps  ) elle crois en fonction de l’age elle varie entre -7 et -1 .
    
  
    
    plot(lca.male$bx, main="Coef. bx sur donnees francaise", ylim=c(0,0.03),xlab="Age", ylab="bx", type="l")
    lines(x=lca.male$age, y=lca.male$bx, main="bx", lty=2)
    legend("bottomright","Male", cex=0.8,  lty=1:2)

 #- βx indique la sensibilité de la mortalité instantanée par rapport à l’évolution générale de la mortalité. Si on se situe à partir de 18 ans, on constate que les âges les plus sensibles à l’évolution temporelle de la mortalité sont ceux entre 70 et 80 ans . On atteint en effet des pics sur ces tranches d’âges.
    
  
    plot(lca.total$ax-lca.male$ax, main="Ecart avec population totale", xlab="Age x", ylab=expression(paste(Delta, " ax")), type="l" , col='green')
    lines(x=lca.male$age, y=lca.male$ax-lca.total$ax, main="delta", lty=2 , col ="blue")
    legend("topright","Male", cex=0.8, lty=1:2)

    
    plot(lca.total$kt, xlab="Year", main="Coef. kt sur donnees francaises",ylab="kt", type="l",ylim=c(-100, 100))
    lines(lca.male$year, y=lca.male$kt, main="kt", lty=2 , col="blue")
    legend("topright", "Male", cex=0.8, lty=1:2)

     # - D’après la figure ci-dessus et comme  kt indique l’évolution générale de la mortalité dans le temps ;
   # On constate une tendance linéaire à la décroissance Cette tendance à la décroissance du paramètre k, qui devient négatif au cours de la période, associée à la positivité moyenne du paramètre β implique d’après la formule de Lee-Carter, une diminution des taux instantanés de mortalité. En conséquence, on assiste à une augmentation  de la probabilité  de la  survie sur la période observée.
    
    

    FR.stmomo.t<-StMoMoData(data=France1_ls_m ,series = "male",type="central")
    #ajustement du model (fitiing) :
    LC1 <- lc(link = "logit" )
    LCfit1 <- fit(LC1, data = central2initial(FR.stmomo.t), ages.fit = ages.fit,  years.fit = years.fit)
    # parametre kt :
    plot(LCfit1$years,LCfit1$kt,type='l')
  
    ## Le résidus du modèle 
    plot(lca.male$residuals)
    
    
# Question4 :La projection centrale est une estimation moyenne sur 20 ans des taux de mortalite. Cette quantite peut donner une idee sur la projection totale.
    
      FR.stmomo.t<-StMoMoData(data=France1_ls_m ,series = "male",type="central")

      LC1 <- lc(link = "logit" )
      LCfit1 <- fit(LC1, data = central2initial(FR.stmomo.t), ages.fit = ages.fit,  years.fit = years.fit)
      
      horizon=20
      LCfor.t <- forecast(LCfit1, h = horizon)
      plot(LCfor.t)
      
      
      rates.t<-cbind(France1$rate$male[0:100,],LCfor.t$rate[0:100,])
      
      
      plot(seq(min(France1$year),max(France1$year)+20),rates.t[60,],xlab="Years",ylab="Death Rates",type="l",main="Taux observes et projetes a un horizon de 20 ans pour x = 60 ans")
      
      abline(v = 2010 , col="red" ,lwd=3, lty=2)
  #Question5
      
      chosen_cohort=2010        #doit appartenir aux years de LCfit
      
      plot(0:8, extractCohort(fitted(LCfit1, type = "rates"),
                              cohort = chosen_cohort),
           type = "l", log = "y", xlab = "age", ylab = "q(x)",
           main = paste(c("Cohort",toString(chosen_cohort),"mortality rates"), collapse = " "),
           xlim = c(0,103), ylim = c(0.000005, 0.007))
      
      #adding fitted projections
      #lines(9:28, extractCohort(LCfor.m$rates, cohort = chosen_cohort), 
           # lty = 2, lwd=2, col="red")
      
 #Question 6
      library(lifecontingencies)
      MaleFrance<-read.table(file="C:/Users/oussama/Downloads/Human-mortality-dashboard-main/Human-mortality-dashboard-main/mltper_1x1.txt", header = TRUE,skip = 1, sep = "", dec = ".")
      MaleFrance1<-MaleFrance[which(MaleFrance$Year == 2000),names(MaleFrance)]
      df_m<-data.frame(MaleFrance1)
      df_m$Age
      x_num <- as.numeric(df_m$Age)
      x_m=as.integer(x_num)
      lx_m<-df_m$lx
      
      df<-data.frame(x_m,lx_m)
      soa08Act=with(df_m, new("actuarialtable",x=x_m,lx=lx_m,name="France_male"))
      #evaluate and life-long annuity for an aged 65
      VAP=axn(soa08Act, x=60,n=30)
      VAP
      
      #Question 7
   
      
      
    
      
      non_smokers_percentage <- 0.3
      smokers_percentage <- 0.7
      
      total_number_of_assurees <- sum(lx_m)
      num_non_smokers <- round(non_smokers_percentage * total_number_of_assurees)
      num_smokers <- round(smokers_percentage * total_number_of_assurees)
      
      table_non_smokers <- new("actuarialtable", x = x_m, lx = lx_m * 0.9, name = "France_male_non_smokers")
      table_smokers <- new("actuarialtable", x = x_m, lx = lx_m * 1.15, name = "France_male_smokers")
      
      VAP_non_smokers <- axn(table_non_smokers, x = 55, n = 30)  
      VAP_smokers <- axn(table_smokers, x = 55, n = 30)  
      VAP_non_smokers
      VAP_smokers
      

      variation1 <- VAP_non_smokers - VAP
      variation1
      variation2 <- VAP_smokers - VAP
      variation2
      
      
