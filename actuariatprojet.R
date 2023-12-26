x<-read.table(file=file.choose(),header = TRUE,blank.lines.skip = TRUE,sep="," )
library(survival)
install.packages("StMoMo")
library(StMoMo)
install.packages("rlang")


# Création du modèle de survie (modèle de Lee-Carter)
model <- StMoMo(formula = ~ Year + Age, data = x, method = "LC")

# Estimation des paramètres du modèle
fit_model <- fit(model)

# Calcul des taux de mortalité estimés
mortality_rates <- rates(fit_model)

# Calcul des intervalles de confiance au seuil de 99% (méthode de Wald)
lower_ci <- mortality_rates$lowerCI
upper_ci <- mortality_rates$upperCI

# Tracé des taux de mortalité et des intervalles de confiance
plot(x$Age, mortality_rates$rates, type = "l", xlab = "Âge", ylab = "Taux de mortalité",
     main = "Estimation des taux de mortalité par âge")
lines(x$Age, lower_ci, lty = 2, col = "red")
lines(x$Age, upper_ci, lty = 2, col = "red")
legend("topright", legend = c("Taux de mortalité", "Intervalle de confiance"),
       lty = c(1, 2), col = c("black", "red"))

data2 <- new("lifetable", x = as.vector(x$Age), lx = as.vector(x$Lx))
head(data2)






