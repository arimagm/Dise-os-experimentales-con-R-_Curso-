# Curso-Taller en Linea: Diseño experimentales con R
# 20 de noviembre del 2023
# Flaviano Godínez Jaimes

library(multcomp)     # Multiple comparisons of means in AN(C)OVA: Dunett
library(agricolae)    #SNK, DHS(Tukey), Duncan, LSD
library(car)          # Diagnosticos y leveneTest
library(nortest)      #Pruebas de normalidad
library(tseries)      # Jarque-Bera

#--------------------------------------------------------------------------#
#- Los datos -#
#--------------------------------------------------------------------------#
setwd("/Users/flaviano/Documents/UAG 2012/Cursos/Diseño de experimentos/DE 2014/DE Wily/Tareas con tex y scrip de R")
datos<-read.csv("Ejercicio1 de DBA.csv")  
datos
summary(datos)
names(datos)  # Nombres de las variables
attach(datos)

#________________________________________________________________________#
#                        Analisis exploratorio                           #
#________________________________________________________________________#

MDureza=mean(Dureza);  round(MDureza, 2);     VDureza=var(Dureza);   round(VDureza, 2)
MPunta = tapply(Dureza, Punta, mean); round(MPunta, 2)
VPunta = tapply(Dureza, Punta, var); round(VPunta, 5)
MPlaca = tapply(Dureza, Placa, mean); round(MPlaca, 2)
VPlaca = tapply(Dureza, Placa, var);  round(VPlaca, 2)

colBlq <- c("tomato", "purple","pink","gold")
colPta <- c("cyan", "royalblue","steelblue","blue")
boxplot(Dureza ~ Placa, col = colBlq,  main="Dureza")
boxplot(Dureza ~ Punta, col = colPta,  main="Dureza")

interaction.plot(Punta, Placa, Dureza, lty=1, lwd=2, col =colBlq ,  main="")
interaction.plot(Placa, Punta, Dureza, lty=1, lwd=2, col = colPta, main="")

png("FDBCAPla.png")
boxplot(Dureza ~ Placa, col = colBlq,  main=" ")
dev.off()

png("FDBCAPta.png")
boxplot(Dureza ~ Punta, col = colPta,  main=" ")
dev.off()

png("FDBCAInt1.png")
interaction.plot(Punta, Placa, Dureza, lty=1, lwd=2, col =colBlq ,  main="")
dev.off()

png("FDBCAInt2.png")
interaction.plot(Placa, Punta, Dureza, lty=1, lwd=2,  col = colPta, main="")
dev.off()



#---------------------------------------------------------------------#
# Aleatorización 
#---------------------------------------------------------------------#

# 4 Tratamientos y 4 bloques
trt<-c("P1","P2","P3","P4")
outdesign <-design.rcbd(trt, 4, serie=1, seed = 986, kinds = "Wichmann-Hill") # seed = 986
outdesign


#________________________________________________________________________#
#                         Ajuste del modelo BA                           #
#________________________________________________________________________#

g <- lm(Dureza ~ Placa+Punta, data=datos)
anova(g);            summary(g)

#________________________________________________________________________#
#            Verificacion de los supuestos del modelo BA                           #
#________________________________________________________________________#
layout(matrix(1:4, 2,2))
plot(g, main="")
dev.off()

plot(g, main="", which=1)
plot(g, main="", which=2)
plot(g, main="", which=3)
plot(g, main="", which=4)

qqnorm(g$res,xlab="Cuantiles teoricos",ylab="Residuos estandarizados",main="")
qqline(g$res, col="blue")

png("FDBCANorm.png")
qqnorm(g$res,xlab="Cuantiles teoricos",ylab="Residuos estandarizados",main="")
qqline(g$res, col="blue")
dev.off()

car::qqPlot(g, id.n=3)         

png("FDBCANormT.png")
car::qqPlot(g, id.n=3)        
dev.off()


h1<-hist(g$res,col="gold",xlim=0.3*c(-1,1))
normal.freq(h1,col="blue")

# Normalidad
# H0: Los errores tienen DN vs 
# H1: Los errores NO tienen DN

  stats::shapiro.test(g$res)        # n<50
nortest::lillie.test(g$res)        # n>=50    
nortest::ad.test(g$res)
nortest::cvm.test(g$res)
nortest::pearson.test(g$res)    
nortest::sf.test(g$res)       
tseries::jarque.bera.test(g$res)    #
# Todos los valores p son mayores a 0.05

#  Pruebas de igualdad de varianzas
plot(g$fit,g$res,xlab="Ajustados",ylab="Residuos",main="")


png("FDBCAHom.png")
plot(g$fit,g$res,xlab="Ajustados",ylab="Residuos",main="")
dev.off()


# Pruebas de hipótesis. Homegeneidad de varianzas
# H0:Las varianzas de los TC de los tratamientos son iguales vs
# H1:Las varianzas de los TC de los tratamientos NO son iguales
stats::bartlett.test(g$res, Punta)   # Confiable si hay Normalidad
  car::leveneTest(g$res, Punta)
stats::fligner.test(g$res, Punta)
# Todos los valores p son mayores a 0.05

#________________________________________________________________________#
#  Todo bien!
# Podemos interpretar los resultados!
anova(g);            summary(g)


modelo <- aov(Dureza ~ Placa+Punta)
model.tables(modelo,"mean")

#________________________________________________________________________#
#                       Comparaciones multiples                          #
#________________________________________________________________________#


gaov<-aov(Dureza ~ Placa+Punta, data=datos)

#Prueba de Tukey: Honest Significative Difference
CompHSD=HSD.test(gaov, "Punta",  main="")
CompHSD


#Prueba de Duncan
CompDun=duncan.test(gaov, "Punta", main="")
CompDun

#Prueba LSD
CompLSD=LSD.test(gaov, "Punta", p.adj="bonferroni", main="")
CompLSD

#Prueba de Student-Newman-Keuls
CompSNK=SNK.test(gaov, "Punta",  main="")
CompSNK

#Prueba de Scheffe
Compscheffe=scheffe.test(gaov, "Punta",  main="")
Compscheffe

# .......................................................
#Prueba de Dunnett
summary(glht(Fgaov, linfct=mcp(FPunta="Dunnett")))


detach()

