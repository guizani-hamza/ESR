########### IMPORTATION DES DONNEES ###############

### changer de repertoire ###

getwd()
#setwd("C:/Users/kateb/Documents/CNAM/DATA_ANALYST/S1/STA101/data") 
setwd("C:/Users/Fatiha Admin/Documents/data_Kateb/perso/CNAM/DATA_ANALYST/S1/STA101/data") 



library(FactoMineR) # n√©cessaire pour les analyses factorielles
library(factoextra) # permet d'extraire et de repr√©senter les sorties relatives √† l'analyse multidimensionnelle
library(missMDA) # n√©cessaire pour d√©terminer le nombre d'axes en AFDM par validation crois√©e
library(BioStatR) # n√©cessaire pour calculer le eta^2 entre deux variables
library(DescTools) # utile pour repr√©senter des matrices graphiquement (matrice de corr√©lation ou autres)
library(stargazer) # tr√®s utile pour les analyses univari√©es sous forme de tableaux  d'indicateurs (pour variable quantitatives)
library(questionr) # idem pour variables qualitatives (fonction freq notamment)
library(car)# notamment utilis√© pour construire des bo√Ætes √† moustaches telles que les individus en dehors des moustaches soient labellis√©s
library(xtable) # permet d'exporter des tableaux dans diff√©rents formats
library(missMDA)

################ LECTURE DONNEES ############ 

data_ESR <- read.csv("data_ESR_2021_uni_ecoles_short.csv",
                     header=TRUE, stringsAsFactors=TRUE, sep=";", na.strings="NA", dec=",", 
                     strip.white=TRUE, row.name=1, fileEncoding="latin1")

### visualiser les donnees ###

View(data_ESR)

### Fusionner les classes
str(data_ESR)
levels(data_ESR$TYPE_U)

#install.packages("tidyverse")
library(tidyverse)
data_ESR2 <- data_ESR %>% 
  mutate(TYPE_U = fct_collapse(TYPE_U,
                               USH=c("UTALLSHS", "UTDEG"), 
                               ECOLES=c("EI", "GE", "AUTRES")))

library(forcats)
fct_count(data_ESR2$TYPE_U)

### identifier les variables quali et quanti
quali <- which(sapply(data_ESR2, is.factor))
quanti <- which(sapply(data_ESR2, is.numeric))


########### EXPLORATION - ANALYSE UNIVARIEE ###############

summary(data_ESR2)

Summary_univ = summary(data_ESR2)
view(Summary_univ)
as.data.frame(Summary_univ)
write.csv(Summary_univ, file = 'summary.csv')

########## stargazer --> plus d'indicateurs que summary

# chargement de la librairie stargazer 

#install.packages("stargazer")   # telecharger la librairie 
library(stargazer)              # importer la librairie 

# affichage de quelques indicateurs statistiques pour variables quantitatives
summary_univ <- stargazer(data_ESR2[,quanti],
                          summary.stat = c("n","min","p25","median","mean","p75","max","sd"),
                          digits = 2,
                          decimal.mark = ".",
                          type = "text", title = "Tableau 1 : Analyse univariee des variables quantitatives " , out = "univariee.txt")

write(Summary_univ, file = 'summary.txt')

########### Representation graphique univariee
# Boite a moutache pour la variable Disque

cex = 0.5
Boxplot( ~ X.AM2D, data = data_ESR2, main = "Histogramme de la variable X.AM2D", cex.lab=1, cex.main=1, label.cex=0.2, outcex=0.8)

########## Histogrammes
library(lattice)
hist(data_ESR2$RESS_AUTRES,main="Histogramme de la variable RESS_autres",
     xlab="Autres ressources",
     ylab = "FrÈquence",
     col="darkgrey",
     freq=TRUE)

densite <- density(data_ESR2$X.AM2D) # estimer la densitÈ que reprÈsente ces diffÈrentes valeurs

lines(densite, col = "red",lwd=2) # Superposer une ligne de densitÈ ‡ l'histogramme

########### FIN - ANALYSE UNIVARIEE ######################

########### EXPLORATION - ANALYSE BIVARIEE ###############

##########VAR QUANTITATIVE##########

### 1- test de normalite : inspection visuelle + test shapiro

plot(data_ESR2[ , c(12,3:4,16)])

library("ggpubr")
# Diagramme de densitÈ
ggdensity(data_ESR2$SCSP_ETU, fill = "lightgray")
# QQ plot
ggqqplot(data_ESR2$SCSP_ETU, title ="q-q plot pour la varibale ETU"  )
shapiro.test(data_ESR2$SCSP_ETU)

### 2- Calcul des coefficients de correlation

mat.cor <- cor(data_ESR2[,quanti], use = "complete.obs")
mat.cor_spear <- cor(data_ESR2[,quanti], method = "spearman", use = "complete.obs")

write.csv(print(round(mat.cor, 2)), file = 'pearson.csv')
write.csv(print(round(mat.cor_spear, 2)), file = 'spearman.csv')


### 3- DÈtermination de la significativitÈ de la correlation

###### Pearson
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j],method = c("pearson") , ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
} 
p.mat <- cor.mtest(data_ESR2[,quanti])
write.csv(print(p.mat), file = "p_test_pearson.csv")

###### Spearman
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j],method = c("spearman") , ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
} 
p.mat_spear <- cor.mtest(data_ESR2[,quanti])
write.csv(print(p.mat_spear), file = "p_test_spearman.csv")

# representer matrice des correlations (heatmap)
#install.packages("corrplot")
library(corrplot)

# representation graphique de la matrice
corrplot(mat.cor_spear, method = "number" , title = "matrice des coefficients de Spearman", mar=c(2,0,2,0), type="upper", order="hclust",tl.col="black", tl.srt=45,  tl.cex = 0.5, number.cex = 0.5)  # cercles
corrplot(mat.cor, title = "matrice des coefficients de Pearson",  mar=c(2,0,1,0), type="upper", order="hclust",tl.col="black", tl.srt=45,  tl.cex = 0.5, number.cex = 0.5)  # cercles

# Avec Indication des correlations non significatives 
# nombres + test significativite (blanc)
corrplot(mat.cor, method="number", type="upper",  
         p.mat = p.mat, sig.level = 0.05, tl.col="black", tl.srt=45,  tl.cex = 0.5,number.cex = 0.6,insig = "blank", title = "Pearson", mar=c(2,0,1,0))
corrplot(mat.cor_spear, method="number", type="upper",  tl.col="black", tl.srt=45,  tl.cex = 0.5,number.cex = 0.5,
         p.mat = p.mat_spear, sig.level = 0.05, insig = "blank", title = "spearman", mar=c(2,0,1,0))  


##########VAR QUALITATIVE##########
# 1- test d'Ègalite des variances
var.test(data_ESR2$MAS_SAL ~ data_ESR2$TYPE_U)

identify_outliers(as.data.frame(data_ESR2$MAS_SAL))

#Determination des variances par groupe
library(rstatix)
data_ESR2 %>%
  group_by(data_ESR2$TYPE_U) %>%
  get_summary_stats(ETU, type = "mean_sd")

# Boxplot en fonction de TYPE_U
install.packages("ggplot2")
#remove.packages(ggplot2)
library(ggplot2)
ggplot(data_ESR2, aes(x=TYPE_U, y=ETU)) + 
  geom_boxplot(aes(fill = factor(TYPE_U)))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method=anova)

########### FIN - ANALYSE BIVARIEE ###############

########### EXPLORATION - ANALYSE MULTIVARIEE ###############

#Imputation +  ACP

# lignes de codes de factoshiny

nb <- missMDA::estim_ncpPCA(data_ESR2,quanti.sup=c(3,15),quali.sup=c(1,2,24))$ncp
dfcompleted <- missMDA::imputePCA(data_ESR2,ncp=nb,quanti.sup=c(3,15),quali.sup=c(1,2,24))$completeObs
res.PCA<-PCA(dfcompleted,quali.sup=c(1,2,24),quanti.sup=c(3,15),graph=FALSE)
summary(res.PCA)

### fin lignes de code factoshiny

write.csv(res.PCA$eig, file = "summary_PCA.csv")
dimdesc(res.PCA)

### diagramme eboulis

eigenvalues <- res.PCA$eig 
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Diagramme des Èboulis", 
        xlab = "Composantes principales", 
        ylab = "Pourcentages de variance", 
        col ="steelblue") 
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], type="b", pch=19, col = "red")


##### GRAPHIQUES : cercles correlation VAR + diagrammes individus

library("factoextra")

# graphe des var colore

fviz_pca_var(res.PCA, axes = c(1,2), col.var = "contrib",labelsize = 3, repel = TRUE,
             gradient.cols = c("white", "blue", "red"),
             ggtheme = theme_minimal())
fviz_pca_var(res.PCA, axes = c(3,4), col.var = "contrib",labelsize = 3, repel = TRUE,
             gradient.cols = c("white", "blue", "red"),
             ggtheme = theme_minimal())

# graphe des indiv colore

fviz_pca_ind(res.PCA, axes = c(1,2), col.ind = "contrib",
             labelsize = 3, repel = FALSE, gradient.cols = c("white", "blue", "red"),
             legend.title = data_ESR$TYPE_U, cex=0.1)
fviz_pca_ind(res.PCA, axes = c(3,4), col.ind = "contrib",gradient.cols = c("white", "blue", "red"),
             labelsize = 3, repel = FALSE,
             legend.title = data_ESR$TYPE_U, cex=0.1)


#indiv + ellipses

fviz_pca_ind(res.PCA, label="none", habillage=data_ESR2$TYPE_U,labelsize = 3, repel = TRUE,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")


#### cos≤ + contributions
dimdesc(res.PCA)
write.csv(res.PCA$ind$cos2, file = "PCA_ind.csv")

####################### FIN PCA ##################

#######################  CLASSIFICATION ##########

# CAH consolidee

res.HCPC<-HCPC(res.PCA,nb.clust=3,consol=TRUE,graph=FALSE)
plot.HCPC(res.HCPC,choice='tree',title='Arbre hiÈrarchique')
plot.HCPC(res.HCPC,choice='map',draw.tree=FALSE,title='Plan factoriel')
plot.HCPC(res.HCPC,choice='3D.map',ind.names=FALSE,centers.plot=FALSE,angle=60,title='Arbre hi?rarchique sur le plan factoriel')
res.HCPC$desc.var$test.chi2
res.HCPC$desc.var$category
res.HCPC$desc.axes
res.HCPC$desc.ind


fviz_dend(res.HCPC, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)


fviz_cluster(res.HCPC, axes = c(3,4),
             repel = FALSE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map", labelsize = 7, 
)


view(res.HCPC$data.clust)
table(res.HCPC$data.clust[,24])
table(data_ESR2$TYPE_U, (res.HCPC$data.clust[,24]))
write.csv(res.HCPC$data.clust, file = "classe.csv")

################ FIN CLASSIFICATION ################
