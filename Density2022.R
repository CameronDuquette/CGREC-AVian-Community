maxcomb=read.csv("Avcom.csv")



####hullyuh includes all species. use for diversity stuff
ordbirds2=maxcomb[2:66]
birdmax2=metaMDS(ordbirds2, distance="bray", k=3, trymax=5000)


#gg_ordiplot(birdmax2,birdenv$TRT.y, kind="sd" )+theme_bw()


plot(birdmax2)

birdenv=maxcomb[,c(67:ncol(maxcomb))]
View(birdenv)
colnames(birdenv)

env4test=birdenv[,c(3,4,5,11,13,15,16,17,18,19,23,24,51,53)]
envtest2=env4test %>% 
  rename(
    KBG = meanpopr.x,
    BROME = meanbrin.x,
    COOL = meanNC3.x,
    FORB = meanNForb.x,
    WOOD = meanNWood.x,
    BARE = meanBare.x,
    DEAD = meanDead.x,
    ROBEL = R.x,
    LITDEP = ld.x,
    TRT = TRT.x,
    SOIL = SiteName.x,
    ROUGH = toposd15_var_2.x,
    WET = Shape.Area,
    PASTYR=pasyr
    
  )


vf <- envfit(birdmax2, envtest2[,c(10:13)], perm = 999)

site.scrs <- as.data.frame(scores(birdmax2, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, TRT = envtest2$TRT) #add grouping variable "Management" to dataframe

env.scores.CAD <- as.data.frame(scores(vf, display = "vectors")) #extracts relevant scores from envifit
env.scores.CAD <- cbind(env.scores.CAD, env.variables = rownames(env.scores.CAD)) #and then gives them their names

env.scores.CAD <- cbind(env.scores.CAD, pval = vf$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.CAD, pval<=0.05) #subset data to show variables significant at 0.05

head(sig.env.scrs)

species.scores <- as.data.frame(scores(birdmax2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data




site.scrs$TRT=as.factor(site.scrs$TRT)
##########lets get plotting
nmds.plot.dune <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  #geom_point(aes(NMDS1, NMDS2, colour = factor(TRT), size = 2))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "TRT")# add legend labels for Management and Landuse
  #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot.dune
##########

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the management factor
df_ell.dune.management <- data.frame() #sets up a data frame before running the function.
for(g in levels(site.scrs$TRT)){
  df_ell.dune.management <- rbind(df_ell.dune.management, cbind(as.data.frame(with(site.scrs [site.scrs$TRT==g,],
                                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,TRT=g))
}

# data for labelling the ellipse
NMDS.mean.dune=aggregate(site.scrs[ ,c("NMDS1", "NMDS2")], 
                         list(group = site.scrs$TRT), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(site.scrs[,c("NMDS1", "NMDS2")], 
                    list(group = site.scrs$TRT), mean)
nmds.plot.dune+ 
  geom_path(data = df_ell.dune.management, aes(x = NMDS1, y = NMDS2, group = TRT, colour=factor(TRT, labels=c("PBG20", "PBG40", "SLG", "MTRG")), size=.5))+#this is the ellipse, seperate ones by Site. 
  geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black",alpha=0.5, max.overlaps = 15)+theme_classic()

######create spider plot
  

birdscrs20<-scores(birdmax2, display='sites')
birdscrs20<-cbind(as.data.frame(birdscrs20), Treatment=maxcomb$TRT.y)

birdspecies.scores20 <- as.data.frame(scores(birdmax2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
birdspecies.scores20$species <- rownames(birdspecies.scores20)  # create a column of species, from the rownames of species.scores

birdcent20=aggregate(cbind(NMDS1, NMDS2)~Treatment, data=birdscrs20, FUN=mean)

birdsegs20 <- merge(birdscrs20, setNames(birdcent20, c('Treatment','oNMDS1','oNMDS2')),
                by = 'Treatment', sort = FALSE)

birdsegs20$lineColor <- ifelse(birdsegs20$Treatment=="PBG20", "red",
                           ifelse(birdsegs20$Treatment=="PBG40", "green",
                                  ifelse(birdsegs20$Treatment=="SLG", "blue","purple")))
birdscrs20$lineColor <- ifelse(birdscrs20$Treatment=="PBG20", "red",
                           ifelse(birdscrs20$Treatment=="PBG40", "green",
                                  ifelse(birdscrs20$Treatment=="SLG", "blue","purple")))

birdspecies.scores20.2=birdspecies.scores20[c(12,17,18,29,51,60),]

##this is vegetation plot. change to birds
birdvgg20=ggplot(birdscrs20, aes(x = NMDS1, y = NMDS2, color=Treatment)) +
  geom_segment(data = birdsegs20,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = birdcent20, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.title.y = element_text(face='bold', size=12),axis.title.x = element_text(face='bold', size=12))+
  geom_text_repel(data=birdspecies.scores20.2,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15)+
  scale_color_manual(values = c("#00BC59","#00A5FF", "#E7B800", "#FC4E07"))+
  theme(axis.text.x= element_text(colour="black",size=11))+
  theme(axis.text.y= element_text(colour="black",size=11))


##end of spider plot business






####ordination with inherent het

maxcomb=read.csv("Avcom.csv")




####hullyuh includes all species. use for diversity stuff
ordbirds2=maxcomb[2:66]
birdmax2=metaMDS(ordbirds2, distance="bray", k=3, trymax=5000)

plot(birdmax2)

birdenv=maxcomb[,c(67:ncol(maxcomb))]
View(birdenv)
colnames(birdenv)

env4test=birdenv[,c(3,4,5,11,13,15,16,17,18,19,23,24,51)]
envtest2=env4test %>% 
  rename(
    KBG = meanpopr.x,
    BROME = meanbrin.x,
    COOL = meanNC3.x,
    FORB = meanNForb.x,
    WOOD = meanNWood.x,
    BARE = meanBare.x,
    DEAD = meanDead.x,
    ROBEL = R.x,
    LITDEP = ld.x,
    TRT = TRT.x,
    SOIL = SiteName.x,
    ROUGHNESS = toposd15_var_2.x,
    WETLAND = Shape.Area
    
  )


vf <- envfit(birdmax2, envtest2[,c(10:13)], perm = 999)

site.scrs <- as.data.frame(scores(birdmax2, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, SOIL = envtest2$SOIL) #add grouping variable "Management" to dataframe

env.scores.CAD <- as.data.frame(scores(vf, display = "vectors")) #extracts relevant scores from envifit
env.scores.CAD <- cbind(env.scores.CAD, env.variables = rownames(env.scores.CAD)) #and then gives them their names

env.scores.CAD <- cbind(env.scores.CAD, pval = vf$vectors$pvals) # add pvalues to dataframe

head(sig.env.scrs)

species.scores <- as.data.frame(scores(birdmax2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data




site.scrs$SOIL=as.factor(site.scrs$SOIL)
##########lets get plotting
######create spider plot


birdscrs20<-scores(birdmax2, display='sites')
birdscrs20<-cbind(as.data.frame(birdscrs20), Soil=maxcomb$SiteName.y)

birdspecies.scores20 <- as.data.frame(scores(birdmax2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
birdspecies.scores20$species <- rownames(birdspecies.scores20)  # create a column of species, from the rownames of species.scores

birdcent20=aggregate(cbind(NMDS1, NMDS2)~Soil, data=birdscrs20, FUN=mean)

birdsegs20 <- merge(birdscrs20, setNames(birdcent20, c('Soil','oNMDS1','oNMDS2')),
                    by = 'Soil', sort = FALSE)

birdspecies.scores20.2=birdspecies.scores20[c(12,17,18,29,51,60),]

##this is vegetation plot. change to birds
birdsoilvgg20=ggplot(birdscrs20, aes(x = NMDS1, y = NMDS2)) +
 geom_segment(data = birdsegs20,
               mapping = aes(xend = oNMDS1, yend = oNMDS2, color=Soil)) + # spiders
  geom_point(data = birdcent20, size = 5, mapping = aes(color=Soil)) +                         # centroids
  geom_point(mapping = aes(color=Soil))+                                         # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.title.y = element_text(face='bold', size=12),axis.title.x = element_text(face='bold', size=12))+
  geom_text_repel(data=birdspecies.scores20.2,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15, size=3)+
  #scale_color_manual(values = c("#00BC59","#00A5FF", "#E7B800", "#FC4E07"))
  geom_segment(data = env.scores.CAD, aes(x = 0, xend = 4*NMDS1, y = 0, yend = 4*NMDS2),
                            colour = "black", size=1.0) +
               geom_text(data = env.scores.CAD, aes(x = 7*NMDS1, y = 7*NMDS2, label = env.variables, fontface='bold'),size = 3, vjust=0.0,hjust=-.02, colour="black")+
  theme(axis.text.x= element_text(colour="black",size=11))+
  theme(axis.text.y= element_text(colour="black",size=11))




#birdsoilenv=ggplot()+geom_segment(data = env.scores.CAD, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#  geom_text(data = env.scores.CAD, aes(x = NMDS1, y = NMDS2, label = env.variables),size = 3)

####end of ordination with inherent het

#betadisperson test
head(envtest2)

dis <- vegdist(ordbirds2,method="bray")
mod <- betadisper(dis, envtest2$PASTYR, bias.adjust = TRUE)
disavg=as.data.frame(mod$group)
disavg$dist=mod$distances

###########STARTUPHERE4_11###disavg is each row's distance to pasture centroid
disavg=read.csv("vegdstden22.csv")
mod=aov(dist~TRT, data=disavg)
(mod.HSD <- TukeyHSD(mod))

adonis.II(ordbirds2 ~ SiteName.x+TRT.x+toposd15_var_2.x+Shape.Area+SiteName.x*TRT.x+toposd15_var_2.x*TRT.x+Shape.Area*TRT.x, data = birdenv, by = NULL)

pairwise.adonis(ordbirds2, factors=birdenv$SiteName.x)

###disavg is each row's distance to pasture centroid
disavg=read.csv("distancetopasture.csv")
mod=aov(dist~TRT, data=disavg)
(mod.HSD <- TukeyHSD(mod))
#####bobodensitysetup

birdgnf=read.csv("Communityalldone.csv")
View(birdgnf)
###C&E for wag transects
dissertveg=read.csv("Final2020vegagg.csv")
###C&E for wag transects

#remove MTORG for relative trt analysis
transectlist=read.csv("transectlist.csv")
###C&E for wag transects
bobogf<-subset(birdgnf,Species=="BOBO") #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES

###merge subset with transect list to add in empty transects

bobogf<-merge(bobogf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
bobogf$ptyr=gsub('.$', '', bobogf$bloop)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

##merge new subset with veg so we can line them up for analysis with sort
bobogf<-merge(bobogf, dissertveg, by.x='ptyr', by.y='ptyr')
View(bobogf)
bobogf$SiteName=factor(bobogf$SiteName.x)

boboveg=bobogf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

bobogfform<-formatDistData(distData = bobogf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
bobogfform <- bobogfform[order(row.names(bobogfform)),]

boboveg=merge(bobogfform, dissertveg, by.x=bobogfform[1,], by.y='ptyr')

boboveg <- boboveg[order(row.names(boboveg)),]

#add levels to categorical covs

boboveg$TRT
boboveg$TRT=factor(boboveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

boboveg$toposcale=scale(boboveg$toposd15_var_2)
boboveg$wetscale= scale(boboveg$Shape.Area)

### make sure to merge after the distance matrix is made
gfbobo<-unmarkedFrameGDS(bobogfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=boboveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
bogfdist1<-gdistsamp(~1,~1,~1,gfbobo,keyfun = "halfnorm",output = "density")
bogfdist2<-gdistsamp(~1,~1,~1,gfbobo,keyfun = "exp",output = "density")
bogfdist3<-gdistsamp(~1,~1,~1,gfbobo,keyfun = "hazard",output = "density")
bogfdist4<-gdistsamp(~1,~1,~1,gfbobo,keyfun = "uniform",output = "density")
sfit<-fitList(bogfdist1,bogfdist2,bogfdist3,bogfdist4) #Put model results in a list ordered by AIC
modSel(sfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

gfbobo@siteCovs

##exp is best

bogfdist1<-gdistsamp(~1,~1,~1,gfbobo,keyfun = "exp",output = "density")
bogfsoiltype<-gdistsamp(~SiteName,~1,~1,gfbobo,keyfun = "exp",output = "density")
bogftopovar<-gdistsamp(~toposcale,~1,~1,gfbobo,keyfun = "exp",output = "density")
bogfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfbobo,keyfun = "exp",output = "density")
bogfTRT<-gdistsamp(~TRT,~1,~1,gfbobo,keyfun = "exp",output = "density")

bobofit<-fitList(bogfdist1, bogfsoiltype, bogftopovar, bogfsoiltypetopo, bogfTRT)
modSel(bobofit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES
bogfsoiltype

bobopre<-predict(bogfsoiltype, type="lambda", appendData=T, SiteName=Sandy)

bobosoil<- data.frame(SiteName=factor(c("Loamy", "Sandy", "Shallow Gravel", "Thin Loamy", "Very Shallow")))

bobotrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(Elambdatrt <- predict(bogfTRT, type="lambda", newdata=bobotrt,
                    appendData=TRUE))

(Elambda <- predict(bogfsoiltype, type="lambda", newdata=bobosoil,
                      appendData=TRUE))
View(bobopre)

#################################
cclogf<-subset(birdgnf,Species=="CCLO") #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES
#View(bobogf)
#View(transectlist2)
#### Create plot
ggplot(cclogf,aes(x=Distance))+geom_histogram(aes(y=..density..),binwidth = 10) 
###merge subset with transect list to add in empty transects

cclogf<-merge(cclogf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
cclogf$ptyr=gsub('.$', '', cclogf$bloop)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

##merge new subset with veg so we can line them up for analysis with sort
cclogf<-merge(cclogf, dissertveg, by.x='ptyr', by.y='ptyr')

cclogf$SiteName=factor(cclogf$SiteName.x)

ccloveg=cclogf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

cclogfform<-formatDistData(distData = cclogf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
cclogfform <- cclogfform[order(row.names(cclogfform)),]

ccloveg=merge(cclogfform, dissertveg, by.x=cclogfform[1,], by.y='ptyr')

ccloveg <- ccloveg[order(row.names(ccloveg)),]

#add levels to categorical covs

ccloveg$TRT=factor(ccloveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

ccloveg$toposcale=scale(ccloveg$toposd15_var_2)
ccloveg$wetscale= scale(ccloveg$Shape.Area)

### make sure to merge after the distance matrix is made
gfcclo<-unmarkedFrameGDS(cclogfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=ccloveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
cclogfdist1<-gdistsamp(~1,~1,~1,gfcclo,keyfun = "halfnorm",output = "density")
cclogfdist2<-gdistsamp(~1,~1,~1,gfcclo,keyfun = "exp",output = "density")
cclogfdist3<-gdistsamp(~1,~1,~1,gfcclo,keyfun = "hazard",output = "density")
cclogfdist4<-gdistsamp(~1,~1,~1,gfcclo,keyfun = "uniform",output = "density")
cclosfit<-fitList(cclogfdist1,cclogfdist2,cclogfdist3,cclogfdist4) #Put model results in a list ordered by AIC
modSel(cclosfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

##exp is best

cclogfdist1<-gdistsamp(~1,~1,~1,gfcclo,keyfun = "hazard",output = "density")
cclogfsoiltype<-gdistsamp(~SiteName,~1,~1,gfcclo,keyfun = "hazard",output = "density")
cclogftopo<-gdistsamp(~toposcale,~1,~1,gfcclo,keyfun = "hazard",output = "density")
cclogfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfcclo,keyfun = "hazard",output = "density")
cclogfTRT<-gdistsamp(~TRT,~1,~1,gfcclo,keyfun = "hazard",output = "density")

cclofit2<-fitList(cclogfdist1, cclogfsoiltype, cclogftopo, cclogfsoiltypetopo, cclogfTRT)
modSel(cclofit2) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES


cclopre<-predict(cclogfdist1, type="lambda", appendData=T)
View(cclopre)


cclotrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(Elambdatrt <- predict(cclogfTRT, type="lambda", newdata=cclotrt,
                       appendData=TRUE))

(Elambda <- predict(bogfsoiltype, type="lambda", newdata=bobosoil,
                    appendData=TRUE))

#######################################################CCSP

ccspgf<-subset(birdgnf,Species=="CCSP") #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES
#### Create plot
ggplot(ccspgf,aes(x=Distance))+geom_histogram(aes(y=..density..),binwidth = 10) #Plots Sora distances in a histogram with density on y axis

###merge subset with transect list to add in empty transects

ccspgf<-merge(ccspgf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
ccspgf$ptyr=gsub('.$', '', ccspgf$bloop)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")


##merge new subset with veg so we can line them up for analysis with sort
ccspgf<-merge(ccspgf, dissertveg, by.x='ptyr', by.y='ptyr')

ccspgf$SiteName=factor(ccspgf$SiteName.x)

ccspveg=ccspgf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

ccspgfform<-formatDistData(distData = ccspgf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
ccspgfform <- ccspgfform[order(row.names(ccspgfform)),]

ccspgfform2=as.data.frame(ccspgfform)
ccspgfform2$ptyr=rownames(ccspgfform2)

ccspveg=merge(ccspgfform2, dissertveg, by.x='ptyr', by.y='ptyr')

ccspveg <- ccspveg[order(ccspveg$ptyr),]

#add levels to categorical covs

ccspveg$TRT=factor(ccspveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

ccspveg$toposcale=scale(ccspveg$toposd15_var_2)
ccspveg$wetscale= scale(ccspveg$Shape.Area)
### make sure to merge after the distance matrix is made
gfccsp<-unmarkedFrameGDS(ccspgfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=ccspveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
ccspgfdist1<-gdistsamp(~1,~1,~1,gfccsp,keyfun = "halfnorm",output = "density")
ccspgfdist2<-gdistsamp(~1,~1,~1,gfccsp,keyfun = "exp",output = "density")
ccspgfdist3<-gdistsamp(~1,~1,~1,gfccsp,keyfun = "hazard",output = "density")
ccspgfdist4<-gdistsamp(~1,~1,~1,gfccsp,keyfun = "uniform",output = "density")
sfit<-fitList(ccspgfdist1,ccspgfdist2,ccspgfdist3,ccspgfdist4) #Put model results in a list ordered by AIC
modSel(sfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

gfbobo@siteCovs

##exp is best

ccspgfdist1<-gdistsamp(~1,~1,~1,gfccsp,keyfun = "exp",output = "density")
ccspgfsoiltype<-gdistsamp(~SiteName,~1,~1,gfccsp,keyfun = "exp",output = "density")
ccspgftopovar<-gdistsamp(~toposcale,~1,~1,gfccsp,keyfun = "exp",output = "density")
ccspgfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfccsp,keyfun = "exp",output = "density")
ccspgfTRT<-gdistsamp(~TRT,~1,~1,gfccsp,keyfun = "exp",output = "density")

ccspfit<-fitList(ccspgfdist1, ccspgfsoiltype, ccspgftopovar, ccspgfsoiltypetopo, ccspgfTRT)
modSel(ccspfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

ccsppre<-predict(ccspgfdist1, type="lambda", appendData=T)
View(ccsppre)

ccsptrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(Elambdatrtccsp <- predict(ccspgfTRT, type="lambda", newdata=ccsptrt,
                       appendData=TRUE))

(Elambda <- predict(bogfsoiltype, type="lambda", newdata=bobosoil,
                    appendData=TRUE))

#################################################SAVS
savsgf<-subset(birdgnf,Species=="SAVS"& Distance <55) #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES
#### Create plot
ggplot(savsgf,aes(x=Distance))+geom_histogram(aes(y=..density..),binwidth = 10) 

###merge subset with transect list to add in empty transects

savsgf<-merge(savsgf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
savsgf$ptyr=gsub('.$', '', savsgf$bloop)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

##merge new subset with veg so we can line them up for analysis with sort
savsgf<-merge(savsgf, dissertveg, by.x='ptyr', by.y='ptyr')

savsgf$SiteName=factor(savsgf$SiteName.x)

savsveg=savsgf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

savsgfform<-formatDistData(distData = savsgf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
savsgfform <- savsgfform[order(row.names(savsgfform)),]

savsgfform2=as.data.frame(savsgfform)
savsgfform2$ptyr=rownames(savsgfform2)

savsveg=merge(savsgfform2, dissertveg, by.x='ptyr', by.y='ptyr')

savsveg <- savsveg[order(savsveg$ptyr),]

#add levels to categorical covs

savsveg$TRT=factor(savsveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

savsveg$toposcale=scale(savsveg$toposd15_var_2)
savsveg$wetscale= scale(savsveg$Shape.Area)
### make sure to merge after the distance matrix is made
gfsavs<-unmarkedFrameGDS(savsgfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=savsveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
savsgfdist1<-gdistsamp(~1,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")
savsgfdist2<-gdistsamp(~1,~1,~1,gfsavs,keyfun = "exp",output = "density")
savsgfdist3<-gdistsamp(~1,~1,~1,gfsavs,keyfun = "hazard",output = "density")
savsgfdist4<-gdistsamp(~1,~1,~1,gfsavs,keyfun = "uniform",output = "density")
sfit<-fitList(savsgfdist1,savsgfdist2,savsgfdist3,savsgfdist4) #Put model results in a list ordered by AIC
modSel(sfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

gfbobo@siteCovs

##exp is best

savsgfdist1<-gdistsamp(~1,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")
savsgfsoiltype<-gdistsamp(~SiteName,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")
savsgftopovar<-gdistsamp(~toposcale,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")
savsgfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")
savsgfTRT<-gdistsamp(~TRT,~1,~1,gfsavs,keyfun = "halfnorm",output = "density")

savsfit<-fitList(savsgfdist1, savsgfsoiltype, savsgftopovar, savsgfsoiltypetopo, savsgfTRT)
modSel(savsfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

savspre<-predict(savsgfdist1, type="lambda", appendData=T)
View(savspre)

savssoil<- data.frame(SiteName=factor(c("Loamy", "Sandy", "Shallow Gravel", "Thin Loamy", "Very Shallow")))

savstrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(savstrt <- predict(savsgfTRT, type="lambda", newdata=savstrt,
                       appendData=TRUE))

(savssoil <- predict(savsgfsoiltype, type="lambda", newdata=savssoil,
                    appendData=TRUE))

##############WEME
wemegf<-subset(birdgnf,Species=="weme"& Distance <55) #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES
#### Create plot
ggplot(wemegf,aes(x=Distance))+geom_histogram(aes(y=..density..),binwidth = 10) 

###merge subset with transect list to add in empty transects

wemegf<-merge(wemegf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
wemegf$ptyr=gsub('.$', '', wemegf$bloop)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

##merge new subset with veg so we can line them up for analysis with sort
wemegf<-merge(wemegf, dissertveg, by.x='ptyr', by.y='ptyr')

wemegf$SiteName=factor(wemegf$SiteName.x)

wemeveg=wemegf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

wemegfform<-formatDistData(distData = wemegf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
wemegfform <- wemegfform[order(row.names(wemegfform)),]

wemegfform2=as.data.frame(wemegfform)
wemegfform2$ptyr=rownames(wemegfform2)

wemeveg=merge(wemegfform2, dissertveg, by.x='ptyr', by.y='ptyr')

wemeveg <- wemeveg[order(wemeveg$ptyr),]

#View(boboveg)
#add levels to categorical covs

wemeveg$TRT=factor(wemeveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

wemeveg$toposcale=scale(wemeveg$toposd15_var_2)
wemeveg$wetscale= scale(wemeveg$Shape.Area)

### make sure to merge after the distance matrix is made
gfweme<-unmarkedFrameGDS(wemegfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=wemeveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
wemegfdist1<-gdistsamp(~1,~1,~1,gfweme,keyfun = "halfnorm",output = "density")
wemegfdist2<-gdistsamp(~1,~1,~1,gfweme,keyfun = "exp",output = "density")
wemegfdist3<-gdistsamp(~1,~1,~1,gfweme,keyfun = "hazard",output = "density")
wemegfdist4<-gdistsamp(~1,~1,~1,gfweme,keyfun = "uniform",output = "density")
sfit<-fitList(wemegfdist1,wemegfdist2,wemegfdist3,wemegfdist4) #Put model results in a list ordered by AIC
modSel(sfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

gfbobo@siteCovs

##exp is best

wemegfdist1<-gdistsamp(~1,~1,~1,gfweme,keyfun = "halfnorm",output = "density")
wemegfsoiltype<-gdistsamp(~SiteName,~1,~1,gfweme,keyfun = "halfnorm",output = "density")
wemegftopovar<-gdistsamp(~toposcale,~1,~1,gfweme,keyfun = "halfnorm",output = "density")
wemegfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfweme,keyfun = "halfnorm",output = "density")
wemegfTRT<-gdistsamp(~TRT,~1,~1,gfweme,keyfun = "halfnorm",output = "density")

wemefit<-fitList(wemegfdist1, wemegfsoiltype, wemegftopovar, wemegfsoiltypetopo, wemegfTRT)
modSel(wemefit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

wemepre<-predict(wemegfdist1, type="lambda", appendData=T)
View(wemepre)

savssoil<- data.frame(SiteName=factor(c("Loamy", "Sandy", "Shallow Gravel", "Thin Loamy", "Very Shallow")))

wemetrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(wemeltrt <- predict(wemegfTRT, type="lambda", newdata=wemetrt,
                    appendData=TRUE))

(savssoil <- predict(savsgfsoiltype, type="lambda", newdata=savssoil,
                     appendData=TRUE))

##############################GRSP
grspgf<-subset(birdgnf,Species=="GRSP"& Distance <55) #THIS WILL NEED TO BE CHANGED FOR EACH SPECIES
#### Create plot
ggplot(grspgf,aes(x=Distance))+geom_histogram(aes(y=..density..),binwidth = 10) #Plots Sora distances in a histogram with density on y axis

###merge subset with transect list to add in empty transects

grspgf<-merge(grspgf, transectlist, by.x="bloop", by.y="bloop1", all=TRUE)
grspgf$ptyr=gsub('.$', '', grspgf$bloop)

#View(bobogf)

dissertveg$ptyr=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

##merge new subset with veg so we can line them up for analysis with sort
grspgf<-merge(grspgf, dissertveg, by.x='ptyr', by.y='ptyr')

grspgf$SiteName=factor(grspgf$SiteName.x)

grspveg=grspgf[,c(1, 22:54)] #species specific veg covs

##format distance data
library(unmarked)

grspgfform<-formatDistData(distData = grspgf,distCol = "Distance",transectNameCol = "ptyr",
                           dist.breaks = c(0,10,20,30,40,50), occasionCol ="Round")
grspgfform <- grspgfform[order(row.names(grspgfform)),]

grspgfform2=as.data.frame(grspgfform)
grspgfform2$ptyr=rownames(grspgfform2)

grspveg=merge(grspgfform2, dissertveg, by.x='ptyr', by.y='ptyr')

grspveg <- grspveg[order(grspveg$ptyr),]

#View(boboveg)

#add levels to categorical covs

grspveg$TRT=factor(grspveg$TRT,levels=c("PBG20","PBG40", "SLG", "WAG"),ordered=FALSE)

View(scaledgrsp)
###scale grsp veg

grspveg$toposcale=scale(grspveg$toposd15_var_2)
grspveg$wetscale= scale(grspveg$Shape.Area)


### make sure to merge after the distance matrix is made
gfgrsp<-unmarkedFrameGDS(grspgfform,dist.breaks = c(0,10,20,30,40,50), siteCovs=grspveg, tlength =rep(150,480),
                         survey = "line",unitsIn = "m",numPrimary = 4) #NOTE HERE 'numPrimary' IS ACTUALLY NUMBER OF SECONDARY SAMPLING OCASSIONS (VISITS)
#### Run constant (i.e., null) models to determine best key function to use in subsequent models
#### ~ REPRESENTS THE EQUATION YOU WANT TO MODEL ON THE PARAMETER
#### FIRST PARAMETER IS ABUNDANCE/DENSITY, SECOND IS AVAILABILITY, THIRD IS DETECTION
grspgfdist1<-gdistsamp(~1,~1,~1,gfgrsp,keyfun = "halfnorm",output = "density")
grspgfdist2<-gdistsamp(~1,~1,~1,gfgrsp,keyfun = "exp",output = "density")
grspgfdist3<-gdistsamp(~1,~1,~1,gfgrsp,keyfun = "hazard",output = "density")
grspgfdist4<-gdistsamp(~1,~1,~1,gfgrsp,keyfun = "uniform",output = "density")
sfit<-fitList(grspgfdist1,grspgfdist2,grspgfdist3,grspgfdist4) #Put model results in a list ordered by AIC
modSel(sfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

##exp is best

grspgfdist1<-gdistsamp(~1,~1,~1,gfgrsp,keyfun = "hazard",output = "density")
grspgfsoiltype<-gdistsamp(~SiteName,~1,~1,gfgrsp,keyfun = "hazard",output = "density")
grspgftopovar<-gdistsamp(~toposcale,~1,~1,gfgrsp,keyfun = "hazard",output = "density")
grspgfsoiltypetopo<-gdistsamp(~wetscale,~1,~1,gfgrsp,keyfun = "hazard",output = "density")
grspgfTRT<-gdistsamp(~TRT,~1,~1,gfgrsp,keyfun = "hazard",output = "density")

grspfit<-fitList(grspgfdist1, grspgfsoiltype, grspgftopovar, grspgfsoiltypetopo, grspgfTRT)
modSel(grspfit) #Review model list to determine which key function fits data best...THIS MAY BE DIFFERENT FOR EACH SPECIES

grsppre<-predict(grspgfdist1, type="lambda", appendData=T)
View(grsppre)

savssoil<- data.frame(SiteName=factor(c("Loamy", "Sandy", "Shallow Gravel", "Thin Loamy", "Very Shallow")))

grsptrt<- data.frame(TRT=factor(c("PBG20", "PBG40", "SLG", "WAG")))

(grspltrt <- predict(grspgfTRT, type="lambda", newdata=grsptrt,
                     appendData=TRUE))

(savssoil <- predict(savsgfsoiltype, type="lambda", newdata=savssoil,
                     appendData=TRUE))


