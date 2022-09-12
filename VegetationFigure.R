remotes::install_github("jfq3/ggordiplots" )
pacman::p_load(ggordiplots)


dissertveg=read.csv("Final2020vegagg1.csv")#this is new file with wag changed to mtrg

dissertveg$TRYR=paste(dissertveg$Transect, dissertveg$YEAR, sep = "")

colnames(dissertveg)

###2017
library(vegan)
library(ggplot2)
dissertveg17=subset(dissertveg,YEAR=="2017")
disvegstr17=dissertveg17[,c(15,16,17,18)]##stopped here

dis17 <- vegdist(disvegstr17,method="bray")
mod17 <- betadisper(dis17, dissertveg17$pasyr, bias.adjust = TRUE)
disavg17=as.data.frame(mod17$group)
disavg17$dist=mod17$distances

dis17avg=read.csv("vegdst17.csv")
mod17=aov(dist~TRT, data=dis17avg)
(mod17.HSD <- TukeyHSD(mod17))

vegstreet17=metaMDS(dissertveg17[,c(15,16,17,18)], distance="bray", k=2, trymax=5000)
plot(vegstreet17)

strucdist17 <- vegdist(dissertveg17[,c(15,16,17,18)],method="bray")
modstruc17 <- betadisper(strucdist17, dissertveg17$TRT)

scrs17<-scores(vegstreet17, display='sites')
scrs17<-cbind(as.data.frame(scrs17), Treatment=dissertveg17$TRT)

species.scores17 <- as.data.frame(scores(vegstreet17, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores17$species <- rownames(species.scores17)  # create a column of species, from the rownames of species.scores

cent17=aggregate(cbind(NMDS1, NMDS2)~Treatment, data=scrs17, FUN=mean)

segs17 <- merge(scrs17, setNames(cent17, c('Treatment','oNMDS1','oNMDS2')),
              by = 'Treatment', sort = FALSE)

segs17$lineColor <- ifelse(segs17$Treatment=="PBG20", "red",
                    ifelse(segs17$Treatment=="PBG40", "green",
                    ifelse(segs17$Treatment=="SLG", "blue","purple")))
scrs17$lineColor <- ifelse(scrs17$Treatment=="PBG20", "red",
                           ifelse(scrs17$Treatment=="PBG40", "green",
                                  ifelse(scrs17$Treatment=="SLG", "blue","purple")))

vgg17=ggplot(scrs17, aes(x = NMDS1, y = NMDS2, color=Treatment)) +
  geom_segment(data = segs17,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent17, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.text.x= element_text(colour="black",size=10))+
  theme(axis.text.y= element_text(colour="black",size=10))+
  
  geom_text_repel(data=species.scores17,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15)+scale_color_manual(values = c("#00A5FF", "#E7B800", "#FC4E07"))
  
####2018
library(vegan)
library(ggplot2)
dissertveg18=subset(dissertveg,YEAR=="2018")
disvegstr18=dissertveg18[,c(15,16,17,18)]##stopped here


vegstreet18=metaMDS(dissertveg18[,c(15,16,17,18)], distance="bray", k=2, trymax=5000)
plot(vegstreet18)

dis18 <- vegdist(disvegstr18,method="bray")
mod18 <- betadisper(dis18, dissertveg18$pasyr, bias.adjust = TRUE)
disavg18=as.data.frame(mod18$group)
disavg18$dist=mod18$distances

dis18avg=read.csv("C:/Users/cduquett/Dropbox/aviancommunity/Revisions/vegdst18.csv")
mod18=aov(dist~TRT, data=dis18avg)
(mod18.HSD <- TukeyHSD(mod18))

strucdist18 <- vegdist(dissertveg18[,c(15,16,17,18)],method="bray")
modstruc18 <- betadisper(strucdist18, dissertveg18$TRT)

scrs18<-scores(vegstreet18, display='sites')
scrs18<-cbind(as.data.frame(scrs18), Treatment=dissertveg18$TRT)

species.scores18 <- as.data.frame(scores(vegstreet18, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores18$species <- rownames(species.scores18)  # create a column of species, from the rownames of species.scores

cent18=aggregate(cbind(NMDS1, NMDS2)~Treatment, data=scrs18, FUN=mean)

segs18 <- merge(scrs18, setNames(cent18, c('Treatment','oNMDS1','oNMDS2')),
                by = 'Treatment', sort = FALSE)

segs18$lineColor <- ifelse(segs18$Treatment=="PBG20", "red",
                           ifelse(segs18$Treatment=="PBG40", "green",
                                  ifelse(segs18$Treatment=="SLG", "blue","purple")))
scrs18$lineColor <- ifelse(scrs18$Treatment=="PBG20", "red",
                           ifelse(scrs18$Treatment=="PBG40", "green",
                                  ifelse(scrs18$Treatment=="SLG", "blue","purple")))
vgg18=ggplot(scrs18, aes(x = NMDS1, y = NMDS2, color=Treatment)) +
  geom_segment(data = segs18,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent18, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.text.x= element_text(colour="black",size=10))+
  theme(axis.text.y= element_text(colour="black",size=10))+
  geom_text_repel(data=species.scores18,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15)+scale_color_manual(values = c("#00BC59","#00A5FF", "#E7B800", "#FC4E07"))


###2019
library(vegan)
library(ggplot2)
dissertveg19=subset(dissertveg,YEAR=="2019")

disvegstr19=dissertveg19[,c(15,16,17,18)]##stopped here

dis19 <- vegdist(disvegstr19,method="bray")
mod19 <- betadisper(dis19, dissertveg19$pasyr, bias.adjust = TRUE)
disavg19=as.data.frame(mod19$group)
disavg19$dist=mod19$distances

dis19avg=read.csv("vegdst19.csv")
mod19=aov(dist~TRT, data=dis19avg)
(mod19.HSD <- TukeyHSD(mod19))

vegstreet19=metaMDS(dissertveg19[,c(15,16,17,18)], distance="bray", k=2, trymax=5000)
plot(vegstreet19)


strucdist19 <- vegdist(dissertveg19[,c(15,16,17,18)],method="bray")
modstruc19 <- betadisper(strucdist19, dissertveg19$TRT)

scrs19<-scores(vegstreet19, display='sites')
scrs19<-cbind(as.data.frame(scrs19), Treatment=dissertveg19$TRT)

species.scores19 <- as.data.frame(scores(vegstreet19, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores19$species <- rownames(species.scores19)  # create a column of species, from the rownames of species.scores

cent19=aggregate(cbind(NMDS1, NMDS2)~Treatment, data=scrs19, FUN=mean)

segs19 <- merge(scrs19, setNames(cent19, c('Treatment','oNMDS1','oNMDS2')),
                by = 'Treatment', sort = FALSE)

segs19$lineColor <- ifelse(segs19$Treatment=="PBG20", "red",
                           ifelse(segs19$Treatment=="PBG40", "green",
                                  ifelse(segs19$Treatment=="SLG", "blue","purple")))
scrs19$lineColor <- ifelse(scrs19$Treatment=="PBG20", "red",
                           ifelse(scrs19$Treatment=="PBG40", "green",
                                  ifelse(scrs19$Treatment=="SLG", "blue","purple")))
vgg19=ggplot(scrs19, aes(x = NMDS1, y = NMDS2, color=Treatment)) +
  geom_segment(data = segs19,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent19, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.text.x= element_text(colour="black",size=10))+
  theme(axis.text.y= element_text(colour="black",size=10))+
  geom_text_repel(data=species.scores19,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15)+scale_color_manual(values = c("#00BC59","#00A5FF", "#E7B800", "#FC4E07"))
###2020

library(vegan)
library(ggplot2)
dissertveg20=subset(dissertveg,YEAR=="2020")
disvegstr20=dissertveg20[,c(15,16,17,18)]##stopped here

dis20 <- vegdist(disvegstr20,method="bray")
mod20 <- betadisper(dis20, dissertveg20$pasyr, bias.adjust = TRUE)
disavg20=as.data.frame(mod20$group)
disavg20$dist=mod20$distances
dis20avg=read.csv("vegdst20.csv")
mod20=aov(dist~TRT, data=dis20avg)
(mod20.HSD <- TukeyHSD(mod20))

vegstreet20=metaMDS(dissertveg20[,c(15,16,17,18)], distance="bray", k=2, trymax=5000)
plot(vegstreet20)

strucdist20 <- vegdist(dissertveg20[,c(15,16,17,18)],method="bray")
modstruc20 <- betadisper(strucdist20, dissertveg20$TRT)

scrs20<-scores(vegstreet20, display='sites')
scrs20<-cbind(as.data.frame(scrs20), Treatment=dissertveg20$TRT)

species.scores20 <- as.data.frame(scores(vegstreet20, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores20$species <- rownames(species.scores20)  # create a column of species, from the rownames of species.scores

cent20=aggregate(cbind(NMDS1, NMDS2)~Treatment, data=scrs20, FUN=mean)

segs20 <- merge(scrs20, setNames(cent20, c('Treatment','oNMDS1','oNMDS2')),
                by = 'Treatment', sort = FALSE)

segs20$lineColor <- ifelse(segs20$Treatment=="PBG20", "red",
                           ifelse(segs20$Treatment=="PBG40", "green",
                                  ifelse(segs20$Treatment=="SLG", "blue","purple")))
scrs20$lineColor <- ifelse(scrs20$Treatment=="PBG20", "red",
                           ifelse(scrs20$Treatment=="PBG40", "green",
                                  ifelse(scrs20$Treatment=="SLG", "blue","purple")))

vgg20=ggplot(scrs20, aes(x = NMDS1, y = NMDS2, color=Treatment)) +
  geom_segment(data = segs20,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent20, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+theme_classic()+
  theme(axis.text.x= element_text(colour="black",size=10))+
  theme(axis.text.y= element_text(colour="black",size=10))+
  geom_text_repel(data=species.scores20,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black", max.overlaps = 15)+scale_color_manual(values = c("#00BC59","#00A5FF", "#E7B800", "#FC4E07"))

veggrid=grid.arrange(vgg17,vgg18,vgg19,vgg20, nrow = 2)

