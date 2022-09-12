library(vegan)
maxcomb2=read.csv("rawcom.csv")

####hullyuh includes all species. use for diversity stuff
ordbirds2=maxcomb2[2:17]
birdmax2=metaMDS(ordbirds2, distance="bray", k=3, trymax=1000)

birdenv=maxcomb2[,c(18:ncol(maxcomb2))]
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
    ROUGH = toposd15_var_2.x,
    WET = Shape.Area
    
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
library(ggplot2)
nmds.plot.dune <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(TRT), size = 2))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "TRT")#+ # add legend labels for Management and Landuse
  #theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
  
  nmds.plot.dune + labs(title = "Basic ordination plot") #displays plot
##########
  library(ggrepel)
nmds.plot.dune+
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="Ordination with environmental vectors")

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
  geom_path(data = df_ell.dune.management, aes(x = NMDS1, y = NMDS2, group = TRT))+#this is the ellipse, seperate ones by Site. 
  geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species, fontface='bold'),color="black",alpha=0.5, max.overlaps = 15)


#betadisperson test
dis <- vegdist(ordbirds2,method="bray")
mod <- betadisper(dis, envtest2$TRT)
anova(mod)
(mod.HSD <- TukeyHSD(mod))

adonis(ordbirds2 ~ SiteName.x+TRT.x+toposd15_var_2.x+Shape.Area+YEAR.x+R.x+ld.x+meanNWood.x+meanNForb.x+meanDead.x+meanpopr.x+meanbrin.x, data = birdenv, by = NULL)
pairwise.adonis(ordbirds2, factors=birdenv$SiteName.x)



