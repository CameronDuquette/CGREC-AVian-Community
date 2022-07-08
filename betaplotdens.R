commbeta=read.csv("distbetas.csv")

plotcom=ggplot(commbeta, aes(x=Var, y=Density)) + geom_point(aes(size=0.8))+
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1) +
  geom_line()+theme_classic()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Treatment")+ theme(legend.position="none")+theme(text = element_text(size=8, face = "bold"))+ facet_wrap( ~ Species, ncol=3,scales="free", drop = TRUE)+theme(strip.text.x = element_text(angle = 0, size = 9, face = "bold"),strip.background = element_blank())+theme(plot.margin=unit(c(.5,1,.5,.5),"cm"))

soilbeta=read.csv("soilbetas.csv")
plotsoil=ggplot(soilbeta, aes(x=Var, y=Density)) + geom_point(aes(size=0.8))+
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1) +
  geom_line()+theme_classic()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("Soil Type")+theme(legend.position="none")+theme(text = element_text(size=9, face = "bold"))+ facet_wrap( ~ Species, ncol=2,scales="free", drop = TRUE)+theme(strip.text.x = element_text(angle = 0, size = 13, face = "bold"),strip.background = element_blank())+theme(plot.margin=unit(c(.5,1,.5,.5),"cm"))
