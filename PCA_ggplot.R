library(ggplot2)
setwd("/Users/robbi/Box/Robbi_PhD/02_Inka_Royal_Ancestry/IncaModern/EAG/") #set working directory

fn = "/Users/robbi/Box/Robbi_PhD/02_Inka_Royal_Ancestry/IncaModern/EAG/pca_1.pca.evec.txt"
evecDat1 = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "Pop"))

#make colour palette
colours <- c("blue", "#E2AC36", "#B3493A", "#47271C", "#BCAB7E","#BCAB7E")
colours <- colours[as.factor(evecDat1$Pop)]  

#define symbols
shapes <- c(0,0,3,0,0)
shapes <- shapes[as.factor(evecDat1$Pop)]  

#base plot 1v2
plot_1v2 <- (
  plot(evecDat1$PC1, evecDat1$PC2, col = colours, xlab="PC1", ylab="PC2", pch = shapes, cex = 2, lwd=1.5,
       cex.axis=1.3, cex.lab=1.2, las=1))
legend("topleft", legend = levels(as.factor(evecDat1$Pop)), col = c("blue", "#E2AC36", "#B3493A", "#47271C", "#BCAB7E","#BCAB7E"),
       pch = c(0,0,3,0,0), cex = 1.2, bty="o", pt.cex=2, pt.lwd=1.5, scale_fill_discrete(labels = labels))
#ggplot 1v2
PC1_2 <- ggplot(evecDat1, #dataset to plot
              aes(x = PC1, #x-axis is petal width
                  y = PC2, #y-axis is petal length
                  color = Pop)) + #each species is represented by a different shape
  geom_point() + #default scatter plot
  theme_light() + #light theme with no grey background and with grey lines and axes
  scale_x_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on top 
  scale_y_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on right
  theme(panel.grid.major = element_blank(), #remove major gridline
        panel.grid.minor = element_blank(), #remove minor gridline
#        legend.justification = c(0, 0), #justification is centered by default, c(1, 0) means bottom right
#        legend.position = c(0.97, 0.01), #position relative to justification
#        legend.background = element_rect(color = "grey"), #legend box with grey lines
#        legend.text = element_text(face = "italic"), #since the legend is species names, display in italics
        axis.title.x.top = element_blank(), #remove top x-axis title
        axis.text.x.top = element_blank(), #remove labels on top x-axis
        axis.title.y.right = element_blank(), #remove right y-axis title
        axis.text.y.right = element_blank()) + #remove labels on right y-axis
  scale_shape(labels = c("Mbuti", "Karitiana", "Han", "Basque", "IncaDescendant"), #edit legend labels
              solid = FALSE) + #force hollow points to increase clarity in case of overlaps
  labs(x = "PC1 (7.633 %)", #edit x-axis label
       y = "PC2 (3.975 %)", #edit y-axis label
       shape = "Population") #edit legend title
PC1_2 #show plot

#ggplot 1v2
PC1_3 <- ggplot(evecDat1, #dataset to plot
                aes(x = PC1, #x-axis is petal width
                    y = PC3, #y-axis is petal length
                    color = Pop)) + #each species is represented by a different shape
  geom_point() + #default scatter plot
  theme_light() + #light theme with no grey background and with grey lines and axes
  scale_x_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on top 
  scale_y_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on right
  theme(panel.grid.major = element_blank(), #remove major gridline
        panel.grid.minor = element_blank(), #remove minor gridline
        legend.justification = c(0, 0), #justification is centered by default, c(1, 0) means bottom right
        #        legend.position = c(0.97, 0.01), #position relative to justification
        legend.background = element_rect(color = "grey"), #legend box with grey lines
        #        legend.text = element_text(face = "italic"), #since the legend is species names, display in italics
        axis.title.x.top = element_blank(), #remove top x-axis title
        axis.text.x.top = element_blank(), #remove labels on top x-axis
        axis.title.y.right = element_blank(), #remove right y-axis title
        axis.text.y.right = element_blank()) + #remove labels on right y-axis
  scale_shape(labels = c("Mbuti", "Karitiana", "Han", "Basque", "IncaDescendant"), #edit legend labels
              solid = FALSE) + #force hollow points to increase clarity in case of overlaps
  labs(x = "PC1 (7.633 %)", #edit x-axis label
       y = "PC3 (2.650 %)", #edit y-axis label
       shape = "Population") #edit legend title
PC1_3 #show plot

#composite 
library(patchwork)
PC1_2_patch <- PC1_2 + #prepare a version of gg3 for patchwork design
      theme(legend.position="none") #remove legend
PC1_2_patch #show plot

patch <- PC1_2_patch + PC1_3 + #assemble patchwork
  plot_layout(widths = c(3, 3)) #set width ratio of the 3 panels
#  plot_annotation(tag_levels = "A") #add plot labels (uppercase Latin letters)
patch #show plot

ggsave("patch.pdf", width = 12, height = 6) #save in pdf format with size 12 x 6 in
