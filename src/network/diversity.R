diversity <- edges[edges$date1 < "2017-01-01" & edges$date2 < "2017-01-01" & edges$fract_sites_IBD > -1,]
diversity$kept <- diversity$fract_sites_IBD > 0.5
diversity$group <- sapply(diversity$fract_sites_IBD, FUN = function(x){
  if(x < 0.5){
    return(-1)
  }
  if(x >= 0.5 & x < 0.9){
    return(0)
  }
  if(x >= 0.9){
    return(1)
  }
})

gg <- ggplot(diversity)+
  geom_bar(aes(x = floor(fract_sites_IBD*20)/20,
               fill = factor(group)),
           show.legend = T, color = "black", size = 0.25, width = 1/20)+
  geom_hline(yintercept = 500, linetype = 2, color = "red", size = 0.25)+
  scale_x_continuous(expand = c(0,0), n.breaks = 10)+
  scale_y_continuous(expand = c(0,0), n.breaks = 10)+
  xlab("IBD")+
  ylab("Number of pairs")+
  scale_fill_manual(name = "",
                    values = c("white", "grey60", "black"),
                    labels = c("Unrelated", "Related", "Highly related"))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 16),
        legend.position = "top")
gg

ggsave("out/diversity-IBD.png", width = 7, height = 3.5)

gg <- ggplot(diversity)+
  geom_bar(aes(x = floor(fract_sites_IBD*20)/20,
               y = after_stat(ifelse(count>=500, 501, count)),
               fill = factor(group)),
           show.legend = T, color = "black", size = 0.25, width = 1/20)+
  geom_hline(yintercept = 500, linetype = 2, color = "red", size = 0.25)+
  scale_x_continuous(expand = c(0,0), n.breaks = 10)+
  scale_y_continuous(expand = c(0,0), n.breaks = 5)+
  xlab("IBD")+
  ylab("Number of pairs")+
  scale_fill_manual(name = "",
                    values = c("white", "grey60", "black"),
                    labels = c("Unrelated", "Related", "Highly related"))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 16),
        legend.position = "top")
gg
ggsave("out/diversity-IBD-zoom.png", width = 7, height = 3.5)

nrow(diversity[diversity$fract_sites_IBD > 0.5,])/
  nrow(diversity)
