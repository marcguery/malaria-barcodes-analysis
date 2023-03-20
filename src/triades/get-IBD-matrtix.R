triades <- read.table("out/triades.tsv", header = T)[,c(1:3)]

lookuptriades <- expand.grid(c(triades[,1], triades[,2], triades[,3]),
                                   c(triades[,1], triades[,2], triades[,3]))
lookuptriades[,1] <- as.character(lookuptriades[,1])
lookuptriades[,2] <- as.character(lookuptriades[,2])
lookuptriades <- data.frame(pmin(lookuptriades[,1], lookuptriades[,2]),
                   pmax(lookuptriades[,1], lookuptriades[,2]))

lookuptriades <- lookuptriades[!duplicated(lookuptriades),]
colnames(lookuptriades) <- c("ID1", "ID2")

lookuptriades <- merge(edges, lookuptriades, all.y = TRUE)
lookuptriades.na <- lookuptriades[is.na(lookuptriades$seq_shared_best_traj),c(1,2)]
lookuptriades.na <- merge(edges, lookuptriades.na, 
                       by.x = c("ID2", "ID1"), by.y = c("ID1", "ID2")
                       )

lookuptriades <- rbind(lookuptriades[!is.na(lookuptriades$seq_shared_best_traj),], lookuptriades.na)

##ORDER according to parents
childrenmanualorder <- c("J0060701_1511", "J0262731_1412", "K0030301_1504", "K0030302_1512",
                       "K0161621_1506", "K0374408_1412","K0374413_1611", "K0060617_1512", 
                        "P0020206_1702", "P0020206_1703", "K0171738_1412")
triades[triades$ID=="K0374413_1611",c(2,3)] <- triades[triades$ID=="K0374413_1611",c(3,2)]
#triades <- triades[-which(triades[,1]%in%c("K0060617_1512", "P0020206_1703")),]
parents <- triades[order(factor(triades[,1], 
                                levels = childrenmanualorder)),c(2,3)]
parents <- data.frame(gdata::interleave(matrix(parents[,1]), matrix(parents[,2])))
colnames(parents) <- "ID"
parents$number <- c(1:nrow(parents))

children <- data.frame(triades[order(factor(triades[,1], 
                                            levels = childrenmanualorder)),1])
colnames(children) <- "ID"
children$number <- c((nrow(parents)+1):(nrow(parents)+nrow(children)))



samples <-  gdata::interleave(matrix(parents$number), 
                              matrix(gdata::interleave(matrix(children$number), matrix(children$number), drop = TRUE)), 
                              drop = TRUE)
samples <- samples[!duplicated(samples)]
samples.df <- data.frame(samples)
colnames(samples.df) <- "number"
samples.df$ID <- NA
samples.df$ID[samples.df$number%in%parents$number] <- parents$ID
samples.df$ID[samples.df$number%in%children$number] <- children$ID

lookuptriades <- merge(lookuptriades, samples.df, by.x = "ID1", by.y = "ID")
lookuptriades <- merge(lookuptriades, samples.df, by.x = "ID2", by.y = "ID", 
                       suffixes = c("1", "2"))

lookuptriades$number1 <- factor(lookuptriades$number1,
                                levels = samples)
lookuptriades$number2 <- factor(lookuptriades$number2,
                                levels = samples)

lookuptriades.short <- lookuptriades[,c(1,2,4,27,ncol(lookuptriades)-1,ncol(lookuptriades))]
colnames(lookuptriades.short) <- c("ID2", "ID1", "seq_shared_best_traj", "family", "number1", "number2")
lookuptriades.short.rev <- lookuptriades[,c(2,1,4,27,ncol(lookuptriades),ncol(lookuptriades)-1)]
colnames(lookuptriades.short.rev) <- c("ID2", "ID1", "seq_shared_best_traj", "family", "number1", "number2")
lookuptriades.short <- rbind(lookuptriades.short, lookuptriades.short.rev)

triades.lines <-  data.frame(seq(1,nrow(triades)*3,3)-0.5)
colnames(triades.lines) <- "rowNumber"
triades.lines$index <- c(1:nrow(triades.lines))


gg <- ggplot(lookuptriades.short)+
  geom_tile(aes(x = number1, 
                y = number2,
                fill = factor(floor((seq_shared_best_traj*0.999999)*5)/5)),
            show.legend = FALSE)+
  geom_text(data = lookuptriades.short[round(lookuptriades.short$seq_shared_best_traj, 2)>0,],
            aes(x = number1, 
                y = number2,
                label = as.character(sprintf("%03s",round(seq_shared_best_traj, 2)*100))),
            size = 4)+
  geom_tile(data = triades.lines,
            aes(x = rowNumber+1.5, y = rowNumber+1.5), fill = "black", width = 1, height = 1)+
  scale_fill_brewer("IBD", type = "seq", palette = 7, labels = c("0-0.2","0.2-0.4","0.4-0.6",
                                                                 "0.6-0.8","0.8-1"))+
  scale_x_discrete(labels = samples.df$ID)+
  scale_y_discrete(labels = samples.df$ID)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90))
gg
ggsave("out/IBD-triades.png", width = 13, height = 13)

gg <- ggplot(lookuptriades.short[lookuptriades.short$ID1%in%c(childrenmanualorder) & lookuptriades.short$ID2%in%c(childrenmanualorder),])+
  geom_tile(aes(x = factor(ID1, levels = childrenmanualorder), 
                y = factor(ID1, levels = childrenmanualorder)),
            fill = "grey80", show.legend = FALSE)+
  geom_tile(aes(x = factor(ID1, levels = childrenmanualorder), 
                y = factor(ID2, levels = childrenmanualorder),
                fill = factor(floor((seq_shared_best_traj*0.999999)*5)/5)),
            show.legend = FALSE)+
  geom_text(data = lookuptriades.short[lookuptriades.short$ID1%in%c(childrenmanualorder) & lookuptriades.short$ID2%in%c(childrenmanualorder) & round(lookuptriades.short$seq_shared_best_traj, 2)>0,],
            aes(x = ID1, 
                y = ID2,
                label = as.character(sprintf("%03s",round(seq_shared_best_traj, 2)*100))),
            size = 8)+
  scale_fill_brewer("IBD", type = "seq", palette = "Reds", labels = c("0-0.2","0.2-0.4","0.4-0.6",
                                                                 "0.6-0.8","0.8-1"))+
  xlab("")+
  ylab("")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        text = element_text(size = 20))
gg
ggsave("out/offsprings/IBD-children.png", width = 13, height = 13)

allchildren$unknownbases <- sum(chromlength$length)-apply(allchildren[,c("misbases", "parent1bases", "parent2bases")], 1, sum)+allchildren[,c("samebases")]
allchildren$parent1onlybases <- allchildren$parent1bases-allchildren$samebases
allchildren$parent2onlybases <- allchildren$parent2bases-allchildren$samebases
allchildren.melt <- melt(allchildren, id.vars = c("ID", "parent1", "parent2", "parent1bases", "parent2bases",
                                                  "parentsIBD", "parent1IBD", "parent2IBD",
                                                  "family"))
allchildren.melt$variable <- factor(allchildren.melt$variable,
                                    levels = rev(c("parent1onlybases",
                                               "samebases",
                                               "parent2onlybases",
                                               "misbases",
                                               "unknownbases")))
gg <- ggplot(allchildren.melt)+
  geom_col_pattern(aes(y = factor(ID, levels = childrenmanualorder), x = value, 
               fill = factor(variable),
               pattern = factor(variable),
               pattern_fill = factor(variable),
               pattern_density = factor(variable)),
           position = "stack", pattern_colour = NA, pattern_spacing = 0.025)+
  geom_text(data = allchildren,
            aes(x = parent1onlybases/2, y = ID, label = parent1),
            color = "white", fontface = "bold", size = 6)+
  geom_text(data = allchildren,
            aes(x = parent1onlybases+samebases+parent2onlybases/2, y = ID, label = parent2),
            color = "black", fontface = "bold", size = 6)+
  scale_pattern_manual("Offspring genome origin",
                       breaks = c("parent1onlybases",
                                  "samebases",
                                  "parent2onlybases",
                                  "misbases",
                                  "unknownbases"),
                       labels = c("Parent 1", "Parents 1 & 2",
                                  "Parent 2",
                                  "Different origin",
                                  "Unknown origin"),
                       values = c("none", "stripe", "none", "none", "none"))+
  scale_pattern_density_manual("Offspring genome origin",
                               breaks = c("parent1onlybases",
                                  "samebases",
                                  "parent2onlybases",
                                  "misbases",
                                  "unknownbases"),
                               labels = c("Parent 1", "Parents 1 & 2",
                                          "Parent 2",
                                          "Different origin",
                                          "Unknown origin"),
                       values = c(0,0.5, 0, 0, 0))+
  scale_pattern_fill_manual("Offspring genome origin",
                            breaks = c("parent1onlybases",
                                  "samebases",
                                  "parent2onlybases",
                                  "misbases",
                                  "unknownbases"),
                            labels = c("Parent 1", "Parents 1 & 2",
                                       "Parent 2",
                                       "Different origin",
                                       "Unknown origin"),
                       values = c(NA, "gold1", "gold1", "pink", "grey80"))+
  scale_fill_manual("Offspring genome origin",
                    breaks = c("parent1onlybases",
                                "samebases",
                                "parent2onlybases",
                                "misbases",
                                "unknownbases"),
                    labels = c("Parent 1", "Parents 1 & 2",
                               "Parent 2",
                               "Different origin",
                               "Unknown origin"),
                    values = c("blue", "blue", "gold1", "pink", "grey80"))+
  scale_x_continuous(n.breaks = 20,
                     labels = function(x){x/1000},
                     expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  xlab("Genomic size (kb)")+
  ylab("Offspring genome")+
  theme(legend.position = "top",
        panel.grid = element_blank(),
        text = element_text(size = 18),
        legend.key.size = unit(1.25, 'cm'))
gg
png("out/offsprings/genomic-share.png", 
    width = 20, height = 8,units = "in", res = 150)
plot(gg)
dev.off()


###########STATS#############
1-round((allchildren$unknownbases[allchildren$parent1=="K0141403_1504" & 
                                 allchildren$parent2=="K0374502_1412"]+
         allchildren$misbases[allchildren$parent1=="K0141403_1504" & 
                                    allchildren$parent2=="K0374502_1412"])/sum(chromlength$length),3)
1-round(allchildren$samebases[allchildren$parent1=="K0141403_1504" & 
                                    allchildren$parent2=="K0374502_1412"]/(sum(chromlength$length)-(allchildren$unknownbases[allchildren$parent1=="K0141403_1504" & 
                                                                                                      allchildren$parent2=="K0374502_1412"]+
                                                                             allchildren$misbases[allchildren$parent1=="K0141403_1504" & 
                                                                                                    allchildren$parent2=="K0374502_1412"])),3)

1-round((allchildren$unknownbases[!(allchildren$parent1=="K0141403_1504" & 
                                    allchildren$parent2=="K0374502_1412")]+
           allchildren$misbases[!(allchildren$parent1=="K0141403_1504" & 
                                  allchildren$parent2=="K0374502_1412")])/sum(chromlength$length),3)
1-round(allchildren$samebases[!(allchildren$parent1=="K0141403_1504" & 
                                allchildren$parent2=="K0374502_1412")]/(sum(chromlength$length)-(allchildren$unknownbases[!(allchildren$parent1=="K0141403_1504" & 
                                                                                                                           allchildren$parent2=="K0374502_1412")]+
                                                                                                  allchildren$misbases[!(allchildren$parent1=="K0141403_1504" & 
                                                                                                                         allchildren$parent2=="K0374502_1412")])),3)
############################