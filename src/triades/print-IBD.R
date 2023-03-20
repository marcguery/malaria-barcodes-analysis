
##################FUNCTIONS###############
sets<-function(start, end, group, overlap=length(unique(group))) {
  dd<-rbind(data.frame(pos=start, event=1), data.frame(pos=end, event=-1))
  dd<-aggregate(event~pos, dd, sum)
  dd<-dd[order(dd$pos),]
  dd$open <- cumsum(dd$event)
  r<-rle(dd$open>=overlap)
  ex<-cumsum(r$lengths-1 + rep(1, length(r$lengths))) 
  sx<-ex-r$lengths+1
  data.frame(start = dd$pos[sx[r$values]],
             end = dd$pos[ex[r$values]+1])
} 

parentsdiff <- function(ibddata, child, parent1, parent2){
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.1$id <- 1
  subdata.1 <- subdata.1[subdata.1$different==1,]
  subdata.2 <- subduo(ibddata, child, parent2)
  subdata.2$id <- 2
  subdata.2 <- subdata.2[subdata.2$different==1,]
  subdata <- rbind(subdata.1, subdata.2)
  if(nrow(subdata)==0){
    return(subdata[,c("chr", "start", "end")])
  }
  subdata%>%
    group_by(chr)%>%
    do(sets(.$start, .$end, .$id, 2))
}
parentssame <- function(ibddata, child, parent1, parent2){
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.1$id <- 1
  subdata.1 <- subdata.1[subdata.1$different==0,]
  subdata.2 <- subduo(ibddata, child, parent2)
  subdata.2$id <- 2
  subdata.2 <- subdata.2[subdata.2$different==0,]
  subdata <- rbind(subdata.1, subdata.2)
  if(nrow(subdata)==0){
    return(subdata[,c("chr", "start", "end")])
  }
  subdata%>%
    group_by(chr)%>%
    do(sets(.$start, .$end, .$id, 2))
}

ibdminus <- function(fragments1, fragments2, chrs){
  newfragments <- data.frame(matrix(ncol = 3))
  colnames(newfragments) <- c("chr", "start", "end")
  newfragments = newfragments[FALSE,]
  
  for (chr in chrs){
    Ints1 <- apply(fragments1[fragments1$chr==chr,], MARGIN = 1,
                  FUN = function(x){sets::interval(l=as.numeric(x[4]), r=as.numeric(x[5]))})
    Ints2 <- apply(fragments2[fragments2$chr==chr,], MARGIN = 1,
                   FUN = function(x){sets::interval(l=as.numeric(x[4]), r=as.numeric(x[5]))})
    diffset <- sets::interval_intersection(Ints1, sets::interval_symdiff(Ints2, Ints1))
    if (length(diffset) == 0){
      next
    }
    
    
    newfragments.current <- data.frame(t(as.matrix(sapply(diffset, range))))
    colnames(newfragments.current) <- c("start", "end")
    newfragments.current$chr <- chr
    newfragments <- rbind(newfragments, newfragments.current)
  }
  return(newfragments)
}
ibdcommon<- function(fragments1, fragments2, chrs){
  newfragments <- data.frame(matrix(ncol = 3))
  colnames(newfragments) <- c("chr", "start", "end")
  newfragments = newfragments[FALSE,]
  
  for (chr in chrs){
    Ints1 <- apply(fragments1[fragments1$chr==chr,], MARGIN = 1,
                   FUN = function(x){sets::interval(l=as.numeric(x[4]), r=as.numeric(x[5]))})
    Ints2 <- apply(fragments2[fragments2$chr==chr,], MARGIN = 1,
                   FUN = function(x){sets::interval(l=as.numeric(x[4]), r=as.numeric(x[5]))})
    sameset <- sets::interval_intersection(Ints2, Ints1)
    if (length(sameset) == 0){
      next
    }
    
    
    newfragments.current <- data.frame(t(as.matrix(sapply(sameset, range))))
    colnames(newfragments.current) <- c("start", "end")
    newfragments.current$chr <- chr
    newfragments <- rbind(newfragments, newfragments.current)
  }
  return(newfragments)
}

subduo <- function(ibddata, above, below){
  subdata <- ibddata[ibddata$ID1==above,]
  subdata <- subdata[subdata$ID2==below,]
  if (dim(subdata)[1]==0){
    subdata <- ibddata[ibddata$ID2==above,]
    subdata <- subdata[subdata$ID1==below,]
  }
  if(dim(subdata)[1]==0){
    stop("There is no correlation between ", above, " and ", below)
  }
  return(subdata)
  
}

parentsmap <- function(ibddata, child, parent1, parent2, chlen){
  print(paste(child, parent1, parent2))
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.2 <- subduo(ibddata, child, parent2)
  
  chs <- chlen
  
  gg <- ggplot()+
    geom_rect(data=chs,
              aes(xmin = 0, xmax=length, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
              color="grey60", fill="grey60")+
    geom_text(data=chs,
              aes(x = 0-200000, y=as.numeric(chr)*2.75, 
                  label=chr), 
              color="grey50", size=3)+
    geom_segment(data=subdata.1,
                 aes(x=start, xend=end, 
                     y=as.numeric(chr)*2.75+0.75, yend=as.numeric(chr)*2.75+0.75, 
                     color=factor(different)), size=1,alpha=0.5)+
    geom_segment(data=subdata.2,
                 aes(x=start, xend=end, 
                     y=as.numeric(chr)*2.75-0.75, yend=as.numeric(chr)*2.75-0.75, 
                     color=factor(different)), size=1,alpha=0.5)+
    scale_color_manual(values=c("0"="red", "1"="blue"))+
    geom_rect(data=subdata.1[subdata.1$different==1,],
              aes(xmin=start, xmax=end, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5, 
                  fill=factor(different)), size=0.1,alpha=0.5)+
    geom_rect(data=subdata.2[subdata.2$different==1,],
              aes(xmin=start, xmax=end, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5, 
                  fill=factor(different)), size=0.1,alpha=0.5)+
    geom_rect(data=subdata.1[subdata.1$different!=1,],
              aes(xmin=start, xmax=end, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5, 
                  fill=factor(different)), size=0.1,alpha=0.5)+
    geom_rect(data=subdata.2[subdata.2$different!=1,],
              aes(xmin=start, xmax=end, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5, 
                  fill=factor(different)), size=0.1,alpha=0.5)+
    scale_fill_manual(values=c("0"="yellow", "1"="cyan"))+
    labs(title = paste(child, "from", parent1, "(top) and", parent2, "(bottom)"))+
    theme(plot.background = element_blank(),panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          legend.position = "none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  return(gg)
}

parentsmapv2 <- function(ibddata, child, parent1, parent2, chlen, mode=0){
  print(paste(child, parent1, parent2))
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.2 <- subduo(ibddata, child, parent2)
  bothdiff <- parentsdiff(ibddata, child, parent1, parent2)
  bothsame <- parentssame(ibddata, child, parent1, parent2)
  
  chs <- chlen
  if(mode==0){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.1[subdata.1$different==1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "gold1", size=0.1)+
      geom_rect(data=subdata.2[subdata.2$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75),
                color = NA, fill = "gold1", size=0.1)+
      geom_rect(data=subdata.2[subdata.2$different==1,],
                aes(xmin=start, xmax=end,
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75),
                color = NA, fill = "blue", size=0.1)+
      geom_rect(data=subdata.1[subdata.1$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75, ymax=as.numeric(chr)*2.75+0.5), 
                fill = "blue", color = NA, size=0.1)+
      geom_rect(data=bothsame,
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                fill = "lightgreen", color = NA, size=0.1)+
      geom_rect(data=bothdiff,
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", size=0.1)+
      labs(title = paste(child, "from", parent1, "(blue) and", parent2, "(yellow)"))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", size = 0.5))
    return(gg)
  }
  
  if(mode==1){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.1[subdata.1$different!=1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "blue", size=0.1)+
      geom_rect(data=subdata.1[subdata.1$different==1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", size=0.1)+
      labs(title = paste(child, "with", parent1))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", size = 0.5))
    return(gg)
    
  }
  
  if(mode==2){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.2[subdata.2$different!=1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "gold1", size=0.1)+
      geom_rect(data=subdata.2[subdata.2$different==1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", size=0.1)+
      labs(title = paste(child, "with", parent2))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", size = 0.5))
    return(gg)
  }
}

parentsmapv3 <- function(ibddata, child, parent1, parent2, chlen, mode=0){
  print(paste(child, parent1, parent2))
  subdata.1 <- subduo(ibddata, child, parent1)
  subdata.2 <- subduo(ibddata, child, parent2)
  bothdiff <- parentsdiff(ibddata, child, parent1, parent2)
  bothsame <- parentssame(ibddata, child, parent1, parent2)
  
  chs <- chlen
  if(mode==0){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.2[subdata.2$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "gold1", size=0.1)+
      geom_rect(data=subdata.1[subdata.1$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                fill = "blue", color = NA, size=0.1)+
      geom_rect_pattern(data=bothsame,
                        aes(xmin=start, xmax=end, 
                            ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                        fill = "blue", pattern_fill ="gold1", pattern_density = 0.5, pattern_spacing = 0.015,
                        color = NA, pattern_color = NA, size=0.1)+
      geom_rect(data=bothdiff,
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", size=0.1)+
      labs(title = paste(child, "from", parent1, "(blue) and", parent2, "(yellow)"))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", size = 0.5))
    return(gg)
  }
  
  if(mode==1){
    gg <- ggplot()+
      geom_rect(data=chs,
                aes(xmin = 0, xmax=length, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
                color=NA, fill="grey80")+
      geom_text(data=chs,
                aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                    label=chr), 
                color="grey50", size=3)+
      geom_rect(data=subdata.1[subdata.1$different!=1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "lightgreen", size=0.1)+
      geom_rect(data=subdata.1[subdata.1$different==1,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = "pink", size=0.1)+
      labs(title = paste(child, "with", parent1))+
      scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                         labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                         minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                         expand = c(0.01,0.01))+
      scale_y_discrete(expand = c(0,0))+
      theme(plot.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
            panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
            panel.background = element_rect(fill = "grey95"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "grey60", size = 0.5))
    return(gg)
    
  }
}

familymap <- function(ibddata, genome, versus, chlen){
  maxother = min(8, length(versus))
  colorlist <- brewer.pal(n = 8, name = "Set1")
  chs <- chlen
  subdata.1 <- subduo(ibddata, genome, versus[1])
  overallibd <- subdata.1
  
  gg <- ggplot()+
    geom_rect(data=chs,
              aes(xmin = 0, xmax=length, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
              color=NA, fill="grey80")+
    geom_text(data=chs,
              aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                  label=chr), 
              color="grey50", size=3)+
    geom_rect(data=subdata.1[subdata.1$different==0,],
              aes(xmin=start, xmax=end, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
              color = NA, fill = colorlist[1], size=0.1)
  
  for (i in c(2:maxother)){
    subdata <- subduo(ibddata, genome, versus[i])
    subdata.same <- subdata[subdata$different == 0,]
    overallibd.same <- overallibd[overallibd$different == 0,]
    newfragments <- ibdminus(subdata.same, overallibd.same, chromlength$chr)
    gg <- gg+
      geom_rect(data=newfragments,
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5),
                color = NA, fill = colorlist[i], size=0.1)
    
    overallibd <- rbind(overallibd, subdata)
  }
  
  gg <- gg+
    scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                       labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                       minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                       expand = c(0.01,0.01))+
    scale_y_discrete(expand = c(0,0))+
    theme(plot.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
          panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
          panel.background = element_rect(fill = "grey95"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(color = "grey60", size = 0.5))
  return(gg)
}

familymapv2 <- function(ibddata, genome, versus, chlen){
  maxother = min(8, length(versus))
  coordsplot <- c(-0.5,c(-0.5)+c(1:maxother)/maxother)
  colorlist <- brewer.pal(n = 8, name = "Set1")
  chs <- chlen
  
  gg <- ggplot()+
    geom_rect(data=chs,
              aes(xmin = 0, xmax=length, 
                  ymin=as.numeric(chr)*2.75-0.5, ymax=as.numeric(chr)*2.75+0.5), 
              color=NA, fill="grey80")+
    geom_text(data=chs,
              aes(x = 0-100000, y=as.numeric(chr)*2.75, 
                  label=chr), 
              color="grey50", size=3)
  
  for (i in c(1:maxother)){
    othergenome <- versus[i]
    if (othergenome == genome){
      selfcomp <- chs
      selfcomp$coord1 <- coordsplot[i]
      selfcomp$coord2 <- coordsplot[i+1]
      gg <- gg+
        geom_rect(data=selfcomp,
                  aes(xmin=0, xmax=length, 
                      ymin=as.numeric(chr)*2.75+coord1, 
                      ymax=as.numeric(chr)*2.75+coord2),
                  color = NA, fill = colorlist[i], size=0.1, 
                  alpha = 0.5)
      next
    }
    
    subdata <- subduo(ibddata, genome, othergenome)
    subdata$coord1 <- coordsplot[i]
    subdata$coord2 <- coordsplot[i+1]
    
    gg <- gg+
      geom_rect(data=subdata[subdata$different==0,],
                aes(xmin=start, xmax=end, 
                    ymin=as.numeric(chr)*2.75+coord1, 
                    ymax=as.numeric(chr)*2.75+coord2),
                color = NA, fill = colorlist[i], size=0.1)
    
  }
  
  gg <- gg+
    scale_x_continuous(breaks = seq(0,3.25*10^6,0.25*10^6), 
                       labels = seq(0,3.25*10^6,0.25*10^6)/10^6,
                       minor_breaks = seq(0,3.25*10^6,0.125*10^6),
                       expand = c(0.01,0.01))+
    scale_y_discrete(expand = c(0,0))+
    theme(plot.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey60", linetype = 2, size = 0.5),
          panel.grid.minor.x = element_line(color = "grey80", linetype = 2, size = 0.5),
          panel.background = element_rect(fill = "grey95"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(color = "grey60", size = 0.5))
  return(gg)
}
#################################

####################READ####################
infogenome <- read.table("rawdata/Core_NonCore_Miles16.txt", h=F)
info10 <- read.table("rawdata/centromere_10.txt", h=F)
infogenome <- rbind(infogenome, info10)
chromlength <- read.table("../read/rawdata/chrom-sizes.tsv")
colnames(chromlength) <- c("ID", "length")
chromlength$chr <- sub("[^_]+_", "", chromlength$ID)
chromlength$chr <- sub("_\\S+", "", chromlength$chr)
chromlength$chr <- as.numeric(chromlength$chr)
chromlength <- chromlength[!is.na(chromlength$chr),]
chrombarcodes <- read_delim(file = "../hmmibd/out/IBD-WGS-snps.hmm.txt.zip",
                            delim = "\t")
chrombarcodes <- merge(chrombarcodes, barcodes[,c("Internal.Sample.ID", "ID")],
                       by.x = "sample1", by.y = "Internal.Sample.ID")
chrombarcodes <- merge(chrombarcodes, barcodes[,c("Internal.Sample.ID", "ID")],
                       by.x = "sample2", by.y = "Internal.Sample.ID",
                       suffixes = c("1", "2"))
chrombarcodes <- chrombarcodes[!is.na(chrombarcodes$ID1) & !is.na(chrombarcodes$ID2),]
########################################
