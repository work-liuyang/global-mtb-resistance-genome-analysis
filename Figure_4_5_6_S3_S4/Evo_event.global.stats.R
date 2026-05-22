library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(corrplot)
library(ggsci)
library(reshape2)


color<-c('#5050FF', '#CE3D32', '#749B58', '#F0E685', 
         '#466983', '#BA6338', '#5DB1DD', '#802268', 
         '#6BD76B', '#D595A7', '#924822', '#837B8D',
         '#C75127', '#D58F5C', '#7A65A5', '#E4AF69',
         '#3B1B53', '#CDDEB7', '#612A79', '#AE1F63',
         '#E7C76F', '#5A655E', '#CC9900', '#99CC00',
         '#A9A9A9', '#33CC00', '#00CC33', '#00CC99', 
         '#0099CC', '#0A47FF', '#4775FF', '#FFC20A', 
         '#FFD147', '#990033', '#991A00', '#996600', 
         '#809900', '#339900', '#00991A', '#009966', 
         '#008099', '#003399', '#1A0099', '#660099', 
         '#990080', '#D60047', '#FF1463', '#00D68F', 
         '#14FFB1' )

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

resacdata<-read.table("global.res_ac.stats.xls", 
                      sep="\t",header = T,check.names = F)
country <- read.table("country.list",
                      sep="\t",header = F,check.names = F)
names(country) <- c('region', 'continent')
gdp <- read.table("gdp.xls",sep = "\t",header = T, 
                  row.names = 1, check.names = F)

country <- arrange(country,country$continent, country$region)
allcolor <- data.frame(country, color = color[1:dim(country)[1]])


#----------------------------------------------------------------------
#global
global<-resacdata
global$event<-factor(global$event,
                     levels = c("Sus2HR","Sus2RR","Sus2MDR","Sus2Pre-XDR",
                                "Sus2XDR","HR2MDR","HR2Pre-XDR","RR2MDR",
                                "RR2Pre-XDR","MDR2Pre-XDR","MDR2XDR","Pre-XDR2XDR"))
#freq <- table(global$event)
#filter <- names(freq[freq < 3])
create_dir("general_stats")
p1 <- ggboxplot(global, x = "event", y = "age", fill = "#9EADD4", legend="none") + 
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1)) +
  labs(y="Date of resistance acquisition(years ago)",x=NULL, 
       title = "The evolution of events leading to MTB acquiring resistance") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(method = "kruskal.test",aes(label = "p.format"),
                     label.x = 2, label.y = 160, size = 4) + 
  scale_x_discrete(labels=c('Sensitive->HR-TB','Sensitive->RR-TB',
                            'Sensitive->MDR-TB','Sensitive->Pre-XDR-TB',
                            'Sensitive->XDR-TB','HR-TB->MDR-TB',
                            'HR-TB->Pre-XDR-TB','RR-TB->MDR-TB',
                            'RR-TB->Pre-XDR-TB','MDR-TB->Pre-XDR-TB',
                            'MDR-TB->XDR-TB','Pre-XDR-TB->XDR-TB'))


ggsave(filename = "general_stats/global.evo_event.stats.png",plot = p1,width = 7.68, height = 5.7,dpi = 300)
ggsave(filename = "general_stats/global.evo_event.stats.pdf",plot = p1,width = 7.68, height = 5.7)

global.stat <- global %>% pairwise_wilcox_test (age ~ event, p.adjust.method = "bonf")
global.sig <- global.stat[global.stat$p.adj<0.05,
                 c('group1','group2','n1','n2','p','p.adj','p.adj.signif')]
write.table(global.sig,file = "general_stats/global.pairwise.sig.xls",quote = F, sep = "\t",
            row.names = F, col.names = T)
write.table(global.stat[,c('group1','group2','n1','n2','p','p.adj','p.adj.signif')],
            file = "general_stats/global.pairwise.all.xls",quote = F, sep = "\t",
            row.names = F, col.names = T)

#----------------------------------------------------------------------
#drug
drug <- resacdata
for (i in 1:dim(resacdata)[1]) {
  if(resacdata$event[i] %in% c("Sus2HR","RR2MDR")){
    drug[i,'drug'] = 'Isoniazid'
  }
  else if(resacdata$event[i] %in% c("Sus2RR","HR2MDR")){
    drug[i,'drug'] = 'Rifampicin'
  }
  else if(resacdata$event[i] %in% c("MDR2Pre-XDR")){
    drug[i,'drug'] = 'Fluoroquinolones'
  }
  else if(resacdata$event[i] %in% c("Pre-XDR2XDR")){
    drug[i,'drug'] = 'Bedaquiline/Linezolid'
  }
  else{
    drug[i,'drug'] = 'Complex'
  }
}
for (i in 1:dim(drug)[1]) {
  if(drug$drug[i]=="Isoniazid"){
    drug[i,'drugdate'] <- 72
  }
  else if(drug$drug[i]=="Rifampicin"){
    drug[i,'drugdate'] <- 59
  }
  else if(drug$drug[i]=="Fluoroquinolones"){
    drug[i,'drugdate'] <- 39
  }
  else if(drug$drug[i]=="Complex"){
    drug[i,'drugdate'] <- 'NA'
  }
  else{
    drug[i,'drugdate'] <- 24
  }
}

write.table(drug,"general_stats/event2drug.stats.xls",sep = "\t",quote = F,row.names = F,col.names = T)

#drug vs GDP
create_dir("drug_vs_gdp")

gdp <- read.table("gdp.xls", 
                  sep="\t",header = T,check.names = F, row.names = 1)

drug_gdp <- function(data, drug, abbr, width, height){
  annox <- 35000
  if(abbr=="FLQ") annox <- 7000
  drugdf <- data[data$drug==drug,] %>% arrange(region)
  freq <- table(drugdf$region)
  filter <- names(freq[freq < 5])
  drugdf <- drugdf[!(drugdf$region %in% filter),]
  median_date <- aggregate(drugdf$age,by=list(drugdf$region),median) %>% arrange(Group.1)
  names(median_date) <- c('region','medianage')
  median_date <- cbind(median_date,gdp[median_date$region,])
  #mod <- summary(lm(medianage~median,median_date %>% arrange(median)))
  mod <- summary(lm(medianage~mean,median_date %>% arrange(mean)))
  pvalue <- round(mod$coefficients[2,4],2)
  R <- round(mod$r.squared,2)
  plm <- ggplot(median_date %>% arrange(median),aes(median,medianage)) + 
    geom_point(size=3,color="#9EADD4") + 
    geom_smooth(method = "lm", se = FALSE, color='black') + 
    labs(title=paste0("Median time of resistance acquisition of ",abbr," vs median GDP per capita"),
         x="Median GDP per capita (US$)", y='Median time of resistance acquisition(years ago)') + theme_classic() + 
    theme(axis.text = element_text(size=9),
          plot.title = element_text(size=9,hjust = 0.5),
          axis.title = element_text(size=9)) + 
    annotate('text',x = annox,y=10,label=paste0("R^2=",R,", F-test p=",pvalue), size=4)

  ggsave(paste0("drug_vs_gdp/",drug,"_vs_gdp.png"),plm,width = width,height = height,dpi = 300)
  ggsave(paste0("drug_vs_gdp/",drug,"_vs_gdp.pdf"),plm,width = width,height = height)
}

drug_gdp(drug, 'Isoniazid', 'INH', 10, 7)
drug_gdp(drug, 'Rifampicin', 'RIF', 10, 7)
drug_gdp(drug, 'Fluoroquinolones', 'FLQ', 10, 7)

#drug vs date introduced
coredrug <- drug[drug$drug!="Complex",]
coredrug$drugdate<-as.numeric(coredrug$drugdate)
coredrug$drug<-factor(coredrug$drug,levels = unique(coredrug$drug))

p2 <- ggboxplot(coredrug, x = "drug", y = "age", fill = "drug", width = 0.5,
                bxp.errorbar=TRUE,bxp.errorbar.width=0.1,outlier.shape = NA) + 
  theme(axis.text = element_text(size=9), axis.title.y = element_text(size=9)) +
  labs(y="Date of resistance acquisition(years ago)",x=NULL, 
       title = NULL) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none') +
  scale_fill_npg()
p2 <- p2 + stat_compare_means(method = "kruskal.test",aes(label = "p.format"),
                              label.x = 1, label.y = 160, size = 3)
drug.test <- coredrug %>% pairwise_wilcox_test (age ~ drug, p.adjust.method = "bonf")
#drug.test<-drug.test%>%add_y_position()
#drug.test$y.position<-drug.test$y.position-100

p2 <- p2 + stat_pvalue_manual(drug.test,label = "p.adj.signif", y.position = 90 * seq(1.1, 1.5, length.out=3),
                              tip.length = 0.01,hide.ns = T)

ggsave(filename = "general_stats/global.drug.stats.png",plot = p2,width = 7.68, height = 5.7,dpi = 300)
ggsave(filename = "general_stats/global.drug.stats.pdf",plot = p2,width = 7.68, height = 5.7)


#coredrug <- coredrug %>%
#  group_by(drug) %>%
#  mutate(
#    group_mean = mean(age),
#    group_sd = sd(age),
#    group_upper_limit = group_mean + 2 * group_sd,
#    group_lower_limit = group_mean - 2 * group_sd,
    # 标记是否为组内极端值
#    is_group_outlier = ifelse(age > group_upper_limit | age < group_lower_limit, 
#                              "outlier", "normal")
#  ) %>% ungroup()

#coredrug_clean_by_group <- coredrug %>% filter(is_group_outlier == "normal")

drug_summary <- coredrug %>%
group_by(drug, drugdate) %>%
summarise(median_resistance_time = median(age),
.groups = 'drop')
drug_summary
summary_model <- lm(median_resistance_time ~ drugdate, data = drug_summary)
summary(summary_model)


#-----------------------------------------------------------------
#all_event by continent
create_dir("by_continent")

by_continent <- resacdata %>% arrange(continent)
by_continent$event<-factor(by_continent$event,
                           levels = c("Sus2HR","Sus2RR","Sus2MDR","Sus2Pre-XDR",
                                      "Sus2XDR","HR2MDR","HR2Pre-XDR","RR2MDR",
                                      "RR2Pre-XDR","MDR2Pre-XDR","MDR2XDR","Pre-XDR2XDR"))

p3 <- ggboxplot(by_continent, x = "continent", y = "age", fill = "event", 
                width = 0.5,outlier.shape = 19, alpha=1) + 
  theme(axis.text.x = element_text(size=12)) +
  labs(y="Date of resistance acquisition(years ago)",x='Continent', 
       title = NULL, fill='') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  scale_fill_brewer(palette = "Set3",
                    labels=c('Sensitive->HR-TB','Sensitive->RR-TB',
                             'Sensitive->MDR-TB','Sensitive->Pre-XDR-TB',
                             'Sensitive->XDR-TB','HR-TB->MDR-TB',
                             'HR-TB->Pre-XDR-TB','RR-TB->MDR-TB',
                             'RR-TB->Pre-XDR-TB','MDR-TB->Pre-XDR-TB',
                             'MDR-TB->XDR-TB','Pre-XDR-TB->XDR-TB')) +
  guides(fill  = guide_legend(ncol = 3))
#p3 <- p3 + stat_compare_means(method = "kruskal.test",aes(label = "p.format"),
#                              label.x = 1, label.y = 160, size = 4)
by_continent.test <- by_continent %>% pairwise_wilcox_test (age ~ continent, p.adjust.method = "bonf")


ggsave(filename = "general_stats/by_continent.evo_event.stats.png",plot = p3,width = 7.68, height = 5.7,dpi = 300)
ggsave(filename = "general_stats/by_continent.evo_event.stats.pdf",plot = p3,width = 7.68, height = 5.7)

by_continent.sig <- by_continent.test[by_continent.test$p.adj<0.05,
                          c('group1','group2','n1','n2','p','p.adj','p.adj.signif')]
write.table(by_continent.sig,file = "by_continent/by_continent.pairwise.sig.xls",quote = F, sep = "\t",
            row.names = F, col.names = T)
write.table(by_continent.test[,c('group1','group2','n1','n2','p','p.adj','p.adj.signif')],
            file = "by_continent/by_continent.pairwise.all.xls",quote = F, sep = "\t",
            row.names = F, col.names = T)


#-----------------------------------------------------------------
#Distribution of date of per drug in each contitent and country/region
stat_by_continent <- function(x, prefix, introdate, delta, width, height){
  
  mycolor <- data.frame(
    row.names = c("Africa","Asia","Europe","North America","Oceania","South America"),
    color = c("#3B4992","#EE0000","#008B45","#631879","#008280","#BB0021")
  )

  labelx <- 6
  if(prefix == "Fluoroquinolones"){
    labelx <- 5
  }
  px <- ggboxplot(x, x = 'continent', y = "age", fill = 'continent', 
                  width = 0.9,outlier.shape = 19, alpha=0.6) + 
    #geom_hline(yintercept = introdate, linetype = "dashed") +
    theme(axis.text.x = element_text(size=12,angle = 45,hjust = 1)) +
    labs(y="Time of resistance acquisition(years ago)",x=NULL, 
         title = prefix) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none') + 
    scale_fill_manual(values = mycolor[unique(x$continent),]) +
    scale_y_continuous(limits = c(0,150))
  
  px <- px + stat_compare_means(method = "kruskal.test",aes(label = "p.format"),
                                label.x = labelx, label.y = delta, size = 4) + coord_flip()
  stat.test <- x %>% pairwise_wilcox_test (age ~ continent, p.adjust.method = "bonf")
  stat.test<- stat.test%>%add_y_position()
  stat.test$y.position<- stat.test$y.position-delta
  
  #px <- px + stat_pvalue_manual( stat.test,label = "p.adj.signif", 
  #                              tip.length = 0.01,hide.ns = T)
  
  ggsave(filename = paste0("by_continent/",prefix,".by_continent.stats.png"),plot = px,width = width, height = height,dpi = 300)
  ggsave(filename = paste0("by_continent/",prefix,".by_continent.stats.pdf"),plot = px,width = width, height = height)
  
  if(prefix == "Fluoroquinolones"){
    return()
  }
  
  diag <- data.frame()
  tmp <- stat.test[,c("group1","group2","p.adj")]
  tmp2 <- unique(c(unique(stat.test$group1),unique(stat.test$group2)))
  for (i in 1:length(tmp2)) {
    diag[i,'group1'] <- tmp2[i]
    diag[i,'group2'] <- tmp2[i]
    diag[i,'p.adj'] <- 1
  }
  md <- acast(rbind(tmp,diag), group1~group2, value.var = "p.adj", fill=0)
  md[lower.tri(md)] <- t(md)[lower.tri(md)]
  
  mdsize <- md
  mdsize[mdsize>=0.05] <- 1
  
  
  pdf(paste0("by_continent/",prefix,".by_continent.heatmap.pdf"), width=width, height=height)
  corrplot(mdsize,is.corr = F,type='upper', method = 'color',
           col = colorRampPalette(c('#9BD5E7','white'))(200),cl.pos = 'n', 
           p.mat = md,insig = 'label_sig',sig.level = c(.001, .01, .05),
           pch.cex = 2,tl.pos="lt", tl.cex=1.3, tl.col="black",
           tl.srt = 45,tl.offset=0.5,addgrid.col = 'black')
  corrplot(mdsize,is.corr = F,type='lower', method = 'number',p.mat = md,
           pch.cex = 1.5,tl.pos="n",add = TRUE, number.cex = 1, number.digits = 5,
           col = colorRampPalette(c('#9BD5E7','white'))(200),
           insig = 'pch',cl.pos = 'n',diag = FALSE)
  
  dev.off() 
  
}

create_dir("by_country")

stat_by_country <- function(x, prefix, introdate, delta, width, height){
  
  freq <- table(x$region)
  filter <- names(freq[freq < 5])
  filter_x <- x[!(x$region %in% filter),]
  for (i in 1:dim(x)[1]) {
    if (x$region[i] %in% filter) {
      x[i,'number'] <- "lt5"
    }else{
      x[i,'number'] <- "ge5"
    }
  }
  
  px <- ggboxplot(x, x = 'region', y = "age", fill = "number", 
                  width = 0.5,outlier.shape = 19) + 
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1)) +
    labs(y="Date of resistance acquisition(years ago)",x=NULL, 
         title = prefix) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = introdate, linetype = "dashed") +
    scale_fill_manual(values=c("grey", "#9EADD4"),
                      breaks=c("lt5","ge5"),
                      labels=c("n<5", "n>=5"),
                      name="")
  
  px <- px + stat_compare_means(method = "kruskal.test",aes(label = "p.format"),
                                label.x = 2, label.y = delta, size = 4)
  
  stat.test <- filter_x %>% pairwise_wilcox_test (age ~ region, p.adjust.method = "bonf")
  stat.test<- stat.test%>%add_y_position()
  stat.test$y.position<- stat.test$y.position-delta
  write.table(stat.test[,c('group1','group2','n1','n2','p','p.adj','p.adj.signif')],
              file = paste0("by_country/",prefix,"_bycountry.pairwise.all.xls"),quote = F, sep = "\t",
              row.names = F, col.names = T)
  
  #px <- px + stat_pvalue_manual( stat.test,label = "p.adj.signif", 
  #                              tip.length = 0.01,hide.ns = T)
  
  ggsave(filename = paste0("by_country/",prefix,".by_country.stats.png"),plot = px,width = 7.68, height = 5.7,dpi = 300)
  ggsave(filename = paste0("by_country/",prefix,".by_country.stats.pdf"),plot = px,width = 7.68, height = 5.7)
  
  diag <- data.frame()
  tmp <- stat.test[,c("group1","group2","p.adj")]
  tmp2 <- unique(c(unique(stat.test$group1),unique(stat.test$group2)))
  for (i in 1:length(tmp2)) {
    diag[i,'group1'] <- tmp2[i]
    diag[i,'group2'] <- tmp2[i]
    diag[i,'p.adj'] <- 1
  }
  md <- reshape2::acast(rbind(tmp,diag), group1~group2, value.var = "p.adj", fill=0)
  md[lower.tri(md)] <- t(md)[lower.tri(md)]
  
  mdsize <- md
  mdsize[mdsize>=0.05] <- 1
  
  pdf(paste0("by_country/",prefix,".by_country.heatmap.pdf"), width=width, height=height)
  corrplot(mdsize,is.corr = F,type='upper', method = 'color',
           col = colorRampPalette(c('#9BD5E7','white'))(200),cl.pos = 'n', 
           p.mat = md,insig = 'label_sig',sig.level = c(.001001, .01, .05),
           pch.cex = 2,tl.pos="lt", tl.cex=1.5, tl.col="black",
           tl.srt = 45,tl.offset=0.7,addgrid.col = 'black')
  corrplot(mdsize,is.corr = F,type='lower', method = 'number',p.mat = md,
           pch.cex = 1.5,tl.pos="n",add = TRUE, number.cex = 1, number.digits = 3,
           col = colorRampPalette(c('#9BD5E7','white'))(200),
           insig = 'pch',cl.pos = 'n',diag = FALSE)
  
  dev.off()  
  
  
}

INH <- drug[drug$drug=="Isoniazid",] %>% arrange(continent,region)
stat_by_continent(INH, "Isoniazid", 72, 100, 6, 5)
stat_by_country(INH, "Isoniazid", 72, 150, 12, 12)

RIF <- drug[drug$drug=="Rifampicin",] %>% arrange(continent,region)
stat_by_continent(RIF, "Rifampicin", 59, 100, 6, 5)
stat_by_country(RIF, "Rifampicin", 59, 150, 12, 12)

FLQ <- drug[drug$drug=="Fluoroquinolones",] %>% arrange(continent,region)
stat_by_continent(FLQ, "Fluoroquinolones", 39, 100, 6, 6)
stat_by_country(FLQ, "Fluoroquinolones", 39, 150, 12, 12)

               
