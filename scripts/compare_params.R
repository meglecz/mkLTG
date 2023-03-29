library(ggplot2)
library(tidyverse)
library(cowplot)
library("psych")  
library(ggExtra)
#install.packages("tidyverse")
#install.packages("psych")       
# Install psych package
#install.packages("ggExtra") 


setwd("~/mkLTG/optimize/graphs")

# read results of the variable pidentity runs
df1<-read.table(file = file.path('~/mkLTG/optimize/series1/LTG_pid_variable/stat_match.tsv'), sep = '\t', header = TRUE)
df2<-read.table(file = file.path('~/mkLTG/optimize/series2/LTG_pid_variable/stat_match.tsv'), sep = '\t', header = TRUE)
df3<-read.table(file = file.path('~/mkLTG/optimize/series3/LTG_pid_variable/stat_match.tsv'), sep = '\t', header = TRUE)
df4<-read.table(file = file.path('~/mkLTG/optimize/series4/LTG_pid_variable/stat_match.tsv'), sep = '\t', header = TRUE)
df5<-read.table(file = file.path('~/mkLTG/optimize/series5/LTG_pid_variable/stat_match.tsv'), sep = '\t', header = TRUE)

# read results of the fixed pidentity runs
dff1<-read.table(file = file.path('~/mkLTG/optimize/series1/LTG_pid_fix/stat_match.tsv'), sep = '\t', header = TRUE)
dff2<-read.table(file = file.path('~/mkLTG/optimize/series2/LTG_pid_fix/stat_match.tsv'), sep = '\t', header = TRUE)
dff3<-read.table(file = file.path('~/mkLTG/optimize/series3/LTG_pid_fix/stat_match.tsv'), sep = '\t', header = TRUE)
dff4<-read.table(file = file.path('~/mkLTG/optimize/series4/LTG_pid_fix/stat_match.tsv'), sep = '\t', header = TRUE)
dff5<-read.table(file = file.path('~/mkLTG/optimize/series5/LTG_pid_fix/stat_match.tsv'), sep = '\t', header = TRUE)

# rename variables df
varnames <-  c("param_setting","pid","pcov","phit","taxn","seqn","refres","ltgres", "taxlevel","taxlevel_index","TP","FP","FN")
names(df1) <- varnames
names(df2) <- varnames
names(df3) <- varnames
names(df4) <- varnames
names(df5) <- varnames


# rename variables dff
names(dff1) <- varnames
names(dff2) <- varnames
names(dff3) <- varnames
names(dff4) <- varnames
names(dff5) <- varnames


# make sure to have the same order in each df
df1 <- df1 %>% arrange(param_setting, taxlevel_index)
df2 <- df2 %>% arrange(param_setting, taxlevel_index)
df3 <- df3 %>% arrange(param_setting, taxlevel_index)
df4 <- df4 %>% arrange(param_setting, taxlevel_index)
df5 <- df5 %>% arrange(param_setting, taxlevel_index)


# make sure to have the same order in each dff
dff1 <- dff1 %>% arrange(param_setting, taxlevel_index)
dff2 <- dff2 %>% arrange(param_setting, taxlevel_index)
dff3 <- dff3 %>% arrange(param_setting, taxlevel_index)
dff4 <- dff4 %>% arrange(param_setting, taxlevel_index)
dff5 <- dff5 %>% arrange(param_setting, taxlevel_index)

# make nex df
df <- df1
# add TP, FP, FN columns of the 5 dataframes
df[,11:13] <- df[,11:13] + df2[,11:13]
df[,11:13] <- df[,11:13] + df3[,11:13]
df[,11:13] <- df[,11:13] + df4[,11:13]
df[,11:13] <- df[,11:13] + df5[,11:13]


# pool 5 dff
dff <- dff1
dff[,11:13] <- dff[,11:13] + dff2[,11:13]
dff[,11:13] <- dff[,11:13] + dff3[,11:13]
dff[,11:13] <- dff[,11:13] + dff4[,11:13]
dff[,11:13] <- dff[,11:13] + dff5[,11:13]

# shorten parmeter settings when possible
df$pcov[df$pcov == "100;100;100;100;100;100"] <- 100
df$pcov[df$pcov == "90;90;90;90;90;90"] <- 90
df$pcov[df$pcov == "80;80;80;80;80;80"] <- 80
df$pcov[df$pcov == "70;70;70;70;70;70"] <- 70

df$phit[df$phit == "100;100;100;100;100;100"] <- 100
df$phit[df$phit == "90;90;90;90;90;90"] <- 90
df$phit[df$phit == "80;80;80;80;80;80"] <- 80
df$phit[df$phit == "70;70;70;70;70;70"] <- 70
# replace ; by . to simplify legend of the graphs
df$pid <- gsub(";", ".", df$pid)
df$refres <- gsub(";", ".", df$refres)
df$ltgres <- gsub(";", ".", df$ltgres)
df$taxn <- gsub(";", ".", df$taxn)
df$seqn <- gsub(";", ".", df$seqn)

# TPR, PPV, F1 => this should be maximized
df$TPR = df$TP/(df$TP + df$FN)
df$PPV = df$TP/(df$TP + df$FP)
df$F1 = 2 * df$TPR * df$PPV /(df$TPR + df$PPV)
df$F1[df$F1 == "NaN"] <- 0
df$PPV[df$PPV == "NaN"] <- 0
df$TPR[df$TPR == "NaN"] <- 0


# add FPFN, TPR, PPV, F1 => this should be maximized
dff$TPR = dff$TP/(dff$TP + dff$FN)
dff$PPV = dff$TP/(dff$TP + dff$FP)
dff$F1 = 2 * dff$TPR * dff$PPV /(dff$TPR + dff$PPV)
dff$F1[dff$F1 == "NaN"] <- 0
dff$PPV[dff$PPV == "NaN"] <- 0
dff$TPR[dff$TPR == "NaN"] <- 0

#############################################################################
# make an overall F1 for each parameter setting in df
# empty df with correct dimension and colnames
df_pool <- data.frame(matrix(ncol = 10, nrow = max(df$param_setting)))
colnames(df_pool) <-  c("param_setting","pid","pcov","phit","taxn","seqn","refres","ltgres", "F1", "TP")

# make mean F1 over species, genus, family and order levels
for(i in seq(1:max(df$param_setting))){ # for each param setting
  #  print(i)
  a <- c() # make a vector with the F1 of the 4 highest resolution taxlevel
  b <- c() # make a vector with the TP of the 4 highest resolution taxlevel
  a[1] <- df[ which(df$param_setting == i & df$taxlevel_index == 1), "F1"]
  a[2] <- df[ which(df$param_setting == i & df$taxlevel_index == 2), "F1"]
  a[3] <- df[ which(df$param_setting == i & df$taxlevel_index == 3), "F1"]
  a[4] <- df[ which(df$param_setting == i & df$taxlevel_index == 4), "F1"]
  a[5] <- df[ which(df$param_setting == i & df$taxlevel_index == 5), "F1"]
  a[6] <- df[ which(df$param_setting == i & df$taxlevel_index == 6), "F1"]
  a[7] <- df[ which(df$param_setting == i & df$taxlevel_index == 7), "F1"]
  a[8] <- df[ which(df$param_setting == i & df$taxlevel_index == 8), "F1"]
  b[1] <- df[ which(df$param_setting == i & df$taxlevel_index == 1), "TP"]
  b[2] <- df[ which(df$param_setting == i & df$taxlevel_index == 2), "TP"]
  b[3] <- df[ which(df$param_setting == i & df$taxlevel_index == 3), "TP"]
  b[4] <- df[ which(df$param_setting == i & df$taxlevel_index == 4), "TP"]
  b[5] <- df[ which(df$param_setting == i & df$taxlevel_index == 5), "TP"]
  b[6] <- df[ which(df$param_setting == i & df$taxlevel_index == 6), "TP"]
  b[7] <- df[ which(df$param_setting == i & df$taxlevel_index == 7), "TP"]
  b[8] <- df[ which(df$param_setting == i & df$taxlevel_index == 8), "TP"]  
  hm <- harmonic.mean(a)
  TP <- sum(b)/(5000*8)
  df_pool[i,] <-  df[ which(df$param_setting == i & df$taxlevel_index == 8), 1:8]
  df_pool[i, "F1"] <- hm
  df_pool[i, "TP"] <- TP
}

max_F1 <- max(df_pool$F1)
print( df_pool[ which(df_pool$F1 == max_F1), ], )
print( df[ which(df$param_setting == 21186), ], )

df_pool <- df_pool %>% arrange(-F1, )
write.csv(df_pool, file="parameter_settings_variable_pid_F1.tsv")
write.csv(df[ which(df$param_setting == 21186), ], file="best_variable_pid.csv")

#############################################################################
# make an overall F1 for each parameter setting in dff
# empty dff with correct dimension and colnames
dff_pool <- data.frame(matrix(ncol = 10, nrow = max(dff$param_setting)))
colnames(dff_pool) <-  c("param_setting","pid","pcov","phit","taxn","seqn","refres","ltgres", "F1", "TP")

# make mean F1 over species, genus, family and order levels
for(i in seq(1:max(dff$param_setting))){ # for each param setting
  a <- c() # make a vector with the F1 of the 4 highest resolution taxlevel
  b <- c() # make a vector with the FPFN of the 4 highest resolution taxlevel
  a[1] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 1), "F1"]
  a[2] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 2), "F1"]
  a[3] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 3), "F1"]
  a[4] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 4), "F1"]
  a[5] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 5), "F1"]
  a[6] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 6), "F1"]
  a[7] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 7), "F1"]
  a[8] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 8), "F1"]
  b[1] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 1), "TP"]
  b[2] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 2), "TP"]
  b[3] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 3), "TP"]
  b[4] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 4), "TP"]
  b[5] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 5), "TP"]
  b[6] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 6), "TP"]
  b[7] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 7), "TP"]
  b[8] <- dff[ which(dff$param_setting == i & dff$taxlevel_index == 8), "TP"] 
  hm <- harmonic.mean(a)
  TP <- sum(b)/(5000*8)
  dff_pool[i,] <-  dff[ which(dff$param_setting == i & dff$taxlevel_index == 8), 1:8]
  dff_pool[i, "F1"] <- hm
  dff_pool[i, "TP"] <- TP
}

max_dff_F1 <- max(dff_pool$F1)
print( dff_pool[ which(dff_pool$F1 == max_dff_F1), ], )
print( dff[ which(df$param_setting == 776), ], )

dff_pool <- dff_pool %>% arrange(-F1, )
write.csv(dff_pool, file="parameter_settings_fixed_pid_F1.tsv")
write.csv(dff[ which(dff$param_setting == 776), ], file="best_fixed_pid.csv")


#########################################################
#Barplots by taxlevel and parameters

# order factors
df_graph <- df
df_graph$taxlevel <- factor(df_graph$taxlevel , levels=c("species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))
df_graph$pcov <- factor(df_graph$pcov , levels=c(100, 90, 80, 70))
df_graph$phit <- factor(df_graph$phit , levels=c(100, 90, 80, 70))

# plots for all tax level F1
p<- ggplot(data=df_graph, aes(x=pcov, y=F1)) +
  geom_col() +
  facet_wrap(~taxlevel, nrow=4)
p = p + labs(title = "F1 sumed over all paramater settings \nfor each %of coverage (pcov)")
p = p + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5))
p = p + labs(x="pcov", y="F1")
p = p + theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
p = p + theme(axis.text.y = element_text(size = 8))
p = p + theme(panel.grid.major = element_line(color = "grey", size = 0.2,linetype = 2))
p = p + theme(panel.border = element_rect(color = "grey",fill = NA, size = 0.5))
p

ggsave("pcov_F1_all.png", width = 119, height = 100, units = "mm", dpi = 300)
dev.off()
ggsave("pcov_F1_all.eps", width = 119, height = 100, units = "mm", dpi = 300)
dev.off()


p<- ggplot(data=df_graph, aes(x=phit, y=F1)) +
    geom_col() +
  facet_wrap(~taxlevel, nrow=4)
p = p + labs(title = "F1 sumed over all paramater settings \nfor each % of hits (phit)")
p = p + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5))
p = p + labs(x="phit", y="F1")
p = p + theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
p = p + theme(axis.text.y = element_text(size = 8))
p = p + theme(panel.grid.major = element_line(color = "grey", size = 0.2,linetype = 2))
p = p + theme(panel.border = element_rect(color = "grey",fill = NA, size = 0.5))
p

ggsave("phit_F1_all.png", width = 119, height = 100, units = "mm", dpi = 300)
dev.off()
ggsave("phit_F1_all.eps", width = 119, height = 100, units = "mm", dpi = 300)
dev.off()


p<- ggplot(data=df_graph, aes(x=taxn, y=F1)) +
  geom_col() +
  facet_wrap(~taxlevel, nrow=4)
p = p + labs(title = "F1 sumed over all paramater settings \nfor each minimum number of taxa (taxn)")
p = p + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5))
p = p + labs(x="Minimum number of taxa", y="F1")
p = p + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1))
p = p + theme(axis.text.y = element_text(size = 8))
p = p + theme(panel.grid.major = element_line(color = "grey", size = 0.2,linetype = 2))
p = p + theme(panel.border = element_rect(color = "grey",fill = NA, size = 0.5))
p

ggsave("taxn_F1_all.png", width = 195, height = 100, units = "mm", dpi = 300)
dev.off()
ggsave("taxn_F1_all.eps", width = 195, height = 100, units = "mm", dpi = 300)
dev.off()


p<- ggplot(data=df_graph, aes(x=refres, y=F1)) +
  geom_col() +
  facet_wrap(~taxlevel, nrow=4)
p = p + labs(title = "F1 sumed over all paramater settings \nfor each minimum resolution of reference sequences (refres)")
p = p + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5))
p = p + labs(x="Minimum resolution of reference sequences", y="F1")
p = p + theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1))
p = p + theme(axis.text.y = element_text(size = 8))
p = p + theme(panel.grid.major = element_line(color = "grey", size = 0.2,linetype = 2))
p = p + theme(panel.border = element_rect(color = "grey",fill = NA, size = 0.5))
p

ggsave("refres_F1_all.png", width = 195, height = 100, units = "mm", dpi = 300)
dev.off()
ggsave("refres_F1_all.eps", width = 195, height = 100, units = "mm", dpi = 300)
dev.off()






