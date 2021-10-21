## CCI analysis
library(circlize)
library(tidyverse)

CCI_HtoM <- read.csv("STCCI_HtoM.csv")
CCI_MtoH <- read.csv("STCCI_MtoH.csv")


Chemo = CCI_HtoM
Chemo2 = CCI_MtoH
uni.chemo = list(Chemo[,1], Chemo[,2], Chemo2[,1], Chemo2[,2])
uni.chemo = unique(unlist(uni.chemo))
CCI_h = ST.combined
CCI_m = mST.combined
CCI_h = CCI_h[uni.chemo,]
CCI_m = CCI_m[uni.chemo,]

unigenes = rownames(CCI_h)
mel_h = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCI_h, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCI_h[, which(x = expr > 2)]
		ge = table(ge@meta.data$group)
		ho = reshape2::melt(ge)
		ho = dplyr::mutate(ho, Gene = unigenes[i])
		mel_h = rbind(mel_h, ho)

	}else{
		print(unigenes[i])
	}
	
}
unigenes = rownames(CCI_m)
mel_m = list()
for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCI_m, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCI_m[, which(x = expr > 2)]
		ge = table(ge@meta.data$group)
		ho = reshape2::melt(ge)
		ho = dplyr::mutate(ho, Gene = unigenes[i])
		mel_m = rbind(mel_m, ho)

	}else{
		print(unigenes[i])
	}
	
}
mel_h <- mel_h[mel_h$value > 0, ]
melmh <- mel_m[mel_m$value > 0, ]


tocell = table(CCI_h@meta.data$group)
Ratiomel_h = list()

for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel_h[mel_h[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel_h = rbind(Ratiomel_h, LS)

}
Ratiome_h = reshape2::dcast(Ratiomel_h, Var1 ~ Gene, value.var = "value")


tocell = table(CCI_m@meta.data$group)
Ratiomel_m = list()
for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel_m[mel_m[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel_m = rbind(Ratiomel_m, LS)

}
Ratiomel_m = reshape2::dcast(Ratiomel_m, Var1 ~ Gene, value.var = "value")


write.table(Ratiomel_h, file = "CCI-hRatioMAT.txt", sep = "	")
write.table(Ratiomel_m, file = "CCI-mRatioMAT.txt", sep = "	")


A = list()
B = list()
CCItotal = list()
# HtoM
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel_h[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel_m[,3])),]
CCI_perlist = dplyr::mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
HtoM = as.character(unique(CCI_perlist$merge))
aa = str_split(HtoM, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]


# MtoH
CCI_perlist2 = Chemo2[Chemo2$AliasA %in% unlist(unique(mel_m[,3])),]
CCI_perlist2 = CCI_perlist2[CCI_perlist2$AliasB %in% unlist(unique(mel_h[,3])),]
CCI_perlist2 = dplyr::mutate(CCI_perlist2, merge = paste(CCI_perlist2$AliasA, CCI_perlist2$AliasB, sep = "&"))
MtoH = as.character(unique(CCI_perlist2$merge))
aa2 = str_split(MtoH, "&", n = 2, simplify = TRUE)
CCI_perlist2 = data.frame(L = aa2[,1], R = aa2[,2], merge = unique(CCI_perlist2$merge))
CCI_perlist2 = CCI_perlist2[!(as.character(CCI_perlist2$L) == as.character(CCI_perlist2$R)),]


Hratiomelt = na.omit(reshape2::melt(Ratiomel_h))
Hratiomelt = Hratiomelt[Hratiomelt$value > 0,]
Mratiomelt = na.omit(reshape2::melt(Ratiomel_m))
Mratiomelt = Mratiomelt[Mratiomelt$value > 0,]


# HtoM
gr = unique(Hratiomelt$Var1)
combihtom = list()
for(i in 1:length(gr)){
	g <- gr[i]
	hogeH <- Hratiomelt[Hratiomelt$Var1 == g,]
	hogeM <- Mratiomelt[Mratiomelt$Var1 == g,]

	geneLlist = unique(hogeH$Gene)
	geneRlist = unique(hogeM$variable)


	for(k in 1: length(geneLlist)){
		geneL = geneLlist[k]
		combi = CCI_perlist[CCI_perlist$L == geneL,]
		combi = combi[combi$R %in% geneRlist,]
		Lscore = hogeH[hogeH$Gene == geneL,]$value
		Rscore = hogeM[hogeM$variable %in% unique(combi$R),]$value
		combi2 = dplyr::mutate(combi, group = g, Lscore = Lscore)
		combi2 = cbind(combi2, Rscore)
		CCIscore = combi2$Lscore * combi2$Rscore
		combi2 = cbind(combi2, CCIscore)
		combihtom = rbind(combihtom, combi2)
	}	
}
# MtoH
gr = unique(Mratiomelt$Var1)
combimtoh = list()
for(i in 1:length(gr)){
	g <- gr[i]
	hogeH <- Hratiomelt[Hratiomelt$Var1 == g,]
	hogeM <- Mratiomelt[Mratiomelt$Var1 == g,]

	geneLlist = unique(hogeM$variable)
	geneRlist = unique(hogeH$Gene)


	for(k in 1: length(geneLlist)){
		geneL = geneLlist[k]
		combi = CCI_perlist2[CCI_perlist2$L == as.character(geneL),]
		combi = combi[combi$R %in% geneRlist,]
		Lscore = hogeM[hogeM$variable == as.character(geneL),]$value
		Rscore = hogeH[hogeH$Gene %in% unique(combi$R),]$value
		combi2 = dplyr::mutate(combi, group = g, Lscore = Lscore)
		combi2 = cbind(combi2, Rscore)
		CCIscore = combi2$Lscore * combi2$Rscore
		combi2 = cbind(combi2, CCIscore)
		combimtoh = rbind(combimtoh, combi2)
	}	
}
write.table(combihtom, file = "CCIhtom.csv", sep = ",")
write.table(combimtoh, file = "CCImtoh.csv", sep = ",")



hHtoM = combihtom[combihtom$CCIscore > 0.1,]
hMtoH = combimtoh[combimtoh$CCIscore > 0.1,]
hHtoMmat = reshape2::dcast(hHtoM, merge ~ group, value.var = "CCIscore")
hMtoHmat = reshape2::dcast(hMtoH, merge ~ group, value.var = "CCIscore")

write.table(hHtoM, file = "CCI0.1htom.csv", sep = ",")
write.table(hMtoH, file = "CCI0.1mtoh.csv", sep = ",")
write.table(table(hHtoM$merge,hHtoM$group), file = "CCI0.1htomTable.csv", sep = ",")
write.table(table(hMtoH$merge,hMtoH$group), file = "CCI0.1mtohTable.csv", sep = ",")
write.table(hHtoMmat, file = "CCImatrix-hHtoM.csv", sep = ",")
write.table(hMtoHmat, file = "CCImatrix-hMtoH.csv", sep = ",")

hHtoMmat[is.na(hHtoMmat)] <- 0
hMtoHmat[is.na(hMtoHmat)] <- 0


### bubble plot HtoM
topHtoM <- combihtom[combihtom$CCIscore > 0.3,]

ggplot(topHtoM,aes(x = group, y = merge, size = CCIscore))+
  geom_point(aes(size = CCIscore), shape = 21, colour = "black", fill = "red") +
  scale_size_area(max_size = 20) +
  theme_bw()


### bubble plot MtoH
topMtoH <- combimtoh[combimtoh$CCIscore > 0.3,]

ggplot(topMtoH,aes(x = group, y = merge, size = CCIscore))+
  geom_point(aes(size = CCIscore), shape = 21, colour = "black", fill = "red") +
  scale_size_area(max_size = 20) +
  theme_bw()


###############
###############
###############
### circle bar plot
library(tidyverse)
 
# Create dataset
data <- hHtoM[,c(3,4,7)]
data <- data[data$CCIscore > 0.3,]
data <- data[order(-data$CCIscore),]

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))
 
# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
 
# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=CCIscore, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(0,1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=CCIscore, label=merge, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3.5, angle= label_data$angle, inherit.aes = FALSE ) 
 
p

