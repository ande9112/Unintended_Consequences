"""This code will take in SNP data from .vcf files that were converted to .csv files with parsed allele calls and calculate the number of homozygous SNPs"

# Written by Jean-Michel Michno


#set working directory
setwd("~/Downloads")

#use ggplot for all graphs in this code
library(ggplot2)
#read in the files
FN <- read.csv("HomoFNCombo.csv", header = TRUE)
WPT <- read.csv("HomoWPTCombo.csv", header = TRUE)
FN2 <- read.csv("FNCombo.csv", header = TRUE)
WPT2 <-read.csv("WPTCombo.csv", header = TRUE)
WPT <- as.data.frame(WPT)

#make all the the important columns variables
WPT$Bert1 <- as.character(WPT$Bert1)
WPT$BertMN2 <- as.character(WPT$Bert2)
WPT$WPT389 <- as.character(WPT$WPT389)
WPT$WPT391 <- as.character(WPT$WPT391)

#test code to see alternative ways (skip)
#HBert <- WPT[WPT$Bert1 == WPT$BertMN2 & WPT$WPT389 != WPT$WPT391,]
#HBert2 <- WPT[WPT$Bert1 == '0/0'&WPT$Bert2 == '0/0'& WPT$WPT389 != WPT$WPT391 & WPT$WPT389!='./.' & WPT$WPT391!='./.',]

#Do hard calls on SNP's and all what the sample ID is on the in a new column on the end
WTest1 <- WPT[WPT$Bert1 == '1/1' & WPT$Bert2 == '0/0'& WPT$WPT389 == '0/0'& WPT$WPT391 == '0/0',]
#make a new vector
Line <- c()
#see how many values we have and append the ID at the end of each row
for (i in 1:length(WTest1[,1])){
  add <- "BertMN01_1"
  Line <- append(Line,add)
}
#combine data and line ID together into one variable
WTest1 <- cbind(WTest1,Line)

#Repeat for all other lines
########
WTest2 <- WPT[WPT$Bert1 == '0/0' & WPT$Bert2 == '1/1' & WPT$WPT389 == '0/0' & WPT$WPT391 == '0/0',]
Line <- c()
for (i in 1:length(WTest2[,1])){
  add <- "BertMN01_2"
  Line <- append(Line,add)
}
WTest2 <- cbind(WTest2,Line)

#########
WTest3 <- WPT[WPT$Bert1 == '0/0' & WPT$Bert2 == '0/0'& WPT$WPT389 == '1/1'& WPT$WPT391 == '0/0',]
Line <- c()
for (i in 1:length(WTest3[,1])){
  add <- "WPT389_2_2"
  Line <- append(Line,add)
}
WTest3 <- cbind(WTest3,Line)

######
WTest4 <- WPT[WPT$Bert1 == '0/0' & WPT$Bert2 == '0/0'& WPT$WPT389 == '0/0'& WPT$WPT391 == '1/1',]
Line <- c()
for (i in 1:length(WTest4[,1])){
  add <- "WPT391_1_6"
  Line <- append(Line,add)
}
WTest4 <- cbind(WTest4,Line)

#Futher combine the data from above
WPTTable <- rbind(WTest1,WTest2,WTest3,WTest4)
#write it to a .csv
write.csv(WPTTable, file = "WPTVariantTable.csv")

#pull out relavent information/data for circos
CWPT <- cbind(WPTTable[,1:2],WPTTable[,2],WPTTable$Line)
row.names(CWPT) <- NULL
colnames(CWPT) <- c("chr","start","stop","id")
write.table(CWPT, file = "WPTCircos.txt", sep = " ", quote = FALSE, row.names = F)






#see how much of the original data is missing
length(which(WPT$Bert1 == './.'))
length(which(WPT$Bert2 == './.'))
length(which(WPT$WPT389 == './.'))
length(which(WPT$WPT391 == './.'))


########### This was a test on het SNP data that was ignored
#WPT2$Bert1 <- as.character(WPT2$Bert1)
#WPT2$BertMN2 <- as.character(WPT2$Bert2)
#WPT2$WPT389 <- as.character(WPT2$WPT389)
#WPT2$WPT391 <- as.character(WPT2$WPT391)

#W2Test1 <- WPT2[WPT2$Bert1 == '0/1' & WPT2$Bert2 == '0/0'& WPT2$WPT389 == '0/0'& WPT2$WPT391 == '0/0',]
#W2Test2 <- WPT2[WPT2$Bert1 == '0/0' & WPT2$Bert2 == '0/1'& WPT2$WPT389 == '0/0'& WPT2$WPT391 == '0/0',]
#W2Test3 <- WPT2[WPT2$Bert1 == '0/0' & WPT2$Bert2 == '0/0'& WPT2$WPT389 == '0/1'& WPT2$WPT391 == '0/0',]
#W2Test4 <- WPT2[WPT2$Bert1 == '0/0' & WPT2$Bert2 == '0/0'& WPT2$WPT389 == '0/0'& WPT2$WPT391 == '0/1',]


################## Original SN calling on a subset of FN data
#FN$M92 <- as.character(FN$M92)
#FN$FN03 <- as.character(FN$FN03)
#FN$FN09 <- as.character(FN$FN09)


#AFN <- FN[(FN$FN03 != FN$FN09) & (FN$M92 == '0/0')& (FN$FN03 != './.')& (FN$FN09 != './.'),]

#AFN12 <- FN[(FN$M92 == '1/1')& (FN$FN03 == '0/0') &(FN$FN09 == '0/0'),]
#AFN13<- FN[(FN$M92 == '0/0')& (FN$FN03 == '1/1') &(FN$FN09 == '0/0'),]
#AFN14 <- FN[(FN$M92 == '0/0')& (FN$FN03 == '0/0') &(FN$FN09 == '1/1'),]

#length(which(AFN$FN03 == './.'))
#length(which(AFN$FN09 == './.'))
#length(which(AFN$M92 == './.'))

#testi <- AFN[AFN$FN03 == '1/1',]

#########
#FN2$M92 <- as.character(FN2$M92)
#FN2$FN03 <- as.character(FN2$FN03)
#FN2$FN09 <- as.character(FN2$FN09)


#AFN3 <- FN2[(FN2$M92 == '0/1')& (FN2$FN03 == '0/0') &(FN2$FN09 == '0/0'),]
#AFN4 <- FN2[(FN2$M92 == '0/0')& (FN2$FN03 == '0/1') &(FN2$FN09 == '0/0'),]
#AFN5 <- FN2[(FN2$M92 == '0/0')& (FN2$FN03 == '0/0') &(FN2$FN09 == '0/1'),]


#length(which(FN2$FN03 == './.'))
#length(which(FN2$FN09 == './.'))
#length(which(FN2$M92 == './.'))

#AFN12$ID <- 1
#AFN13$ID <- 1
#AFN14$ID <- 1
#library(ggplot2) 
#ggplot(AFN12, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('M92-220')
#ggplot(AFN13, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN03')
#ggplot(AFN14, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN09')


#WTest1$ID <- 1
#WTest2$ID <- 1
#WTest3$ID <- 1
#WTest4$ID <- 1

#ggplot(WTest1, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('Bert1')
#ggplot(WTest2, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('Ber2')
#ggplot(WTest3, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('WPT389')
#ggplot(WTest4, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('WPT391')

#FNsub <- rbind(AFN12,AFN13,AFN14)
#FNsub$ID <- NA
#Tester <- rbind(AFN13,FNsub)
#Tester2 <- rbind(AFN12,FNsub)
#Tester3 <- rbind(AFN14,FNsub)
#ggplot(Tester2, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('M92-220')
#ggplot(Tester, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN03')
#ggplot(Tester3, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN09')




#FNsub2 <- rbind(WTest1,WTest2,WTest3,WTest4)
#FNsub2$ID <- NA

#T1 <- rbind(WTest1,FNsub2)
#T2 <- rbind(WTest2,FNsub2)
#T3 <- rbind(WTest3,FNsub2)
#T4 <- rbind(WTest4,FNsub2)

#ggplot(T1, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('Bert1')
#ggplot(T2, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('Ber2')
#ggplot(T3, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('WPT389')
#ggplot(T4, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('WPT391')


######################################################
###################################################
#FN code
###################################################
######################################################

#set working directory
#call ggplot for graphs
setwd("~/Downloads")
library(ggplot2)

#read in the .csv from GATK vcf that was converted into .csv form
FNAll <- read.csv("FNOnethruElevenHOMO.csv", header = TRUE)

#call hard SNPs if one columns is 1/1 and the rest have to be 0/0
FN01 <- FNAll[(FNAll$FN01 == '1/1')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
#make an empty vector
Line <- c()
#find out how many rows of info we have
for (i in 1:length(FN01[,1])){
  #make the line ID
  add <- "FN01"
  #add it to the end of the row
  Line <- append(Line,add)
}
#combine everything into a variable
FN01 <- cbind(FN01,Line)


#repeat for all of the other FN lines
#######
FN02 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '1/1') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN02[,1])){
  add <- "FN02"
  Line <- append(Line,add)
}
FN02 <- cbind(FN02,Line)
########
FN03 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '1/1')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN03[,1])){
  add <- "FN03"
  Line <- append(Line,add)
}
FN03 <- cbind(FN03,Line)
########
FN04 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '1/1')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN04[,1])){
  add <- "FN04"
  Line <- append(Line,add)
}
FN04 <- cbind(FN04,Line)
#########
FN05 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '1/1')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN05[,1])){
  add <- "FN05"
  Line <- append(Line,add)
}
FN05 <- cbind(FN05,Line)
#########
FN06 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '1/1')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN06[,1])){
  add <- "FN06"
  Line <- append(Line,add)
}
FN06 <- cbind(FN06,Line)
###########
FN07 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '1/1')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN07[,1])){
  add <- "FN07"
  Line <- append(Line,add)
}
FN07 <- cbind(FN07,Line)

FN08 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '1/1')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN08[,1])){
  add <- "FN08"
  Line <- append(Line,add)
}
FN08 <- cbind(FN08,Line)
#######
FN09 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '1/1')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN09[,1])){
  add <- "FN09"
  Line <- append(Line,add)
}
FN09 <- cbind(FN09,Line)
########
FN10 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '1/1')& (FNAll$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN10[,1])){
  add <- "FN10"
  Line <- append(Line,add)
}
FN10 <- cbind(FN10,Line)
########
FN11 <- FNAll[(FNAll$FN01 == '0/0')& (FNAll$FN02 == '0/0') &(FNAll$FN03 == '0/0')& (FNAll$FN04 == '0/0')& (FNAll$FN05 == '0/0')& (FNAll$FN06 == '0/0')& (FNAll$FN07 == '0/0')& (FNAll$FN08 == '0/0')& (FNAll$FN09 == '0/0')& (FNAll$FN10 == '0/0')& (FNAll$FN11 == '1/1'),]
Line <- c()
for (i in 1:length(FN11[,1])){
  add <- "FN11"
  Line <- append(Line,add)
}
FN11 <- cbind(FN11,Line)

#compile all of the data into FNlist
FNlist <- rbind(FN01,FN02,FN03,FN04,FN05,FN06,FN07,FN08,FN09,FN10,FN11)
#filter out the region of hetergenaity
FNlist <- FNlist[-c(which((FNlist$Chromosome == "Gm12") & (FNlist$Position > 10000000)& (FNlist$Position < 23000000) & (FNlist$Line == "FN07"))),]


#remove scaffolds
FNlist <- FNlist[- grep("scaffold", FNlist$Chromosome),]
#write to file
write.csv(FNlist, file = "FNVariantTableFinal.csv")

#pull out relecvant data for circos
CFN <- cbind(FNlist[,1:2],FNlist[,2],FNlist$Line)
row.names(CFN) <- NULL
colnames(CFN) <- c("chr","start","stop","id")
write.table(CFN, file = "FNCircos.txt", sep = " ", quote = FALSE, row.names = F)


##get circos SV information generation
FNSV <- read.csv("FN1through11SV.csv", header=F)
FNSV <- FNSV[,1:5]
row.names(FNSV) <- NULL
colnames(FNSV) <- c("chr","start","stop","SVtype","Line")
write.table(FNSV, file = "FNSVCircos.txt", sep = " ", quote = FALSE, row.names = F)

#plot graphs using ggplot of SNP distribution within FN lines
ggplot(FNlist, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + 
  facet_grid(Line~Chromosome) + theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  ylab("FN Line")+ggtitle('SNP distribution within FN Lines')
ggplot(FNlist, aes(Line, fill=Chromosome)) + geom_bar(position="dodge")+
  xlab('FN Line')+ylab('Number of SNPs')+ggtitle('SNP Distribution Box Plot')


#get a list of SNPs per chomosome for one line
a1 <- length(which(FN01$Chromosome == 'Gm01'))
#repeat for all the chromosomes and all of the other lines
a2 <- length(which(FN01$Chromosome == 'Gm02'))
a3 <- length(which(FN01$Chromosome == 'Gm03'))
a4 <- length(which(FN01$Chromosome == 'Gm04'))
a5 <- length(which(FN01$Chromosome == 'Gm05'))
a6 <- length(which(FN01$Chromosome == 'Gm06'))
a7 <- length(which(FN01$Chromosome == 'Gm07'))
a8 <- length(which(FN01$Chromosome == 'Gm08'))
a9 <- length(which(FN01$Chromosome == 'Gm09'))
a10 <- length(which(FN01$Chromosome == 'Gm10'))
a11 <- length(which(FN01$Chromosome == 'Gm11'))
a12 <- length(which(FN01$Chromosome == 'Gm12'))
a13 <- length(which(FN01$Chromosome == 'Gm13'))
a14 <- length(which(FN01$Chromosome == 'Gm14'))
a15 <- length(which(FN01$Chromosome == 'Gm15'))
a16 <- length(which(FN01$Chromosome == 'Gm16'))
a17 <- length(which(FN01$Chromosome == 'Gm17'))
a18 <- length(which(FN01$Chromosome == 'Gm18'))
a19 <- length(which(FN01$Chromosome == 'Gm19'))
a20 <- length(which(FN01$Chromosome == 'Gm20'))
#combine into a row of data
SumFN1 <- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN02$Chromosome == 'Gm01'))
a2 <- length(which(FN02$Chromosome == 'Gm02'))
a3 <- length(which(FN02$Chromosome == 'Gm03'))
a4 <- length(which(FN02$Chromosome == 'Gm04'))
a5 <- length(which(FN02$Chromosome == 'Gm05'))
a6 <- length(which(FN02$Chromosome == 'Gm06'))
a7 <- length(which(FN02$Chromosome == 'Gm07'))
a8 <- length(which(FN02$Chromosome == 'Gm08'))
a9 <- length(which(FN02$Chromosome == 'Gm09'))
a10 <- length(which(FN02$Chromosome == 'Gm10'))
a11 <- length(which(FN02$Chromosome == 'Gm11'))
a12 <- length(which(FN02$Chromosome == 'Gm12'))
a13 <- length(which(FN02$Chromosome == 'Gm13'))
a14 <- length(which(FN02$Chromosome == 'Gm14'))
a15 <- length(which(FN02$Chromosome == 'Gm15'))
a16 <- length(which(FN02$Chromosome == 'Gm16'))
a17 <- length(which(FN02$Chromosome == 'Gm17'))
a18 <- length(which(FN02$Chromosome == 'Gm18'))
a19 <- length(which(FN02$Chromosome == 'Gm19'))
a20 <- length(which(FN02$Chromosome == 'Gm20'))
SumFN2 <- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN04$Chromosome == 'Gm01'))
a2 <- length(which(FN04$Chromosome == 'Gm02'))
a3 <- length(which(FN04$Chromosome == 'Gm03'))
a4 <- length(which(FN04$Chromosome == 'Gm04'))
a5 <- length(which(FN04$Chromosome == 'Gm05'))
a6 <- length(which(FN04$Chromosome == 'Gm06'))
a7 <- length(which(FN04$Chromosome == 'Gm07'))
a8 <- length(which(FN04$Chromosome == 'Gm08'))
a9 <- length(which(FN04$Chromosome == 'Gm09'))
a10 <- length(which(FN04$Chromosome == 'Gm10'))
a11 <- length(which(FN04$Chromosome == 'Gm11'))
a12 <- length(which(FN04$Chromosome == 'Gm12'))
a13 <- length(which(FN04$Chromosome == 'Gm13'))
a14 <- length(which(FN04$Chromosome == 'Gm14'))
a15 <- length(which(FN04$Chromosome == 'Gm15'))
a16 <- length(which(FN04$Chromosome == 'Gm16'))
a17 <- length(which(FN04$Chromosome == 'Gm17'))
a18 <- length(which(FN04$Chromosome == 'Gm18'))
a19 <- length(which(FN04$Chromosome == 'Gm19'))
a20 <- length(which(FN04$Chromosome == 'Gm20'))
SumFN4<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN03$Chromosome == 'Gm01'))
a2 <- length(which(FN03$Chromosome == 'Gm02'))
a3 <- length(which(FN03$Chromosome == 'Gm03'))
a4 <- length(which(FN03$Chromosome == 'Gm04'))
a5 <- length(which(FN03$Chromosome == 'Gm05'))
a6 <- length(which(FN03$Chromosome == 'Gm06'))
a7 <- length(which(FN03$Chromosome == 'Gm07'))
a8 <- length(which(FN03$Chromosome == 'Gm08'))
a9 <- length(which(FN03$Chromosome == 'Gm09'))
a10 <- length(which(FN03$Chromosome == 'Gm10'))
a11 <- length(which(FN03$Chromosome == 'Gm11'))
a12 <- length(which(FN03$Chromosome == 'Gm12'))
a13 <- length(which(FN03$Chromosome == 'Gm13'))
a14 <- length(which(FN03$Chromosome == 'Gm14'))
a15 <- length(which(FN03$Chromosome == 'Gm15'))
a16 <- length(which(FN03$Chromosome == 'Gm16'))
a17 <- length(which(FN03$Chromosome == 'Gm17'))
a18 <- length(which(FN03$Chromosome == 'Gm18'))
a19 <- length(which(FN03$Chromosome == 'Gm19'))
a20 <- length(which(FN03$Chromosome == 'Gm20'))
SumFN3<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN05$Chromosome == 'Gm01'))
a2 <- length(which(FN05$Chromosome == 'Gm02'))
a3 <- length(which(FN05$Chromosome == 'Gm03'))
a4 <- length(which(FN05$Chromosome == 'Gm04'))
a5 <- length(which(FN05$Chromosome == 'Gm05'))
a6 <- length(which(FN05$Chromosome == 'Gm06'))
a7 <- length(which(FN05$Chromosome == 'Gm07'))
a8 <- length(which(FN05$Chromosome == 'Gm08'))
a9 <- length(which(FN05$Chromosome == 'Gm09'))
a10 <- length(which(FN05$Chromosome == 'Gm10'))
a11 <- length(which(FN05$Chromosome == 'Gm11'))
a12 <- length(which(FN05$Chromosome == 'Gm12'))
a13 <- length(which(FN05$Chromosome == 'Gm13'))
a14 <- length(which(FN05$Chromosome == 'Gm14'))
a15 <- length(which(FN05$Chromosome == 'Gm15'))
a16 <- length(which(FN05$Chromosome == 'Gm16'))
a17 <- length(which(FN05$Chromosome == 'Gm17'))
a18 <- length(which(FN05$Chromosome == 'Gm18'))
a19 <- length(which(FN05$Chromosome == 'Gm19'))
a20 <- length(which(FN05$Chromosome == 'Gm20'))
SumFN5<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN06$Chromosome == 'Gm01'))
a2 <- length(which(FN06$Chromosome == 'Gm02'))
a3 <- length(which(FN06$Chromosome == 'Gm03'))
a4 <- length(which(FN06$Chromosome == 'Gm04'))
a5 <- length(which(FN06$Chromosome == 'Gm05'))
a6 <- length(which(FN06$Chromosome == 'Gm06'))
a7 <- length(which(FN06$Chromosome == 'Gm07'))
a8 <- length(which(FN06$Chromosome == 'Gm08'))
a9 <- length(which(FN06$Chromosome == 'Gm09'))
a10 <- length(which(FN06$Chromosome == 'Gm10'))
a11 <- length(which(FN06$Chromosome == 'Gm11'))
a12 <- length(which(FN06$Chromosome == 'Gm12'))
a13 <- length(which(FN06$Chromosome == 'Gm13'))
a14 <- length(which(FN06$Chromosome == 'Gm14'))
a15 <- length(which(FN06$Chromosome == 'Gm15'))
a16 <- length(which(FN06$Chromosome == 'Gm16'))
a17 <- length(which(FN06$Chromosome == 'Gm17'))
a18 <- length(which(FN06$Chromosome == 'Gm18'))
a19 <- length(which(FN06$Chromosome == 'Gm19'))
a20 <- length(which(FN06$Chromosome == 'Gm20'))
SumFN6<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN07$Chromosome == 'Gm01'))
a2 <- length(which(FN07$Chromosome == 'Gm02'))
a3 <- length(which(FN07$Chromosome == 'Gm03'))
a4 <- length(which(FN07$Chromosome == 'Gm04'))
a5 <- length(which(FN07$Chromosome == 'Gm05'))
a6 <- length(which(FN07$Chromosome == 'Gm06'))
a7 <- length(which(FN07$Chromosome == 'Gm07'))
a8 <- length(which(FN07$Chromosome == 'Gm08'))
a9 <- length(which(FN07$Chromosome == 'Gm09'))
a10 <- length(which(FN07$Chromosome == 'Gm10'))
a11 <- length(which(FN07$Chromosome == 'Gm11'))
a12 <- length(which(Cake$Chromosome == 'Gm12'& Cake$Line == 'FN07'))
a13 <- length(which(FN07$Chromosome == 'Gm13'))
a14 <- length(which(FN07$Chromosome == 'Gm14'))
a15 <- length(which(FN07$Chromosome == 'Gm15'))
a16 <- length(which(FN07$Chromosome == 'Gm16'))
a17 <- length(which(FN07$Chromosome == 'Gm17'))
a18 <- length(which(FN07$Chromosome == 'Gm18'))
a19 <- length(which(FN07$Chromosome == 'Gm19'))
a20 <- length(which(FN07$Chromosome == 'Gm20'))
SumFN7<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN08$Chromosome == 'Gm01'))
a2 <- length(which(FN08$Chromosome == 'Gm02'))
a3 <- length(which(FN08$Chromosome == 'Gm03'))
a4 <- length(which(FN08$Chromosome == 'Gm04'))
a5 <- length(which(FN08$Chromosome == 'Gm05'))
a6 <- length(which(FN08$Chromosome == 'Gm06'))
a7 <- length(which(FN08$Chromosome == 'Gm07'))
a8 <- length(which(FN08$Chromosome == 'Gm08'))
a9 <- length(which(FN08$Chromosome == 'Gm09'))
a10 <- length(which(FN08$Chromosome == 'Gm10'))
a11 <- length(which(FN08$Chromosome == 'Gm11'))
a12 <- length(which(FN08$Chromosome == 'Gm12'))
a13 <- length(which(FN08$Chromosome == 'Gm13'))
a14 <- length(which(FN08$Chromosome == 'Gm14'))
a15 <- length(which(FN08$Chromosome == 'Gm15'))
a16 <- length(which(FN08$Chromosome == 'Gm16'))
a17 <- length(which(FN08$Chromosome == 'Gm17'))
a18 <- length(which(FN08$Chromosome == 'Gm18'))
a19 <- length(which(FN08$Chromosome == 'Gm19'))
a20 <- length(which(FN08$Chromosome == 'Gm20'))
SumFN8<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN09$Chromosome == 'Gm01'))
a2 <- length(which(FN09$Chromosome == 'Gm02'))
a3 <- length(which(FN09$Chromosome == 'Gm03'))
a4 <- length(which(FN09$Chromosome == 'Gm04'))
a5 <- length(which(FN09$Chromosome == 'Gm05'))
a6 <- length(which(FN09$Chromosome == 'Gm06'))
a7 <- length(which(FN09$Chromosome == 'Gm07'))
a8 <- length(which(FN09$Chromosome == 'Gm08'))
a9 <- length(which(FN09$Chromosome == 'Gm09'))
a10 <- length(which(FN09$Chromosome == 'Gm10'))
a11 <- length(which(FN09$Chromosome == 'Gm11'))
a12 <- length(which(FN09$Chromosome == 'Gm12'))
a13 <- length(which(FN09$Chromosome == 'Gm13'))
a14 <- length(which(FN09$Chromosome == 'Gm14'))
a15 <- length(which(FN09$Chromosome == 'Gm15'))
a16 <- length(which(FN09$Chromosome == 'Gm16'))
a17 <- length(which(FN09$Chromosome == 'Gm17'))
a18 <- length(which(FN09$Chromosome == 'Gm18'))
a19 <- length(which(FN09$Chromosome == 'Gm19'))
a20 <- length(which(FN09$Chromosome == 'Gm20'))
SumFN9<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN10$Chromosome == 'Gm01'))
a2 <- length(which(FN10$Chromosome == 'Gm02'))
a3 <- length(which(FN10$Chromosome == 'Gm03'))
a4 <- length(which(FN10$Chromosome == 'Gm04'))
a5 <- length(which(FN10$Chromosome == 'Gm05'))
a6 <- length(which(FN10$Chromosome == 'Gm06'))
a7 <- length(which(FN10$Chromosome == 'Gm07'))
a8 <- length(which(FN10$Chromosome == 'Gm08'))
a9 <- length(which(FN10$Chromosome == 'Gm09'))
a10 <- length(which(FN10$Chromosome == 'Gm10'))
a11 <- length(which(FN10$Chromosome == 'Gm11'))
a12 <- length(which(FN10$Chromosome == 'Gm12'))
a13 <- length(which(FN10$Chromosome == 'Gm13'))
a14 <- length(which(FN10$Chromosome == 'Gm14'))
a15 <- length(which(FN10$Chromosome == 'Gm15'))
a16 <- length(which(FN10$Chromosome == 'Gm16'))
a17 <- length(which(FN10$Chromosome == 'Gm17'))
a18 <- length(which(FN10$Chromosome == 'Gm18'))
a19 <- length(which(FN10$Chromosome == 'Gm19'))
a20 <- length(which(FN10$Chromosome == 'Gm20'))
SumFN10<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

a1 <- length(which(FN11$Chromosome == 'Gm01'))
a2 <- length(which(FN11$Chromosome == 'Gm02'))
a3 <- length(which(FN11$Chromosome == 'Gm03'))
a4 <- length(which(FN11$Chromosome == 'Gm04'))
a5 <- length(which(FN11$Chromosome == 'Gm05'))
a6 <- length(which(FN11$Chromosome == 'Gm06'))
a7 <- length(which(FN11$Chromosome == 'Gm07'))
a8 <- length(which(FN11$Chromosome == 'Gm08'))
a9 <- length(which(FN11$Chromosome == 'Gm09'))
a10 <- length(which(FN11$Chromosome == 'Gm10'))
a11 <- length(which(FN11$Chromosome == 'Gm11'))
a12 <- length(which(FN11$Chromosome == 'Gm12'))
a13 <- length(which(FN11$Chromosome == 'Gm13'))
a14 <- length(which(FN11$Chromosome == 'Gm14'))
a15 <- length(which(FN11$Chromosome == 'Gm15'))
a16 <- length(which(FN11$Chromosome == 'Gm16'))
a17 <- length(which(FN11$Chromosome == 'Gm17'))
a18 <- length(which(FN11$Chromosome == 'Gm18'))
a19 <- length(which(FN11$Chromosome == 'Gm19'))
a20 <- length(which(FN11$Chromosome == 'Gm20'))
SumFN11<- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)

FNTable <- rbind(SumFN1,SumFN2,SumFN3,SumFN4,SumFN5,SumFN6,SumFN7,SumFN8,SumFN9,SumFN10,SumFN11)
#relabel the row names and columns and write to a file
colnames(FNTable) <- c("Gm01","Gm02","Gm03","Gm04","Gm05","Gm06","Gm07","Gm08","Gm09","Gm10","Gm11","Gm12","Gm13","Gm14","Gm15","Gm16","Gm17","Gm18","Gm19","Gm20")
row.names(FNTable) <- c("FN01","FN02","FN03","FN04","FN05","FN06","FN07","FN08","FN09","FN10","FN11")
write.csv(FNTable, file = "PerLinePerChromosomeFNTableFinal.csv")

##this is just dealing with the region of heterogenaity for comparison
temp <- which(FN07$Chromosome == 'Gm12')
temp2 <- which(FN07$Chromosome == 'Gm13')
FN7GM12 <- FN07[temp,]
FN7GM13 <- FN07[temp2,]
p1 <- ggplot(FN7GM12, aes(Position, ID)) + geom_point()+ggtitle("FN07 GM12 SNPs")
p2 <- ggplot(FN7GM13, aes(Position, ID)) + geom_point()+ ggtitle("FN07 GM13 SNPs")


#call on the multiplot function from the R cookbook
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



multiplot(p1,p2)

#use this as a what to plot the SNP's across al FN lines
#make sure FN list have no values
T1 <- rbind(FN01,FNlist)
T2 <- rbind(FN02,FNlist)
T3 <- rbind(FN03,FNlist)
T4 <- rbind(FN04,FNlist)
T5 <- rbind(FN05,FNlist)
T6 <- rbind(FN06,FNlist)
T7 <- rbind(FN07,FNlist)
T8 <- rbind(FN08,FNlist)
T9 <- rbind(FN09,FNlist)
T10 <- rbind(FN10,FNlist)
T11 <- rbind(FN11,FNlist)

a <- ggplot(T1, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN01')
b <- ggplot(T2, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN02')
c <- ggplot(T3, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN03')
d <- ggplot(T4, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN04')
e <- ggplot(T5, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN05')
f <- ggplot(T6, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN06')
g <- ggplot(T7, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN07')
h <- ggplot(T8, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN08')
i <- ggplot(T9, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN09')
j <- ggplot(T10, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN10')
k <- ggplot(T11, aes(Position, ID,colour = factor(Chromosome))) + geom_point() + facet_grid(.~Chromosome)+ylab('FN11')
multiplot(a,b,c,d,e,f,g,h,i,j,k, cols=1)
#this is for het data that we didnt work on
#FNAll2 <- read.csv("FNOnethruElevenHET.csv", header = TRUE)
#FNAll2 <- as.data.frame(FNAll2)


#FNAll2$FN03 <- as.character(FNAll2$FN03)
#FNAll2$FN06 <- as.character(FNAll2$FN06)
#FNTest <- FNAll2[(FNAll2$FN06 == '0/0')& (FNAll2$FN03 == '1/1'),]



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#####Transitions and Transversions
setwd("~/Downloads")
#read in data
WPT <- read.csv("WPTVariantTable-2.csv")
FN <- read.csv("FNVariantTableFinal.csv")
#conpare how many bases change to what for all possibilities for all lines
A2C <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_1")))
A2G <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_1")))
A2T <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_1")))

C2A <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_1")))
C2G <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_1")))
C2T <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_1")))

G2A <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_1")))
G2C <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_1")))
G2T <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_1")))

T2A <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_1")))
T2C <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_1")))
T2G <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_1")))

Bert1 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_2")))
A2G <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_2")))
A2T <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_2")))

C2A <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_2")))
C2G <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_2")))
C2T <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_2")))

G2A <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_2")))
G2C <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_2")))
G2T <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "T")& (WPT$Line == "BertMN01_2")))

T2A <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "A")& (WPT$Line == "BertMN01_2")))
T2C <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "C")& (WPT$Line == "BertMN01_2")))
T2G <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "G")& (WPT$Line == "BertMN01_2")))

Bert2 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

A2C <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT389_2_2")))
A2G <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT389_2_2")))
A2T <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT389_2_2")))

C2A <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT389_2_2")))
C2G <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT389_2_2")))
C2T <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT389_2_2")))

G2A <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT389_2_2")))
G2C <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT389_2_2")))
G2T <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT389_2_2")))

T2A <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT389_2_2")))
T2C <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT389_2_2")))
T2G <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT389_2_2")))

WPT389 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

A2C <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT391_1_6")))
A2G <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT391_1_6")))
A2T <- length(which((WPT$Ref_Base == "A") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT391_1_6")))

C2A <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT391_1_6")))
C2G <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT391_1_6")))
C2T <- length(which((WPT$Ref_Base == "C") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT391_1_6")))

G2A <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT391_1_6")))
G2C <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT391_1_6")))
G2T <- length(which((WPT$Ref_Base == "G") & (WPT$Alt_Base == "T")& (WPT$Line == "WPT391_1_6")))

T2A <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "A")& (WPT$Line == "WPT391_1_6")))
T2C <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "C")& (WPT$Line == "WPT391_1_6")))
T2G <- length(which((WPT$Ref_Base == "T") & (WPT$Alt_Base == "G")& (WPT$Line == "WPT391_1_6")))

WPT391 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN01")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN01")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN01")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN01")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN01")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN01")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN01")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN01")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN01")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN01")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN01")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN01")))

FN1 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN02")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN02")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN02")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN02")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN02")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN02")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN02")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN02")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN02")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN02")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN02")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN02")))

FN2 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN03")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN03")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN03")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN03")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN03")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN03")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN03")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN03")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN03")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN03")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN03")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN03")))

FN3 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN04")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN04")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN04")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN04")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN04")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN04")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN04")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN04")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN04")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN04")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN04")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN04")))

FN4 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN05")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN05")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN05")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN05")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN05")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN05")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN05")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN05")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN05")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN05")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN05")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN05")))

FN5 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN06")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN06")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN06")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN06")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN06")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN06")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN06")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN06")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN06")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN06")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN06")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN06")))

FN6 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN07")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN07")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN07")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN07")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN07")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN07")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN07")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN07")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN07")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN07")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN07")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN07")))

FN7 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN08")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN08")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN08")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN08")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN08")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN08")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN08")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN08")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN08")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN08")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN08")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN08")))

FN8 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN09")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN09")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN09")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN09")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN09")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN09")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN09")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN09")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN09")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN09")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN09")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN09")))

FN9 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN10")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN10")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN10")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN10")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN10")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN10")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN10")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN10")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN10")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN10")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN10")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN10")))

FN10 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)


A2C <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "C")& (FN$Line == "FN11")))
A2G <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "G")& (FN$Line == "FN11")))
A2T <- length(which((FN$Ref_Base == "A") & (FN$Alt_Base == "T")& (FN$Line == "FN11")))

C2A <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "A")& (FN$Line == "FN11")))
C2G <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "G")& (FN$Line == "FN11")))
C2T <- length(which((FN$Ref_Base == "C") & (FN$Alt_Base == "T")& (FN$Line == "FN11")))

G2A <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "A")& (FN$Line == "FN11")))
G2C <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "C")& (FN$Line == "FN11")))
G2T <- length(which((FN$Ref_Base == "G") & (FN$Alt_Base == "T")& (FN$Line == "FN11")))

T2A <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "A")& (FN$Line == "FN11")))
T2C <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "C")& (FN$Line == "FN11")))
T2G <- length(which((FN$Ref_Base == "T") & (FN$Alt_Base == "G")& (FN$Line == "FN11")))

FN11 <- cbind(A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G)

#combine all of the data and write it to a file
TransitionTransversion <- rbind(Bert1,Bert2,WPT389,WPT391,FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,FN11)
row.names(TransitionTransversion) <- c("BertMN01-1","Bert1MN01-2","WPT389","WPT391","FN1","FN2","FN3","FN4","FN5","FN6","FN7","FN8","FN9","FN10","FN11")
write.csv(TransitionTransversion, "TransitionTransversionbyLineFinal.csv")


##############ARCHER, MINSOY, NOIR
#compare SNP calling to other lines
setwd("~/Downloads")
Archer <- read.csv("Archer.csv")
T2A <- length(which(Archer$Archer == "1/1"))

Noir1 <- read.csv("Noir1.csv")
T2A2 <- length(which(Noir1$Archer == "1/1"))


####### Exon and Gene analysis
setwd("~/Downloads")
#read in the data and screen for the words gene or exon
WPTExon <- read.csv("HomoWPTComboExon.csv")
WPTExon <- WPTExon[which(WPTExon$ID == "exon"),]
WPTGene <- read.csv("HomoWPTComboGene.csv")
WPTGene <- WPTGene[which(WPTGene$ID == "gene"),]
WPTExon1 <- WPTExon[WPTExon$Bert1 == '1/1' & WPTExon$Bert2 == '0/0' & WPTExon$WPT389 == '0/0' & WPTExon$WPT391 == '0/0',]
WPTExon2 <- WPTExon[WPTExon$Bert1 == '0/0' & WPTExon$Bert2 == '1/1' & WPTExon$WPT389 == '0/0' & WPTExon$WPT391 == '0/0',]
WPTExon3 <- WPTExon[WPTExon$Bert1 == '0/0' & WPTExon$Bert2 == '0/0' & WPTExon$WPT389 == '1/1' & WPTExon$WPT391 == '0/0',]
WPTExon4 <- WPTExon[WPTExon$Bert1 == '0/0' & WPTExon$Bert2 == '0/0' & WPTExon$WPT389 == '0/0' & WPTExon$WPT391 == '1/1',]

WPTGene1 <- WPTGene[WPTGene$Bert1 == '1/1' & WPTGene$Bert2 == '0/0' & WPTGene$WPT389 == '0/0' & WPTGene$WPT391 == '0/0',]
WPTGene2 <- WPTGene[WPTGene$Bert1 == '0/0' & WPTGene$Bert2 == '1/1' & WPTGene$WPT389 == '0/0' & WPTGene$WPT391 == '0/0',]
WPTGene3 <- WPTGene[WPTGene$Bert1 == '0/0' & WPTGene$Bert2 == '0/0' & WPTGene$WPT389 == '1/1' & WPTGene$WPT391 == '0/0',]
WPTGene4 <- WPTGene[WPTGene$Bert1 == '0/0' & WPTGene$Bert2 == '0/0' & WPTGene$WPT389 == '0/0' & WPTGene$WPT391 == '1/1',]

write.csv(WPTGene3, "WTP389_2_2_in_Gene.csv")

FNExon <- read.csv("HomoFNComboExon.csv")
FNGene <- read.csv("HomoFNComboGene.csv")
FNExon <- FNExon[c(which(FNExon$ID == 'exon')),]
FNGene <- FNGene[c(which(FNGene$ID == 'gene')),]

#same as before add the line name to the end of the row
FN1Exon <- FNExon[(FNExon$FN1 == '1/1')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN1Exon[,1])){
  add <- "FN01"
  Line <- append(Line,add)
}
FN01 <- cbind(FN1Exon,Line)

FN2Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '1/1') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN2Exon[,1])){
  add <- "FN02"
  Line <- append(Line,add)
}
FN02 <- cbind(FN2Exon,Line)

FN3Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '1/1')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN3Exon[,1])){
  add <- "FN03"
  Line <- append(Line,add)
}
FN03 <- cbind(FN3Exon,Line)


FN4Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '1/1')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN4Exon[,1])){
  add <- "FN04"
  Line <- append(Line,add)
}
FN04 <- cbind(FN4Exon,Line)

FN5Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '1/1')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN5Exon[,1])){
  add <- "FN05"
  Line <- append(Line,add)
}
FN05 <- cbind(FN5Exon,Line)


FN6Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '1/1')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN6Exon[,1])){
  add <- "FN06"
  Line <- append(Line,add)
}
FN06 <- cbind(FN6Exon,Line)


FN7Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '1/1')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN7Exon[,1])){
  add <- "FN07"
  Line <- append(Line,add)
}
FN07 <- cbind(FN7Exon,Line)

FN8Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '1/1')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN8Exon[,1])){
  add <- "FN08"
  Line <- append(Line,add)
}
FN08 <- cbind(FN8Exon,Line)


FN9Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '1/1')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN9Exon[,1])){
  add <- "FN09"
  Line <- append(Line,add)
}
FN09 <- cbind(FN9Exon,Line)

FN10Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '1/1')& (FNExon$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN10Exon[,1])){
  add <- "FN10"
  Line <- append(Line,add)
}
FN10 <- cbind(FN10Exon,Line)


FN11Exon <- FNExon[(FNExon$FN1 == '0/0')& (FNExon$FN2 == '0/0') &(FNExon$FN3 == '0/0')& (FNExon$FN4 == '0/0')& (FNExon$FN5 == '0/0')& (FNExon$FN6 == '0/0')& (FNExon$FN7 == '0/0')& (FNExon$FN8 == '0/0')& (FNExon$FN9 == '0/0')& (FNExon$FN10 == '0/0')& (FNExon$FN11 == '1/1'),]
Line <- c()
for (i in 1:length(FN11Exon[,1])){
  add <- "FN11"
  Line <- append(Line,add)
}
FN11 <- cbind(FN11Exon,Line)
#combine the data and write to a file
FNV3 <- rbind(FN01,FN02,FN03,FN04,FN05,FN06,FN07,FN08,FN09,FN10,FN11)

write.csv(FNV3, "FN11ExoThrough11HomoExonSNPs.csv")

##Run the same thing for genes
FN1Gene <- FNGene[(FNGene$FN1 == '1/1')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN2Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '1/1') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN3Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '1/1')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN4Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '1/1')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN5Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '1/1')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN6Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '1/1')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN7Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '1/1')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN8Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '1/1')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN9Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '1/1')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
FN10Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '1/1')& (FNGene$FN11 == '0/0'),]
FN11Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '1/1'),]


FN1Gene <- FNGene[(FNGene$FN1 == '1/1')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN1Gene[,1])){
  add <- "FN01"
  Line <- append(Line,add)
}
FN01 <- cbind(FN1Gene,Line)

FN2Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '1/1') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN2Gene[,1])){
  add <- "FN02"
  Line <- append(Line,add)
}
FN02 <- cbind(FN2Gene,Line)

FN3Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '1/1')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN3Gene[,1])){
  add <- "FN03"
  Line <- append(Line,add)
}
FN03 <- cbind(FN3Gene,Line)


FN4Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '1/1')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN4Gene[,1])){
  add <- "FN04"
  Line <- append(Line,add)
}
FN04 <- cbind(FN4Gene,Line)

FN5Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '1/1')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN5Gene[,1])){
  add <- "FN05"
  Line <- append(Line,add)
}
FN05 <- cbind(FN5Gene,Line)


FN6Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '1/1')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN6Gene[,1])){
  add <- "FN06"
  Line <- append(Line,add)
}
FN06 <- cbind(FN6Gene,Line)


FN7Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '1/1')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN7Gene[,1])){
  add <- "FN07"
  Line <- append(Line,add)
}
FN07 <- cbind(FN7Gene,Line)

FN8Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '1/1')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN8Gene[,1])){
  add <- "FN08"
  Line <- append(Line,add)
}
FN08 <- cbind(FN8Gene,Line)


FN9Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '1/1')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN9Gene[,1])){
  add <- "FN09"
  Line <- append(Line,add)
}
FN09 <- cbind(FN9Gene,Line)

FN10Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '1/1')& (FNGene$FN11 == '0/0'),]
Line <- c()
for (i in 1:length(FN10Gene[,1])){
  add <- "FN10"
  Line <- append(Line,add)
}
FN10 <- cbind(FN10Gene,Line)


FN11Gene <- FNGene[(FNGene$FN1 == '0/0')& (FNGene$FN2 == '0/0') &(FNGene$FN3 == '0/0')& (FNGene$FN4 == '0/0')& (FNGene$FN5 == '0/0')& (FNGene$FN6 == '0/0')& (FNGene$FN7 == '0/0')& (FNGene$FN8 == '0/0')& (FNGene$FN9 == '0/0')& (FNGene$FN10 == '0/0')& (FNGene$FN11 == '1/1'),]
Line <- c()
for (i in 1:length(FN11Gene[,1])){
  add <- "FN11"
  Line <- append(Line,add)
}
FN11 <- cbind(FN11Gene,Line)

FNV4 <- rbind(FN01,FN02,FN03,FN04,FN05,FN06,FN07,FN08,FN09,FN10,FN11)

write.csv(FNV4, "FN1Through11HomoGeneSNPs.csv")








