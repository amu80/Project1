#load
library(biomaRt)          #Contains annotations and seqeunces 
library(seqinr)           #Write fasta file 
library(ComplexHeatmap)   #for complex heatmap
library(circlize)         #for adapting complex heatmap
library(RColorBrewer)     #for complex heatmap colouring 
library(readxl)           #for reading excel files
library(scanMiR)          #predicting pseudogene/miRNA targets 
library(Biostrings)       #needed for using TargetSeed function in ScanMiR
library(dplyr)            #needed for counting across row matches in network
library(igraph)           #for network visualisation 
library(circlize)         #for circular network visualisation
library(reshape2)         #for reshaping data for ggplot box plot 
library(ggplot2)          #for multiple box plot 
library(calecopal)        #for correlation expression plot 

#set the directory 
setwd("/Users/amugrg/Documents/P1/proj1")

#install.packages("Biostrings")
#########
################
#################################
#Obtaining names of protein coding genes and Pseudogenes found in patient expression count data 

#connect to human ensembl dataset
myMarts <- listMarts()     
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
myDatasets <- listDatasets(ensembl)  #can see the file where human gene annotation is stored
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")  #will list you all genes
                                                                          #pseudogenes is under attributes ->  biotypes

################################Patient ProteinCodingGenes:
genes <- getBM(mart= ensembl,
               attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","chromosome_name", "start_position",
                              "end_position"),
               filters = c("biotype", "chromosome_name"),
               values = list("protein_coding", c(1:22, "X", "Y")))           #filter only protein coding annotations for chromosomes only
colnames(genes) <- c("ID", "Name", "Type", "Chromosome", "Start",  "End")     

#patient data
counts <- read.table("counts.tab")    #Patient ALL gene expression counts 
columns <- read.table("columns.tab")  #patient extra information

#subset the patient count data to only contain protein coding genes 
genecounts <- counts[rownames(counts) %in% genes$Name, ] 

#generate heat map expression 
#heatmap(as.matrix(genecounts))

################################Patient pseudogenes:
#Obtain all gene annotations for all chromosomes 
genes1 <- getBM(mart= ensembl,
                attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","chromosome_name",
                               "start_position","end_position"),
                filters = c("chromosome_name"),
                values = list(c(1:22, "X", "Y")))
colnames(genes1) <- c("ID", "Name", "Type", "Chromosome", "Start",  "End")           #63,150 rows 

#filter human gene annotation dataset to only contain pseudogenes
pseudogenes <- genes1[grepl("pseudogene", genes1$Type, ignore.case = T), ]          #15,216 rows 

#subset to only transcribed Processed pseudogenes
transcribedPseudogenes <- pseudogenes[grepl("transcribed", pseudogenes$Type, ignore.case = T), ]  #1633 rows

#subset patient specific transcribed processed pseudogenes
transcribedPseudogenesCounts <- counts[rownames(counts) %in% transcribedPseudogenes$Name, ] #669 rows 

##################################Visualise: Heatmap:

#####Simple Heatmap for pseudogenes 
heatmap(as.matrix(transcribedPseudogenesCounts))

#####Complex HeatMap for pseudogenes
#log heatmap as huge variation in count values
mat <- transcribedPseudogenesCounts
mat <- mat + 1                      #cant log2 transform 0 values
mat <- as.matrix(mat)
mat <- log2(mat)

col_fun <- brewer.pal(9,"YlOrBr")   #pick colourway 
heatmap <- Heatmap(mat, name = "mat",
                   column_title = "Pseudogenes Expression Heatmap",
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                   row_title = "Genes",
                   row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                   row_names_gp = gpar(fontsize = 5),
                   col = col_fun,
#                   column_order = order(as.numeric(gsub("Tl6_"," ", colnames(mat)))), <- loses the clustering
                   column_labels = gsub("Tl6_","", colnames(mat)))
draw(heatmap)


##########
###################
###############################Obtaining sequences for both patient protein coding genes and pseudogenes

#################protein coding cDNA seqeunces 
#first obtain list for all protein coding gene ID and transcript ID 
transcripts <- getBM(mart = ensembl,
                     attributes = c("ensembl_gene_id", "ensembl_transcript_id","external_gene_name"),
                     filters = c("biotype", "chromosome_name"),
                     values = list("protein_coding", c(1:22, "X", "Y")))
colnames(transcripts) <- c("GeneId", "TranscriptId","Name")                        #170,463

#obtain cDNA sequences with their transcriptID  for all protein coding gens 
cDNAs <- biomaRt::getSequence(mart = ensembl, type = "ensembl_transcript_id", id = transcripts$TranscriptId, seqType = "cdna")
colnames(cDNAs) <- c("Sequence", "TranscriptId")                                   #170,463 rows 

#Add Gene ID
genecDNAs <- merge(transcripts, cDNAs, by = "TranscriptId", all.x = T)             #170,463 rows 

#save for patient specific protein coding cDNA seqeunces 
genecDNAs <- genecDNAs[genecDNAs$Name %in% rownames(counts), ]  #subset to patients only 
genecDNAs <- genecDNAs[order(genecDNAs$GeneId, genecDNAs$TranscriptId), c(2, 1, 3, 4)]
rownames(genecDNAs) <- NULL
genecDNAs <- genecDNAs[genecDNAs$Sequence != "Sequence unavailable", ]        #168,166 rows compared to inital 19,194 names

#save
write.fasta(sequences = as.list(genecDNAs$Sequence), names = sprintf("%s|%s|%s", genecDNAs$GeneId, genecDNAs$TranscriptId, genecDNAs$Name), file.out = "Patient_Protein_cDNAs.fa")

############Pseudogenes cDNA sequences
transcripts1 <- getBM(mart= ensembl,
                attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","gene_biotype"),
                filters = c("chromosome_name"),
                values = list(c(1:22, "X", "Y")))
colnames(transcripts1) <- c("GeneId","TranscriptId","Name","Type")           #252,893 rows

#subset to only pseudogenes then transcribed only 
transcripts1 <- transcripts1[grepl("pseudogene", transcripts1$Type, ignore.case = T), ]  #15,228
transcripts1 <- transcripts1[grepl("transcribed", transcripts1$Type, ignore.case = T), ]  #1,626 rows 
transcripts1$Type <- NULL

#cDNAs for processed_transcribed_pseudogenes
cDNAs1 <- biomaRt::getSequence(mart = ensembl, type = "ensembl_transcript_id", id = transcripts1$TranscriptId, seqType = "cdna") 
colnames(cDNAs1) <- c("Sequence", "TranscriptId")

#add gene IDS
genecDNAs1 <- merge(transcripts1, cDNAs1, by = "TranscriptId", all.x = T)

#save for patient specific processed_transcribed_pseudogenes cDNA seqeunces 
genecDNAs1 <- genecDNAs1[genecDNAs1$Name %in% rownames(transcribedPseudogenesCounts), ] 
genecDNAs1 <- genecDNAs1[order(genecDNAs1$GeneId, genecDNAs1$TranscriptId), c(2, 1, 3, 4)]
rownames(genecDNAs1) <- NULL
genecDNAs1 <- genecDNAs1[genecDNAs1$Sequence != "Sequence unavailable", ]       #671 rows

#save all list
write.fasta(sequences = as.list(genecDNAs1$Sequence), names = sprintf("%s|%s|%s", genecDNAs1$GeneId, genecDNAs1$TranscriptId, genecDNAs1$Name), file.out = "Patient_Pseudo_Seq.fa")


############
####################
###################################Finding ProteinCoding/miRNA targets in patients
#miRTarBase

#read in the file
mTarAll <- read_excel("hsa_MTI.xlsx")     #contains all homo sapiens version 9.0 #502652 rows
mTarPatients <- mTarAll[mTarAll$'Target Gene' %in% rownames(genecounts), ] #subset it the genes in patients #480929 rows
proteinmiRNA <- mTarPatients[mTarPatients$'Support Type' == 'Functional MTI',]     #10,575
proteinmiRNA <- proteinmiRNA[c('miRNA','Target Gene')]
colnames(proteinmiRNA) <- c('miRNA','TargetGene')
proteinmiRNA <- unique(proteinmiRNA) #7997 entries -> because same miRNA/protein found in multiple experiments 
#proteinmiRNA <- proteinmiRNA[!duplicated(proteinmiRNA),]


##############
#####################
#################################Finding pseudogenes/miRNA targets in patients
###need sequences for this 

########list of all miRNA seqeunces from miRBase
miRBaseAll <- read.fasta(file = "mature.fa")  #list with name and sequence 
miRBaseRef <- miRBaseAll[grep("^hsa", names(miRBaseAll), ignore.case = T)] #humans only miRNA reference sequences #2656

#save file 
write.fasta(sequences = as.list(miRBaseRef), names = names(miRBaseRef), file.out = "miRNAhsaref.fa") #2656 

#############Predicting targets using scanMiR

#read in PATIENT pseudogene seqeunces as DNAstringSet required by prediction findSeedMatches function 
SampleTranscript <- Biostrings::readDNAStringSet("Patient_pseudo_Seq.fa")  #667 rows

miRNA <- Biostrings::readRNAStringSet("miRNAhsaref.fa")   #read in miRNA from mirBase sequences as RNAstringset 
miRNA <- as.character(miRNA)
miRNA <- miRNA[which(nchar(miRNA) <= 26)]        #checkSeedsInput function can only take in length = 26 or less
miRNA <- unname(miRNA)
type(miRNA)

#prediction function 
matches <- findSeedMatches(SampleTranscript, miRNA, verbose = FALSE)
#outputs matches  GRanges object 
matches      #938097 

##################
###############################
##########################################Creating proteinCoding/miRNA/pseudogene Network 

#####Building pseudogene/miRNA network first 
#export GRange object to a dataframe 
df <- as.data.frame(matches)    
duplicated <- df[duplicated(df),]

#turning miRNA sequence from DNA back to RNA
a <- DNAStringSet(df$miRNA)     
a <- as.character(a)             #dataframe cant hold DNAstringSEt 
a <- gsub("T", "U", a)           #convert seq from DNA to RNA
df$miRNA <- a

##Check couple manually for matching seqs-> in word document -> okay

###Need to add miRNA names to pseudogenes/miRNA matches 
miRNA1 <- Biostrings::readRNAStringSet("miRNAhsaref.fa")         #contains ALL miRNA sequences and their names 
miRNA1 <- as.character(miRNA1)
miRNA2 <- unname(miRNA1)  #contains just sequences 
miRNA1 <- as.list(miRNA1) #contains names and sequences 

length(miRNA2) #check unnaming done correctly 
length(miRNA1)

df1 <- data.frame(Name = names(miRNA1), Sequence = miRNA2, stringsAsFactors = FALSE) #contains ALL miRNA names with their seqeunces 
colnames(df1) <- c("Name", "miRNA")

#Merge the dataframe of matches and dataframe of names/sequences by their sequences
pseudogenemiRNA <- merge(df, df1, by="miRNA")     

#Need to separate the pseudogene gene name 
x <- strsplit(as.character(pseudogenemiRNA$seqnames), "\\|")
x <- sapply(x, `[`, 3) #get the gene names 
pseudogenemiRNA$Gene <- x    #add the gene name to the matches dataframe 
pseudogenemiRNA #contains final pseudogene miRNA name pseudogene names genes seqeunces mathces list 
#951,561

write.csv(pseudogenemiRNA, "pseudogenematcheswithgenenames.csv", row.names = FALSE)
#read.csv("pseudogenematcheswithgenenames.csv")

#################################merge protein-coding/pseudogene/miRNA 
#rename both separate dataframe so merging columns names are the same

#protein coding/miRNA
df3 <- proteinmiRNA #7,997 entries
df3 <- df3[c('miRNA','TargetGene')]
sum(duplicated(df3)) #none
colnames(df3) <- c("miRNA", "ProteinCodingGene")

#pseudogene/miRNA
df4 <- pseudogenemiRNA #951,561
df4 <- df4[c('Name','Gene')]      
colnames(df4) <- c("miRNA","Pseudogene")
sum(duplicated(df4)) #37,274
duplicated2 <- df4[duplicated(df4),]
df4 <- unique(df4)     #556,362               #only interested in one unique miRNA/pseudogene match -> not whether one is 6-mer/7-mer whatever

#merge to get protein coding/ miRNA/ Pseudogene
network <- merge(df3, df4, by="miRNA")
network <- network[,c('ProteinCodingGene','miRNA','Pseudogene')] #1,646,631 entries
networkkkk <- unique(network) #1,646,631 entries

#FER1L4/PTEN pair


#protein coding and pseudogen competition
CodingPseudo <- network[,c('ProteinCodingGene','Pseudogene')]   #230, 537 entries 
#un <- unique(CodingPseudo) #124 870 entries

##############
#########################
########################################Summary stats + visualisation

#################################################summary stats 
#######top individual hits 
summary(network)
network$ProteinCodingGene <- factor(network$ProteinCodingGene)
network$miRNA <- factor(network$miRNA)
network$Pseudogene <- factor(network$Pseudogene)
d <- summary(network)
d <- summary(network, maxsum = 21)
d <- as.data.frame.matrix(d, row.names = FALSE) 


#######figure for paper
#gave up not enough time to research how to split a table
summary <- read_excel("summary.xlsx")
as.data.frame(summary)
summary[] <- lapply(summary, trimws)       #trim the white space in the whole dataframe


cgc <- read.delim('Cosmic_CancerGeneCensus_v99_GRCh38.tsv', sep="\t")   #743 entries 

# our genes dataset is filtered to chromosomes 1:22/X/Y but some of the CGC are on alternative scaffolds 
# so filter those out since they'll never match our other analyses
cgc <- cgc[cgc$GENE_SYMBOL %in% genes$Name, ]       #742 entries

# expand roles in cancer
cgc$Oncogene <- grepl("oncogene", cgc$ROLE_IN_CANCER)
cgc$TumourSuppressor <- grepl("TSG", cgc$ROLE_IN_CANCER)
cgc$Fusion <- grepl("fusion", cgc$ROLE_IN_CANCER)

summary <- merge(summary, cgc, by.x = 'ProteinCodingGene', by.y = 'GENE_SYMBOL', all.x = TRUE)   
summary <- summary[c('ProteinCodingGene','No...2','Oncogene','TumourSuppressor','Fusion','miRNA','No...4','Pseudogene','No...6')]
View(summary)
as.numeric(summary$No...2)
type(summary$No...2)

#read in miRNA data
miRCan <- read.delim('miRCancerJune2020.txt', sep="\t")
miRCan$BladderCancer <- grepl("bladder cancer", miRCan$Cancer)
miRCan <- miRCan[c('mirId','BladderCancer')]

#merge summary with miRCan
summary <- merge(summary, miRCan, by.x ="miRNA", by.y = "mirId", all.x = TRUE)


#miRCancerJune2020.txt

#save
install.packages("writexl")
library(writexl)
# Add a existing iris  data set
write_xlsx(summary, "summary1.xlsx")





#######top 3-ways interaction hits -> not possible since there are no duplicates in the network

#######top pseudogene/protein coding gene PAIRS that are targetted by the same miRNA (regardless of what miRNA it is)
e <- CodingPseudo
e <- e %>%                           #Add count number to duplicate pairs 
  add_count(across(everything()))
e <- unique(e)                       #remove duplicate rows but keep row count number of matching pairs 
colnames(e) <- c('ProteinCodingGene','Pseudogenes','Number')
e <- e[order(e$Number, decreasing = TRUE), ]
e <- e[1:50, ]

#########################
#####################################
###################################################network analysis 

#################################overall for everything where Protein/Pseudo/miRNA are all nodes 
link <- network          #1,646,631
link <- unique(link)    #1,646,631      #there are no duplicates 
net <- graph_from_data_frame(d=link, directed=F)     #undirected graph 
V(net)$name[degree(net)==max(degree(net))] #hsa-miR-145-5p



#deg <- degree(net, mode="all")                       #node deg
#V(net)$size <- deg*3
#V(net)$color <- ifelse(V(net)$name %in% link$Pseudogene,"gold","skyblue") #need to change the colours too 
#plot(net) 
#too big does not plot -> didn't know what to subset it by 



###########################overall for top 50 protein/pseudogene matches, miRNA number is edge thickness
#edge list -> protein/pseudo with weighting(miRNA numbers)
link1 <- CodingPseudo
link1 <- link1 %>%                           #Add count number to duplicate pairs 
  add_count(across(everything()))
link1 <- unique(link1)
colnames(link1) <- c('ProteinCodingGene','Pseudogene','Weight')
f <- link1                                  #for later function
link1 <- link1[order(link1$Weight, decreasing = TRUE), ]
link1 <- link1[1:50, ]
rownames(link1) <- NULL

#Convert raw data to igraph object and plot 
net1 <- graph_from_data_frame(d=link1, directed=F)     #undirected graph 
E(net1)$width <- E(net1)$Weight/4                      #Determine edge widgth by weight/6 for visualisation     
# node degrees 
deg1 <- degree(net1, mode="all")
V(net1)$size <- deg1*1.5
#setting colour
V(net1)$color <- ifelse(V(net1)$name %in% link1$Pseudogene,"skyblue","pink")
#plot(net1, vertex.label.color="blue")  
plot(net1, vertex.label.color = ifelse(V(net1)$name %in% link1$Pseudogene, "blue", "red"))


  
#save plot
png("TOP50network.png", width = 25, height = 25, unit = "cm", res = 600)
plot(net1) 
dev.off()

###########################Network For individual ProteinCoding or pseudogene

f #contains Pseudogene/protein targets with weightings for the whole list  

subsett <- function(i,y){
  if (i %in% f$ProteinCodingGene){
    g <- f[grep(i, f$ProteinCodingGene), ]
    g <- g[order(g$Weight, decreasing = TRUE), ]
    g <- g[1:y, ]
    net <- graph_from_data_frame(d=g, directed=F)     #undirected graph 
    E(net)$width <- E(net)$Weight/4
    V(net)$color <- ifelse(V(net)$name %in% g$ProteinCodingGene,"gold","skyblue")
  } else {
    g <- f[grep(i, f$Pseudogene), ]
    g <- g[order(sx$Weight, decreasing = TRUE), ]
    g <- g[1:y, ]
    net <- graph_from_data_frame(d=g, directed=F)
    E(net)$width <- E(net)$Weight/4
    V(net)$color <- ifelse(V(net)$name %in% g$Pseudogene, "skyblue","gold")
  }
  plot(net)
}
subsett("PTEN", 20)

#save plot
png("PTENTOP20network.png", width = 25, height = 25, unit = "cm", res = 600)
subsett("PTEN", 20)
dev.off()

#sx <- se[grep("PTEN", se$ProteinCodingGene), ]
#sx <- sx[order(sx$Weight, decreasing = TRUE), ]
#sx <- sx[1:10, ]
#sx

####pseudogene only nodes connected by miRNA edges network 
edges <- network[, c("Pseudogene", "miRNA")] #225,898 rows 
edges <- unique(edges)
#merge() 
?graph_from_data_frame

##################################Circilize plot for top 50 network

#saveplot
png("top50Circos.png", width = 25, height = 25, unit = "cm", res = 600)
?png()

#generate a diagram with empty titles
chordDiagram(link1, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))),
             big.gap = 10)

#highlighting the different sectors -> has to be before adding the titles 
highlight.sector(link1$ProteinCodingGene, track.index = 1, col = "gold", 
                 text = "ProteinCodingGene", cex = 1, text.col = "black", niceFacing = TRUE,
                 padding = c(-.8, 0, -.3, 0), text.vjust = -2.5)                               #bottom, x, top, x 

highlight.sector(link1$Pseudogene, track.index = 1, col = "skyblue",
                 text = "Pseudogenes", cex = 1, text.col = "black", niceFacing = TRUE,
                 padding = c(-.8, 0, -.3, 0), text.vjust = -2.5)

#then add titles at an angle 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

#save plot
dev.off()

#par(cex = 0.5, mar = c(0, 0, 0, 0))    #text size 
circos.clear()

#####################
###############################
#########################################Potential figures for paper

########################dataneeded -> TRANSCRIBED pseudogene
columns  #contains the pre-post patient information
transcribedPseudogenes <- pseudogenes[grepl("transcribed", pseudogenes$Type, ignore.case = T), ] #1633 rows 
transcribedPseudogenesCounts <- counts[rownames(counts) %in% transcribedPseudogenes$Name, ] #669 rows

# pre/post-NAC expression -> for pseudogene
patientPseudogeneExpressionChanges1 <- lapply(
  sort(unique(columns$Patient)), function(patient) {
    lapply(sort(rownames(transcribedPseudogenesCounts)), function(pseudogene) {
      list(Patient = patient, Pseudogene = pseudogene, 
           Pre = transcribedPseudogenesCounts[pseudogene, rownames(columns)[columns$Patient == patient & columns$Treatment == "Pre"]],
           Post = transcribedPseudogenesCounts[pseudogene, rownames(columns)[columns$Patient == patient & columns$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(abs(Post - Pre)))

# devtools::install_github("an-bui/calecopal")
colours <- cal_palette("desert", max(round(patientPseudogeneExpressionChanges1$Delta)) + 2, type = "continuous")
png("Pseudogene expression changes.png", res = 600, width = 25, height = 15, units = "cm")
par(mfrow = c(1, 2))
plot(patientPseudogeneExpressionChanges1$Pre, patientPseudogeneExpressionChanges1$Post, 
     col = colours[as.factor(round(patientPseudogeneExpressionChanges1$Delta))], 
     cex = ifelse(patientPseudogeneExpressionChanges1$Pseudogene == "FER1L4", 1.2, 1), pch = ifelse(patientPseudogeneExpressionChanges1$Pseudogene == "FER1L4", 17, 20), 
     xlim = c(0, max(transcribedPseudogenesCounts)), ylim = c(0, max(transcribedPseudogenesCounts)),
     main = "Pseudogene expression changes after NAC")
legend("topleft", legend = 0:max(round(patientPseudogeneExpressionChanges1$Delta)), fill = colours, 
       title = "log2(abs(Post - Pre))", ncol = 3)
legend("topright", legend = c("FER1L4", "Other"), pch = c(17, 20), title = "Pseudogene")
plot(patientPseudogeneExpressionChanges1$Pre, patientPseudogeneExpressionChanges1$Post, 
     col = colours[as.factor(round(patientPseudogeneExpressionChanges1$Delta))], 
     cex = ifelse(patientPseudogeneExpressionChanges1$Pseudogene == "FER1L4", 1.2, 1), pch = ifelse(patientPseudogeneExpressionChanges1$Pseudogene == "FER1L4", 17, 20), 
     xlim = c(0, 5000), ylim = c(0, 5000),
     main = "Pseudogene expression changes after NAC\n(zoomed)")
dev.off()

# pre/post-NAC expression
patientProteinExpressionChanges1 <- lapply(
  sort(unique(columns$Patient)), function(patient) {
    lapply(sort(rownames(genecounts)), function(protein) {
      list(Patient = patient, ProteinCodingGene = protein, 
           Pre = genecounts[protein, rownames(columns)[columns$Patient == patient & columns$Treatment == "Pre"]],
           Post = genecounts[protein, rownames(columns)[columns$Patient == patient & columns$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(abs(Post - Pre)))

#box plot 
patientPseudogeneExpressionChanges1
box_plot1 <- patientPseudogeneExpressionChanges1[patientPseudogeneExpressionChanges1$Pseudogene %in% c('FER1L4','TXLNGY','ANKRD20A11P'), ] 
box_plot1 <- subset(box_plot1, select = -Delta) #96 rows 
box_plot1 <- melt(box_plot1, id.vars = c("Patient","Pseudogene"))    #96X2 192 rows to put preandpost in same column as required by ggplot2

box_plott1 <- ggplot(box_plot1, 
                    aes(x=factor(variable),           #x-axis value
                        y=value,                      #y-axis value 
                        fill=factor(variable))) +     #what the individual box plot is representing
  geom_boxplot() +                                    
  facet_wrap(~Pseudogene)                             #What is separating boxplots from each other on the x-axis

box_plott1

#save
png("top3boxplotall.png", width = 25, height = 25, unit = "cm", res = 600)
box_plott1         
dev.off()

#t-test
BP_FER1L4 <- patientPseudogeneExpressionChanges1[patientPseudogeneExpressionChanges1$Pseudogene == 'FER1L4', ]
t.test(BP_FER1L4$Pre, BP_FER1L4$Post, paired = TRUE)

#SMPD4BP
BP_TXLNGY <- patientPseudogeneExpressionChanges1[patientPseudogeneExpressionChanges1$Pseudogene == 'TXLNGY', ]
t.test(BP_TXLNGY$Pre, BP_TXLNGY$Post, paired = TRUE)

#RPLP0P2
BP_ANKRD20A11P <- patientPseudogeneExpressionChanges1[patientPseudogeneExpressionChanges1$Pseudogene == 'ANKRD20A11P', ]
t.test(BP_ANKRD20A11P$Pre, BP_ANKRD20A11P$Post, paried = TRUE)

##############
#######################
#######################################Ways to summarise the entire network 
link3 <- network     #the overall network 
link3 <- networklink3 <- unique(link3)
length(unique(link3$ProteinCodingGene))
length(unique(link3$miRNA))

##########histogram for miRNA:
link4 <- CodingPseudo       #1,646,631
link4 <- link4 %>%                           #Add count number to duplicate pairs 
  add_count(across(everything()))
link4 <- unique(link4)     #817, 447 matches
colnames(link4) <- c('ProteinCodingGene','Pseudogene','Weight')

# Plot a histogram of miRNA frequencies
options(scipen=999)            #get rid of the e notations

png("miRNADistrribution.png", width = 25, height = 25, units = "cm", res = 600)
hist(link4$Weight,
     main = "Distribution of miRNAs in the network",
     xlab = "miRNA Match Number",
     ylab = "Frequency",
     col = "grey",
     border = "black")
dev.off()

mean(link4$Weight)
max(link4$Weight) #76
min(link4$Weight) #1
median(link4$Weight)


####################
###################################
#################################################Correlation analysis
#Only have pseudogenes and protein coding gene expression data -> so can do correlation analysis between predicted genes. 
#see if the correlation went the same way depending on if they got better or worse 

#need to see whether treatment better or not -1 better. +1 worse. 0 no difference 
#######staging 
columns1 <- columns
columns1$rownames <- rownames(columns1)
columns1$stage <- strsplit(as.character(columns1$Barcode), "")
columns1$stage <- sapply(columns1$stage, `[`, 9) #get the gene names 

#difference in pre and post staging
treatmentScoresPre <- columns1[columns1$Treatment == 'Pre', ]
treatmentScoresPre <- treatmentScoresPre[order(treatmentScoresPre$Patient), ]

treatmentScoresPost <- columns1[columns1$Treatment == 'Post', ]
treatmentScoresPost <- treatmentScoresPost[order(treatmentScoresPost$Patient), ]

deltaScores <- as.numeric(treatmentScoresPost$stage) - as.numeric(treatmentScoresPre$stage)

#merge ordered patient names and ordered deltascores together to original data frame 
jf <- cbind(treatmentScoresPre$Patient, deltaScores)
colnames(jf) <- c('Patient', 'deltaScores')
columns1 <- merge(columns1, jf, by = "Patient")

###looks like they all got better?


#####################
##############################
###############################################Correlation analysis

#the top 50 matches that were predicted and what we are trying to match the patient data to

#first obtain final merged table with delta scores for pre and post for every patient for every match 
corr <- link1[1:50,]
corr$Weight <- NULL

#this is the protein changes for proteins for all patients in corr top 50 matches 
#subset the genecounts to top 50 matches protein network 
CorrProtCounts <- genecounts[rownames(genecounts) %in% corr$ProteinCodingGene, ]

CorrProtTargetChanges <- lapply(
  sort(unique(columns$Patient)), function(patient) {
    lapply(sort(rownames(CorrProtCounts)), function(protein) {
      list(Patient = patient, Protein = protein, 
           Pre = CorrProtCounts[protein, rownames(columns)[columns$Patient == patient & columns$Treatment == "Pre"]],
           Post = CorrProtCounts[protein, rownames(columns)[columns$Patient == patient & columns$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))         #128 rows #to maintain the sign change for later on

#this is the pseudogene changes for pseudogenes for all patients in corr top 50 matches
#subset the pseudogene counts to top 50 matches pseudogene network 
CorrPseudoCounts <- transcribedPseudogenesCounts[rownames(transcribedPseudogenesCounts) %in% corr$Pseudogene, ]

CorrPseudoTargetChanges <- lapply(
  sort(unique(columns$Patient)), function(patient) {
    lapply(sort(rownames(CorrPseudoCounts)), function(pseudogene) {
      list(Patient = patient, Pseudogene = pseudogene, 
           Pre = CorrPseudoCounts[pseudogene, rownames(columns)[columns$Patient == patient & columns$Treatment == "Pre"]],
           Post = CorrPseudoCounts[pseudogene, rownames(columns)[columns$Patient == patient & columns$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))      #1,376 rows 

#merging -> contains the pre and post for every pateint for the top 50 match
finalmerge <- merge(corr, CorrProtTargetChanges, by.x = "ProteinCodingGene", by.y = "Protein")               #1,600 rows
finalmerge1 <- merge(finalmerge, CorrPseudoTargetChanges, by = c("Pseudogene","Patient"))                    #1,600 rows 
colnames(finalmerge1) <- c("Pseudogene","Patient","ProteinCodingGene","ProteinPre","ProteinPost",
                           "ProteinDelta","PseudoPre","PseudoPost","PseudoDelta")
finalmerge1 <- finalmerge1[ ,c("Patient","ProteinCodingGene","ProteinPre","ProteinPost", 
                              "ProteinDelta","Pseudogene","PseudoPre","PseudoPost","PseudoDelta")]

###calculating the correlation 
finalmerge1 <- finalmerge1[order(finalmerge1$ProteinCodingGene, finalmerge1$Pseudogene, finalmerge1$Patient), ]
rownames(finalmerge1) <- NULL

n <- 32
PearCor <- data.frame(
  Protein = unique(finalmerge1[, c("ProteinCodingGene","Pseudogene")]),
  ProteinMean = aggregate(finalmerge1$ProteinDelta, list(rep(1:(nrow(finalmerge1) %/% n + 1),
                                                             each = n, len = nrow(finalmerge1))), mean)[-1],
  PseudoMean = aggregate(finalmerge1$PseudoDelta, list(rep(1:(nrow(finalmerge1) %/% n + 1),
                                                               each = n, len = nrow(finalmerge1))), mean)[-1],
  Correlation = by(finalmerge1[c('ProteinDelta', 'PseudoDelta')],
                   as.integer(gl(nrow(finalmerge1), 32, nrow(finalmerge1))), FUN = function(x) cor(x$ProteinDelta, x$PseudoDelta)),
  pValue = by(finalmerge1[c('ProteinDelta', 'PseudoDelta')],
              as.integer(gl(nrow(finalmerge1), 32, nrow(finalmerge1))), FUN = function(x) cor.test(x$ProteinDelta, x$PseudoDelta)$p.value)
  )

colnames(PearCor) <- c("ProteinCodingGene","Pseudogene","ProteinMean","PseudoMean", "Correlation", "pValue")


# Check -> matches for FER1L4
k <- c(0.81984587, -1.56429686, 0.62639829, 1.00548554, -0.35633139, -0.49057013, -1.24997825, -0.52645878, -1.46858697, 0.15832082, 
       1.86649841, 0.30873445, 0.81323149, -1.47393119, 0.72531915, -1.05889369, -0.76921979, -0.27499351, 0.03613379, 
       -0.10299399, 0.55053867, -0.73535274, -0.06127014, -0.29081197, -0.42247261, 0.44814581, 1.31410859, 0.40186211, 
       -0.79822768, 0.64672567, 0.17483380, 0.40637570)

l <- c(-0.13532634, -0.25940384, 0.39317316, -0.01868760, 1.31937239, -0.03329444, 0.38187064, 1.21055728, -1.29299055, 0.18982456,
       1.12212331, 2.30054770, 0.55483107, -2.57778847, 1.42146377, -0.26662543, 4.13845821, 2.84245872, -0.45852054, -0.59621896, 
       -0.95808636, -0.31757169, 6.66615571, 0.45567948, -2.91541467, 1.43666598, 0.63640194, 0.26303441, 5.67796140, 1.13038893, 
       -0.77741334, -1.64550404)

# Calculate the correlation -> yeah matches what the data frame gave me
correlation <- cor(k, l)
correlation

#################
########################
##############################Correlatioon for the whole predicted network matches
#separate patient data into good and bad outcomes
#and then find the correlation between protein and pseudogene for the entire predicted matches network

#overall patient Data 
columns2 <- columns

#outcome data
outcome <- read.delim('Samples.tsv', sep="\t")                       #read in outcome data
outcome <- outcome[order(outcome$Patient, outcome$Treatment), ]      #order based on patient number and then treatment
columns2 <- columns2[order(columns2$Patient, columns2$Treatment), ]  #order so its same format as outcome
row.names(outcome) <- row.names(columns2)                            #add missing sample name to outcome data

#separate patient data into good and bad outcomes 
outcomeB <- outcome[outcome$TGroup == '>=T2', ]        #50 
outcomeG <- outcome[outcome$TGroup == '<T2', ]         #7

#predicted matches -> add miRNA numbers to it
link2 <- CodingPseudo       #1,646,631
link2 <- link2 %>%                           #Add count number to duplicate pairs 
  add_count(across(everything()))
link2 <- unique(link2)     #817, 447 matches
colnames(link2) <- c('ProteinCodingGene','Pseudogene','Weight')
#then network with predicted matches to use for correlation
corr <- link2
corr$Weight <- NULL     #817,446 entries

#subset the protein genecounts to matches protein network 
CorrProtCounts <- genecounts[rownames(genecounts) %in% corr$ProteinCodingGene, ]
#subset the pseudogene genecounts to matches pseudogene network 
CorrPseudoCounts <- transcribedPseudogenesCounts[rownames(transcribedPseudogenesCounts) %in% corr$Pseudogene, ]

#bad outcome
#Prot expression delta changes for Bad outcome patient list 
CorrProtTargetChangesB <- lapply(
  sort(unique(outcomeB$Patient)), function(patient) {
    lapply(sort(rownames(CorrProtCounts)), function(protein) {
      list(Patient = patient, Protein = protein, 
           Pre = CorrProtCounts[protein, rownames(outcomeB)[outcomeB$Patient == patient & outcomeB$Treatment == "Pre"]],
           Post = CorrProtCounts[protein, rownames(outcomeB)[outcomeB$Patient == patient & outcomeB$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))   #66,150 entries

#this is the pseudogene changes for whole matched network for bad outcome patients
#pseudo delta change for Bad outcome patient list 
CorrPseudoTargetChangesB <- lapply(
  sort(unique(outcomeB$Patient)), function(patient) {
    lapply(sort(rownames(CorrPseudoCounts)), function(pseudogene) {
      list(Patient = patient, Pseudogene = pseudogene, 
           Pre = CorrPseudoCounts[pseudogene, rownames(outcomeB)[outcomeB$Patient == patient & outcomeB$Treatment == "Pre"]],
           Post = CorrPseudoCounts[pseudogene, rownames(outcomeB)[outcomeB$Patient == patient & outcomeB$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))   #16,725 entries


#########good outcome

#Prot Changes for Good outcome patient list 
CorrProtTargetChangesG <- lapply(
  sort(unique(outcomeG$Patient)), function(patient) {
    lapply(sort(rownames(CorrProtCounts)), function(protein) {
      list(Patient = patient, Protein = protein, 
           Pre = CorrProtCounts[protein, rownames(outcomeG)[outcomeG$Patient == patient & outcomeG$Treatment == "Pre"]],
           Post = CorrProtCounts[protein, rownames(outcomeG)[outcomeG$Patient == patient & outcomeG$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))  #18, 000 entries 

#pseudo changes for Good outcome patient list 
CorrPseudoTargetChangesG <- lapply(
  sort(unique(outcomeG$Patient)), function(patient) {
    lapply(sort(rownames(CorrPseudoCounts)), function(pseudogene) {
      list(Patient = patient, Pseudogene = pseudogene, 
           Pre = CorrPseudoCounts[pseudogene, rownames(outcomeG)[outcomeG$Patient == patient & outcomeG$Treatment == "Pre"]],
           Post = CorrPseudoCounts[pseudogene, rownames(outcomeG)[outcomeG$Patient == patient & outcomeG$Treatment == "Post"]])
    })
  }) %>% 
  bind_rows() %>%
  mutate(Delta = log2(Pre + 1) - log2(Post + 1))    #4, 683 entries


#######################merging -> (Bad outcome)  protein and pseudogene correlation 
#merge matches network with counts for protein changes

#convert to data table for faster merging time 
library(data.table)
setDT(corr)
setDT(CorrProtTargetChangesB)
setDT(CorrPseudoTargetChangesB)
setDT(CorrProtTargetChangesG)
setDT(CorrPseudoTargetChangesG)


mergeB <- merge(corr, CorrProtTargetChangesB, by.x = "ProteinCodingGene", by.y = "Protein", allow.cartesian=TRUE)         
mergeB <- merge(mergeB, CorrPseudoTargetChangesB, by = c("Pseudogene","Patient"), allow.cartesian=TRUE) 


colnames(mergeB) <- c("Pseudogene","Patient","ProteinCodingGene","ProteinPre","ProteinPost",
                          "ProteinDelta","PseudoPre","PseudoPost","PseudoDelta")
mergeB <- mergeB[ ,c("Patient","ProteinCodingGene","ProteinPre","ProteinPost",
                     "ProteinDelta","Pseudogene","PseudoPre","PseudoPost","PseudoDelta")]    #20,436,150 entries

#need to order it first before calculating so that the number is correct
mergeB <- mergeB[order(mergeB$ProteinCodingGene, mergeB$Pseudogene, mergeB$Patient), ]   #20,436,150 entries
rownames(mergeB) <- NULL

#the below function doesn't like data.table so convert back into dataframe 
mergeB <- as.data.frame(mergeB)

n <- 25    #change repetition if needed for real outcome data
PearCorB <- data.frame(
  Protein = unique(mergeB[, c("ProteinCodingGene","Pseudogene")]),
  ProteinMean = aggregate(mergeB$ProteinDelta, list(rep(1:(nrow(mergeB) %/% n + 1),
                                                            each = n, len = nrow(mergeB))), mean)[-1],
  PseudoMean = aggregate(mergeB$PseudoDelta, list(rep(1:(nrow(mergeB) %/% n + 1),
                                                          each = n, len = nrow(mergeB))), mean)[-1],
  Correlation = by(mergeB[c('ProteinDelta', 'PseudoDelta')],
                   as.integer(gl(nrow(mergeB), 25, nrow(mergeB))), FUN = function(x) cor(x$ProteinDelta, x$PseudoDelta)),
  pValue = by(mergeB[c('ProteinDelta', 'PseudoDelta')],
                 as.integer(gl(nrow(mergeB), 25, nrow(mergeB))), FUN = function(x) cor.test(x$ProteinDelta, x$PseudoDelta)$p.value)
  )

colnames(PearCorB) <- c("ProteinCodingGene","Pseudogene","ProteinMean","PseudoMean", "Correlation","pValue")     #817,446 entries

##############multiple testing:
# Extract p-values
p_values <- PearCorB$pValue

# Adjust p-values using Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to PearCorB data frame
PearCorB$Adjusted_pValue <- adjusted_p_values

#View(PearCorB)     #817,446 entries

############################merging -> Good outcome  protein and pseudogene correlation 
mergeG <- merge(corr, CorrProtTargetChangesG, by.x = "ProteinCodingGene", by.y = "Protein", allow.cartesian=TRUE)         
mergeG <- merge(mergeG, CorrPseudoTargetChangesG, by = c("Pseudogene","Patient"), allow.cartesian=TRUE)

colnames(mergeG) <- c("Pseudogene","Patient","ProteinCodingGene","ProteinPre","ProteinPost",
                      "ProteinDelta","PseudoPre","PseudoPost","PseudoDelta")
mergeG <- mergeG[ ,c("Patient","ProteinCodingGene","ProteinPre","ProteinPost", 
                     "ProteinDelta","Pseudogene","PseudoPre","PseudoPost","PseudoDelta")]  

mergeG <- mergeG[order(mergeG$ProteinCodingGene, mergeG$Pseudogene, mergeG$Patient), ]   
rownames(mergeG) <- NULL               #5,722,122 entries

##the below function doesn't like data.table so convert back into dataframe 
mergeG <- as.data.frame(mergeG)

#finding the correlation 
n <- 7
PearCorG <- data.frame(
  Protein = unique(mergeG[, c("ProteinCodingGene","Pseudogene")]),
  ProteinMean = aggregate(mergeG$ProteinDelta, list(rep(1:(nrow(mergeG) %/% n + 1),
                                                            each = n, len = nrow(mergeG))), mean)[-1],
  PseudoMean = aggregate(mergeG$PseudoDelta, list(rep(1:(nrow(mergeG) %/% n + 1),
                                                          each = n, len = nrow(mergeG))), mean)[-1],
  Correlation = by(mergeG[c('ProteinDelta', 'PseudoDelta')],
                   as.integer(gl(nrow(mergeG), 7, nrow(mergeG))), FUN = function(x) cor(x$ProteinDelta, x$PseudoDelta)),
  pValue = by(mergeG[c('ProteinDelta', 'PseudoDelta')],
              as.integer(gl(nrow(mergeG), 7, nrow(mergeG))), FUN = function(x) cor.test(x$ProteinDelta, x$PseudoDelta)$p.value)
  )
colnames(PearCorG) <- c("ProteinCodingGene","Pseudogene","ProteinMean","PseudoMean", "Correlation","pValue")   #817,446 entries


# Extract p-values
p_values1 <- PearCorG$pValue

# Adjust p-values using Benjamini-Hochberg method
adjusted_p_values1 <- p.adjust(p_values1, method = "BH")

# Add adjusted p-values to PearCorB data frame
PearCorG$Adjusted_pValue <- adjusted_p_values1


#############################
###keep only >0.05 p-adjusted values for the correlation table
PearCorB <- PearCorB[PearCorB$Adjusted_pValue < 0.05, ]        # 817,446 entries to #4,017 entries
PearCorB <- PearCorB[!is.na(PearCorB$Adjusted_pValue), ]

PearCorG <- PearCorG[PearCorG$Adjusted_pValue < 0.05, ]        #4 entries
PearCorG <- PearCorG[!is.na(PearCorG$Adjusted_pValue), ]       #after 38 entries???
                                                               #195 entries 

#################################
######################################adding in the miRNA numbers

#turn into data.table for faster merging
setDT(PearCorB)
setDT(PearCorG)

PearCorB1 <- merge(PearCorB, link2, by = c('ProteinCodingGene','Pseudogene'))

PearCorG1 <- merge(PearCorG, link2, by = c('ProteinCodingGene','Pseudogene'))

#############################colouring the table 

####ProteinCodingGene cancer type colourway

#COSMIC Cancer Gene Census file 
cgc <- read.delim('Cosmic_CancerGeneCensus_v99_GRCh38.tsv', sep="\t")   #743 entries 

# our genes dataset is filtered to chromosomes 1:22/X/Y but some of the CGC are on alternative scaffolds 
# so filter those out since they'll never match our other analyses
cgc <- cgc[cgc$GENE_SYMBOL %in% genes$Name, ]       #742 entries

# expand roles in cancer
cgc$Oncogene <- grepl("oncogene", cgc$ROLE_IN_CANCER)
cgc$TumourSuppressor <- grepl("TSG", cgc$ROLE_IN_CANCER)
cgc$Fusion <- grepl("fusion", cgc$ROLE_IN_CANCER)

# convert effectively logical fields to logical
#cgc$Hallmark <- grepl("Yes", cgc$hallmark)
cgc$Somatic <- grepl("yes", cgc$SOMATIC)
cgc$Germline <- grepl("yes", cgc$GERMLINE)
cgc$OtherGermlineMut <- grepl("yes", cgc$OTHER_GERMLINE_MUT)

#merge this with the table and just drop the unnecessary columns

#concert to data.table
setDT(cgc)
setDT(PearCorB1)
setDT(PearCorG1)

#merge
#keep unmatches rows too because not everything will be a cancer gene obvi
PearCorB11 <- merge(PearCorB1, cgc, by.x = 'ProteinCodingGene', by.y = 'GENE_SYMBOL', all.x = TRUE)   
PearCorG11 <- merge(PearCorG1, cgc, by.x = 'ProteinCodingGene', by.y = 'GENE_SYMBOL', all.x = TRUE)

#format
PearCorB11 <- PearCorB11[ , c('ProteinCodingGene', 'Oncogene', 'TumourSuppressor', 'Fusion', 'Pseudogene', 'ProteinMean',
                              'PseudoMean', 'Correlation','pValue','Adjusted_pValue', 'Weight')]
PearCorG11 <- PearCorG11[ , c('ProteinCodingGene', 'Oncogene', 'TumourSuppressor', 'Fusion', 'Pseudogene', 'ProteinMean',
                              'PseudoMean', 'Correlation','pValue','Adjusted_pValue', 'Weight')]

#write.csv(PearCorB11, "PearCorB11.csv")
#write.csv(PearCorG11, "PearCorG11.csv")



######################
##################################
################################################Heatmap rendition three

######Heatmap for the top 3 Delta pseudogene 

#Heatmap with only response 
posterHeatmap <- function(treatment = c("Pre", "Post"), showPatientNames = T)
{
  png(sprintf("%s-top3.png", treatment), width = 45, height = 23, units = "cm", res = 600)
  
  #heatmap matrix division
  heatmapMatrix <- transcribedPseudogenesCounts[rownames(transcribedPseudogenesCounts) %in% c("FER1L4", "TXLNGY","ANKRD20A11P"), ] #contains only top 3 pseudogenes
  heatmapMatrix <- heatmapMatrix[colnames(heatmapMatrix) %in% rownames(outcome[outcome$Treatment == treatment, ])]
  
  #sample outcome Division
  samples <- outcome[outcome$Treatment == treatment, ]
  samples$TGroup <- factor(samples$TGroup)
  
  #top annotation
  ta <- HeatmapAnnotation(
    Response = samples$TGroup,
    col = list(
      Response = setNames(brewer.pal(3, "Set1")[c(1,3)], levels(samples$TGroup))
    ),
    show_legend = c(T),
    annotation_legend_param = list(
      Response = list(
        labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold"))
    )
  )
  
  #split according to response
  columnSplit <- samples$TGroup
  
  # show patients instead of samples
  colnames(heatmapMatrix) <- samples$Patient[match(colnames(heatmapMatrix), rownames(samples))]
  
  
  #heatmap
  set.seed(123)
  print(draw(ComplexHeatmap::Heatmap((scale(heatmapMatrix, center = T, scale = T)), name = "Z-score", 
                                     column_title = sprintf("%s-Top 3 Expression Scores", treatment), column_title_side = "top", column_title_gp = gpar(fontsize = 20),
                                     cluster_columns = T, cluster_rows = T,
                                     width = unit(30, "cm"), height = unit(16, "cm"),
                                     top_annotation = ta,
                                     column_split = columnSplit,
                                     heatmap_legend_param = list(
                                       labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
                                     )), 
             merge_legend = T))
  
  dev.off()
}

posterHeatmap("Pre")
posterHeatmap("Post")


####################################

#heatmap for top 50 matches protein coding and pseudogenes. 

#contain the top 50 matches
top50 <- link1[1:50,]

#subset the expression matrix to contain only these genes
View(counts)

#patient data
View(outcome)
outcome <- outcome[order(outcome$Patient), ]

#heatmap with only outcome
posterHeatmap1 <- function(treatment = c("Pre", "Post"), showPatientNames = T)
{
  png(sprintf("%s-top50.png", treatment), width = 45, height = 23, units = "cm", res = 600)
  
  #heatmap matrix count division 
  heatmapMatrix <- counts[rownames(counts) %in% top50$ProteinCodingGene, ]
  heatmapMatrix1 <- counts[rownames(counts) %in% top50$Pseudogene, ]
  heatmapMatrix <- rbind(heatmapMatrix, heatmapMatrix1)
  
  #log transformation
  heatmapMatrix <- heatmapMatrix + 1
  heatmapMatrix <- log2(heatmapMatrix)
  
  #heatmap matrix treatment division into pre and post
  heatmapMatrix <- heatmapMatrix[colnames(heatmapMatrix) %in% rownames(outcome[outcome$Treatment == treatment, ])]
  
  #sample outcome Division
  samples <- outcome[outcome$Treatment == treatment, ]
  samples$TGroup <- factor(samples$TGroup)
  
  #top annotation
  ta <- HeatmapAnnotation(
    Response = samples$TGroup,
    col = list(
      Response = setNames(brewer.pal(3, "Set1")[c(3,1)], levels(samples$TGroup))
    ),
    show_legend = c(T),
    annotation_legend_param = list(
      Response = list(
        labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold"))
    )
  )
  
  #split according to response
  columnSplit <- factor(samples$TGroup, levels = c('>=T2','<T2'))
  
  # show patients instead of samples
  colnames(heatmapMatrix) <- samples$Patient[match(colnames(heatmapMatrix), rownames(samples))]
  
  #heatmap
  set.seed(123)
  print(draw(ComplexHeatmap::Heatmap((scale(heatmapMatrix, center = T, scale = T)), name = "Z-score", 
                                     column_title = sprintf("%s-Top 50 Expression Scores", treatment), column_title_side = "top", column_title_gp = gpar(fontsize = 20),
                                     cluster_columns = T, cluster_rows = T,
                                     width = unit(30, "cm"), height = unit(16, "cm"),
                                     top_annotation = ta,
                                     column_split = columnSplit,
                                     heatmap_legend_param = list(
                                       labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
                                     )), 
             merge_legend = T))
  
  dev.off()
}

posterHeatmap1("Pre")
posterHeatmap1("Post")


####################
#######################################
###########################################################

#############heat map for the expression changes 

CorrProtTargetChanges
rownames(CorrProtTargetChanges) <- CorrProtTargetChanges$Protein
CorrPseudoTargetChanges

posterHeatmap2 <- function(treatment = c("Pre", "Post"), showPatientNames = T)
{
  png(sprintf("%s-top50.png", treatment), width = 45, height = 23, units = "cm", res = 600)
  
  #change matrix for protein coding genes 
  
  
  #heatmap matrix count division 
  heatmapMatrix <- counts[rownames(counts) %in% top50$ProteinCodingGene, ]
  heatmapMatrix1 <- counts[rownames(counts) %in% top50$Pseudogene, ]
  heatmapMatrix <- rbind(heatmapMatrix, heatmapMatrix1)
  
  #log transformation
  heatmapMatrix <- heatmapMatrix + 1
  heatmapMatrix <- log2(heatmapMatrix)
  
  #heatmap matrix treatment division into pre and post
  heatmapMatrix <- heatmapMatrix[colnames(heatmapMatrix) %in% rownames(outcome[outcome$Treatment == treatment, ])]
  
  #sample outcome Division
  samples <- outcome[outcome$Treatment == treatment, ]
  samples$TGroup <- factor(samples$TGroup)
  
  #top annotation
  ta <- HeatmapAnnotation(
    Response = samples$TGroup,
    col = list(
      Response = setNames(brewer.pal(3, "Set1")[c(1,3)], levels(samples$TGroup))
    ),
    show_legend = c(T),
    annotation_legend_param = list(
      Response = list(
        labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold"))
    )
  )
  
  #split according to response
  columnSplit <- factor(samples$TGroup, levels = c('>=T2','<T2'))
  
  # show patients instead of samples
  colnames(heatmapMatrix) <- samples$Patient[match(colnames(heatmapMatrix), rownames(samples))]
  
  #heatmap
  set.seed(123)
  print(draw(ComplexHeatmap::Heatmap((scale(heatmapMatrix, center = T, scale = T)), name = "Z-score", 
                                     column_title = sprintf("%s-Top 50 Expression Scores", treatment), column_title_side = "top", column_title_gp = gpar(fontsize = 20),
                                     cluster_columns = T, cluster_rows = T,
                                     width = unit(30, "cm"), height = unit(16, "cm"),
                                     top_annotation = ta,
                                     column_split = columnSplit,
                                     heatmap_legend_param = list(
                                       labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
                                     )), 
             merge_legend = T))
  
  dev.off()
}



posterHeatmap2 <- function()
{
  png("change-in-expression-levels.png", width = 45, height = 23, units = "cm", res = 600)
  
  # Calculate change in expression levels for protein-coding genes
  changeMatrixProtein <- counts[rownames(counts) %in% top50$ProteinCodingGene, "Post"] - counts[rownames(counts) %in% top50$ProteinCodingGene, "Pre"]
  
  # Calculate change in expression levels for pseudogenes
  changeMatrixPseudogene <- counts[rownames(counts) %in% top50$Pseudogene, "Post"] - counts[rownames(counts) %in% top50$Pseudogene, "Pre"]
  
  # Combine change matrices
  changeMatrix <- rbind(changeMatrixProtein, changeMatrixPseudogene)
  
  # Log transformation
  changeMatrix <- log2(changeMatrix + 1)
  
  # Sample outcome division
  samples <- outcome[outcome$Treatment %in% c("Pre", "Post"), ]
  samples$TGroup <- factor(samples$TGroup)
  
  # Top annotation
  ta <- HeatmapAnnotation(
    Response = samples$TGroup,
    col = list(
      Response = setNames(brewer.pal(3, "Set1")[c(1,3)], levels(samples$TGroup))
    ),
    show_legend = c(T),
    annotation_legend_param = list(
      Response = list(
        labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold"))
    )
  )
  
  # Split according to response
  columnSplit <- factor(samples$TGroup, levels = c('>=T2','<T2'))
  
  # Show patients instead of samples
  colnames(changeMatrix) <- samples$Patient[match(colnames(changeMatrix), rownames(samples))]
  
  # Heatmap
  set.seed(123)
  print(draw(ComplexHeatmap::Heatmap((scale(changeMatrix, center = T, scale = T)), name = "Z-score", 
                                     column_title = "Change in Expression Levels", column_title_side = "top", column_title_gp = gpar(fontsize = 20),
                                     cluster_columns = T, cluster_rows = T,
                                     width = unit(30, "cm"), height = unit(16, "cm"),
                                     top_annotation = ta,
                                     column_split = columnSplit,
                                     heatmap_legend_param = list(
                                       labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
                                     )), 
             merge_legend = T))
  
  dev.off()
}

posterHeatmap2()










######T-cell filtrated or not 

#T_Inf <- data.frame(
#  Patient = c('W005862','H002426','H000021','H001574','W000294','H001091','H001856','W000425','H001637','H001147','H001436','C003582','H002104',
#              'C000495','W000799','C000225','H000198','W000654','H001750','W000323','W006091','H000116','H000902','W001005','W006176',
#              'C000727','H000337','W006206','C000537','H000945','H002037','C005923'),
#  T_Inf = c('Y','Y','Y','Y','Y','N','Y','Y','Y','N','Y','Y','Y','N','N','Y','N','N','N','N','N','Y','Y','N','N','N','N','N','N','N','N','N')
#)

#merge the T_nf with columns3 and get pre and post 
#columns3 <- merge(columns3, T_Inf, by = "Patient")
#columns3Pre <- columns3[columns3$Treatment == "Pre", ]
#columns3Post <- columns3[columns3$Treatment == "Post", ]



heatmapMatrix <- counts[rownames(counts) %in% top50$ProteinCodingGene, ]
heatmapMatrix1 <- counts[rownames(counts) %in% top50$Pseudogene, ]
heatmapMatrix <- rbind(heatmapMatrix, heatmapMatrix1)

changeMatrixProtein <- counts[rownames(counts) %in% top50$ProteinCodingGene, "Post"] - counts[rownames(counts) %in% top50$ProteinCodingGene, "Pre"]





