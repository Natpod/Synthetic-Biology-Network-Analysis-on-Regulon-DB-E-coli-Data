#-----------------------------------------------------
#title: 'Synthetic Biology Course Practice: Exercise 1'
#author: "Natalia García Sánchez, 06632584M"
#date: "October 6th, 2022"
#-----------------------------------------------------


# INTRODUCTION

#The present work explores three questions subjacent to
# the properties of the regulatory network interaction between
# **Transcriptional Factors (TF)** and **genes** 
# gathered in RegulonDB[@tierrafría2022] experimental dataset 
# `./network_tf_gene.txt` 
# from experimental evidence *Escherichia coli*.

# OBJECTIVES

#1.  *Distribution of Feed Forward Loop motif frequency in comparison to random networks*

#2.  *Distribution of Autorregulatory motif frequency in comparison to random networks*

#3.  *Longest and shortest cascade length distribution in comparison to random networks*

# DATA PROCESSING

#First we will import the `./network_tf_gene.txt` dataset present in the environment directory

#------------------------------------------------------------------------
library('igraph')
library(showtext)

wt<-read.table("network_tf_gene.txt", sep="\t")
names(wt)<-c("TF_ID", "TF_name", "g_ID", "g_name", "Reg_eff", "ev", "evT")


#The `./network_tf_gene.txt` dataset consists of 8 columns:

#-   TF_ID --\> Transcription Factor (TF) ID

#-   TF_name --\>Transcription Factor (TF) Name

#-   g_ID --\> Gene ID regulated by the TF (regulated gene)

#-   g_name --\> Gene regulated by the TF (regulated gene)

#-   Reg_eff --\> Regulatory effect of the TF on the regulated gene (+ activator, - repressor, +- dual, ? unknown)

#-   ev --\> Evidence that supports the existence of the regulatory interaction

#-   evT --\>Evidence Type [Weak, Strong, Confirmed]


summary(wt)


#Character conversion in TF name column into all undercase letters.


rownum<-nrow(wt)

for (i in 1:rownum){
  wt[i,"TF_name"]<-tolower(wt[i,"TF_name"])
  wt[i,"g_name"]<-tolower(wt[i,"g_name"])
}


## Data exploration

#We will explore the representation of regulatory effects and evidence 
# support levels to define the scope of our analysis

#Check respective evidence levels

rows_regtype <- (wt$evT == "CONFIRMED")
sprintf("Number of instances with confirmed and  evidence support : %s", sum(rows_regtype))

rows_regtype <- (wt$evT == "WEAK")
sprintf("Number of instances with Weak evidence support : %s", sum(rows_regtype))

rows_regtype <- (wt$evT == "STRONG")
sprintf("Number of instances with Strong  evidence support : %s", sum(rows_regtype))


#As there are no cross-validated instances and Weak evidence constitutes 
#a considerable part of the dataset, we will take both evidence 
#instances classes into account for our analysis.

#Check if any of the TF have dual (+-) regulatory effect in any of the genes of the dataset**

rows_regtype <- which(wt$Reg_eff == "+-")
sprintf("Number of instances with dual regulatory effect : %s", sum(wt[rows_regtype, ]))


#0 rows are described by dual regulatory effect (+-), so we will not create the distinction in the 
#adjacency matrix for this regulatory effect


sprintf("Number of instances with positive regulatory effect : %s", sum(wt[,5] == "+"))
sprintf("Number of instances with negative regulatory effect : %s", sum(wt[,5] == "-"))
sprintf("Number of instances with unknown regulatory effect : %s", sum(wt[,5] == "?"))


#Now we will see if there are duplicate TF - gene regulatory tuples


duplicates_specialRE = c(duplicated(wt_PR[, c(2,4)]), duplicated(wt_NR[, c(2,4)]), duplicated(wt_XR[, c(2,4)]))
sprintf("TF-gene Duplicates found for specific regulatory effects :%s", sum(duplicates_specialRE))

duplicates_totalRE = duplicated(wt_TotR[, c(2,4)])
sprintf("TF-gene Duplicates found for total regulatory effects :%s", sum(duplicates_totalRE))


#We will expect adjacency matrices creation with 231 elements with a value over 1 
#in the **Total** regulation adjacency matrices, and matrices of 1s and 0s in the other regulatory effect subset matrices.

## Preprocessing

#Next step is creating a dataframe with logic row indexing of the regulatory effects.
# This will facilitate the creation of the adjacency matrix for each regulatory effect.

#4 different dataframes are created with each regulatory effect and the sum of the 3 former dataframes


wt_PR<- wt[which(wt$Reg_eff == '+'), ]
wt_NR<- wt[which(wt$Reg_eff == '-'), ]
wt_XR<- wt[wt[,5] == "?", ]

wt_TotR<- wt


#Simple network analysis


# listing TOTAL matrix tf/genes
listtf<-sort(unique(wt_TotR[,2]))
listgene<-sort(unique(wt_TotR[,4]))

total_genes<-sort(unique(c(listtf, listgene)))


print("Number of total genes in dataset : ")
length(total_genes)
print("Number of total regulated genes in dataset : ")
length(listgene)
print("Number of total tf genes in dataset : ")
length(listtf)

# ------------------------------------------------------------

# ADJACENCY MATRIX CREATION

#Function to create the adjacency matrix, TF will be represented by rows and genes by columns

#Iterate by tf then index and get genes


# passing df with wtdf$TF_name in index 2 and wtdf_$g_name in index 4

cr_adj_matrix <- function(wtdf = NULL, drop =TRUE){
  
  # adj_matrix using total_genes (sorted gene index') length
  
  adj_matrix<-matrix(0, nrow=length(total_genes), ncol=length(total_genes))
  
  
  # two step iteration : for each row j (representing tf) of adj_matrix:
  
  # indexing genes regulated by the tf
  
  # making a table of gene name occurrencies
  
  # indexing gene in total_gene sorted index to get adj_matrix column j
  
  # assigning wt(tf, gene) count value to adj_matrix[i,j]
  
  for (i in 1:length(total_genes)){
    
    genes_select<- wtdf[which(wtdf[,2] == total_genes[i]),]
    
    gene_list<-table(sort(genes_select$g_name))
    
    for (gene in names(gene_list)){
      
      
      j<-match(gene, total_genes)
      
      adj_matrix[i,j]<-gene_list[gene]
      
      
    }
  }
  
  return(adj_matrix)
  
}

adj_PR<-cr_adj_matrix(wt_PR)
adj_NR<-cr_adj_matrix(wt_NR)
adj_XR<-cr_adj_matrix(wt_XR)
adj_TotR<-cr_adj_matrix(wt_TotR)


### Quantitative analysis of regulation


print("Maximum frequency of Positive Regulation values in adj_matrix")
max(c(adj_PR))
print("Maximum frequency of Negative Regulation values in adj_matrix")
max(c(adj_NR))
print("Maximum frequency of Unknown Regulation values in adj_matrix")
max(c(adj_XR))
print("Maximum frequency of total Regulation values in adj_matrix")
max(c(adj_TotR))
sprintf("Number of values in TotR with accumulated regulatory effects : %s", sum(adj_TotR==2))


#All values in `adj_PR` , `adj_NR`, `adj_XR` take values 1 or 0

#However the Total regulation matrix has the same tf regulating a certain gene. 
#This is due to the fact that this gene is subjected to both positive and negative /unknown 
#regulations alike, and implicit accummulated regulatory effects are seen in `adj_TotR` 
#so 231 values above 1 (2) are observed..

#Now we will add the dimensions names to the matrixes


dimnames(adj_TotR)<-list(total_genes, total_genes)
dimnames(adj_PR)<-list(total_genes, total_genes)
dimnames(adj_NR)<-list(total_genes, total_genes)
dimnames(adj_XR)<-list(total_genes, total_genes)

print("The dimensions of the final adjacency matrices is :")
print(dim(adj_TotR))
print(dim(adj_PR))
print(dim(adj_NR))
print(dim(adj_XR))


# RANDOMIZED NETWORK CREATION

#A function is created using node number and edge number from values in the wild-type 
#adjacency matrices to sample columns randomly and assign corresponding 1 values into a 
#null square matrix with the same dimensions using the `sample` R built-in function. 
#If the argument wild-type adjacency matrix is from the Total Regulation subset, 
#1 or 2 values will be assigned discriminating the 1 value number of edges from 
#the 2 value number of edges.


randomize_network_qualitative<- function(adj_mat, typemat){
  
  num_nodes = dim(adj_mat)[1]
  
  rnd<-matrix(0,num_nodes,num_nodes)
  
  if (typemat %in% c("+", "-", "?")){
    
    # Creates empty rnd with dim = (numnodes, numnodes) and chooses 
	#a n_edge number of edges between indexing numbers (number of elements or genes) with replacement
    
    num_edges = sum(adj_mat==1)
    
    rndrow<-sample(1:num_nodes, num_edges, replace=TRUE) 
    rndcol<-sample(1:num_nodes, num_edges, replace=TRUE)
    rnd[cbind(rndrow,rndcol)]<-1
    
  } else{
    
    # In this case edges can have two values. Calculate edges of 1 oe 2 occurring 
	#and randomly assign the correspondent approximate num_edges expected for each value
    
    num_edges_2 = sum(adj_mat==2)
    num_edges_1 = sum(adj_mat==1)
    
    rndrow_1<-sample(1:num_nodes, num_edges_1, replace=TRUE) 
    rndcol_1<-sample(1:num_nodes, num_edges_1, replace=TRUE)
    rnd[cbind(rndrow_1,rndcol_1)]<-1
    
    rndrow_2<-sample(1:num_nodes, num_edges_2, replace=TRUE) 
    rndcol_2<-sample(1:num_nodes, num_edges_2, replace=TRUE)
    rnd[cbind(rndrow_2,rndcol_2)]<-2
    
  }
  return(rnd)
}


#Most of our adjacency matrixes are binary, so we will expect the values of our randomized network 
#to take values 1 or 0 as well.

#If this were not to happen, an analysis of the distribution of values in the adj_matrixes 
#over 1 would have to be carried out. Then, sampling natural numbers from values 
#1 to `max_interaction_level` using the WT distribution would be done by the 
#following iterative function.

# PLOTTING

### Total WT MATRIX CONVERSION TO 1 AND 0s


adj_TotR_convertedto1 = adj_TotR
adj_TotR_convertedto1[adj_TotR_convertedto1==2] = 1


#Plotting the total subnetworks of the entire network is tricky. 
#We can construct the plottable graph objects from the library igraph with our own adjacency matrices


dimnames(adj_TotR_convertedto1)<-list(total_genes, total_genes)

#reads the adjacency matrix and puts it in a graph object

net<-graph.adjacency(adj_TotR_convertedto1, mode="directed")

graphtot<-graph_from_adjacency_matrix(adj_TotR_convertedto1, mode= "directed")


# We should specify the gene instances in the dataset that are involved 
#in positive regulation and specify the edges involved in that interaction 
#through `positive` vector variable (0 : no positive interaction, 1: positive interaction)

# For randomization, a figure of the number of genes in the dataset are involved
# in positive regulation should be known.


listpos = sort(unique(c(as.vector(wt_PR[,2]), as.vector(wt_PR[,4]))))
positive<-V(graphtot)$name %in% listpos

sprintf("A total of %s genes in the dataset are involved in positive regulation" , sum(positive==1))


# We should specify the gene instances in the dataset that are involved in negative 
#regulation and specify the edges involved in that interaction through `negative` 
#vector variable (0 : no negative interaction, 2: negative interaction)

# For randomization, a figure of the number of genes in the dataset are involved in 
#negative regulation should be known.


listneg = sort(unique(c(as.vector(wt_NR[,2]), as.vector(wt_NR[,4]))))
negative<-V(graphtot)$name %in% listneg + 1

# transform the 1 to 0 for following sum

negative[negative==1]<-0


sprintf("A total of %s genes in the dataset are involved in negative regulation" , sum(negative==2))

# We should specify the gene instances in the dataset that are involved 
#in unknown regulation and specify the edges involved in that interaction
# through `negative` vector variable (0 : no negative interaction, 2: negative interaction)

# For randomization, a figure of the number of genes in the dataset are involved in unknown
# regulation should be known.


listx = sort(unique(c(as.vector(wt_XR[,2]), as.vector(wt_XR[,4]))))
unknown<-V(graphtot)$name %in% listx + 3

# transform the 1 to 0 for following sum
unknown[unknown==3]<-0


sprintf("A total of %s genes in the dataset are involved in unknown regulation" , sum(unknown==4))


# As elements of positive vector are 0 and 1 and negative vector elements
# go from 0 to 2, we can sum up the two former vectors in a
# dual vector +1 
#(1: no specific/known regulation, 2: positive regulation, 3: negative regulation, 4: dual regulation)


total<-positive+negative+unknown
# adds 1 as vectors in r start in 1
total<-total+1


#Star layout to appreciate autorregulation motifs

#Positive and negative regulation will be plotted respectively as 
#red and blue. Negative regulation with unknown regulation and will 
#be plotted as dark blue, similarly to positive with unknown regulation (darkred)

#Unknown regulation will be plotted as black

#Dual positive and negative regulation will be plotted as green

#------

#Graphopt layout


plot(net, 
	vertex.label = NA, 
	vertex.color= c("grey", 
	"tomato2", "royalblue", 
	"yellow3", "black", 
	"darkred", "darkblue","green")[total],
	vertex.size=1, 
	edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total], 
	layout=layout.graphopt, edge.width=1)


#Tree layout

deg <- degree(net, mode="all")

plot(net, 
	vertex.label = NA, 
	vertex.color= c("grey", 
	"tomato2", "royalblue", 
	"yellow3", "black", 
	"darkred", "darkblue","green")[total], 
	vertex.size=round(deg/60), 
	edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total], 
	layout=layout_as_tree, edge.width=1)

#Star layout

plot(net, 
vertex.label = NA, 
vertex.color= c("grey", 
"tomato2", "royalblue", 
"yellow3", "black", 
"darkred", "darkblue",
"green")[total], 
vertex.size=4, 
edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total], 
layout=layout_as_star, edge.width=1)


### Random matrix

#We will use the data for the random regulatory interactions for graph visualization.


rnd_adjmat <- randomize_network_qualitative(adj_TotR_convertedto1, "TOT")
rnd_net<-graph.adjacency(rnd_adjmat, mode="directed")

# generate random regulatory effect data

positive<-vector(length = 1933)
pos_rowsample<-sample(1:1933, 1476, replace=FALSE)
positive[pos_rowsample]<-1

negative<-vector(length = 1933)
neg_rowsample<-sample(1:1933, 1354, replace=FALSE)
negative[neg_rowsample]<-2

unknown<-vector(length = 1933)
neg_rowsample<-sample(1:1933, 23, replace=FALSE)
unknown[neg_rowsample]<-4

total<- positive + negative + unknown +1

#Same plotting total variable with regulatory effect distinction like WT matrix plotting

# Graphopt layout

plot(rnd_net, 
vertex.label = NA, 
vertex.color= c("grey", 
"tomato2", "royalblue", 
"yellow3", "black", 
"darkred", "darkblue","green")[total], 
vertex.size=1, 
edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total]
, layout=layout.graphopt, edge.width=1)


# Tree layout

deg <- degree(rnd_net, mode="all")



plot(rnd_net, 
vertex.label = NA, 
vertex.color= c("grey", "tomato2", 
"royalblue", "yellow3", 
"black", "darkred", "darkblue","green")[total], 
vertex.size=round(deg/60), 
edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total],
 layout=layout_as_tree, edge.width=1)


# Star layout


plot(rnd_net, 
vertex.label = NA, 
vertex.color = c("grey", "tomato2", 
"royalblue", "yellow3", 
"black", "darkred", "darkblue","green")[total], 
vertex.size=4, edge.color=c("grey", "tomato2", "royalblue", "yellow3", "black", "darkred", "darkblue","green")[total],
layout=layout_as_star, edge.width=1)


# 1. Feed Forward Loop (FFL) FREQUENCY

------------------------------------------------------------------------
  
## Total triad count
  
#count_triangles igraph function implements this count


library("igraph")
colnames(adj_TotR)<-total_genes
rownames(adj_TotR)<-total_genes
graphTot<-graph_from_adjacency_matrix(adj_TotR_convertedto1, 
                                      mode="directed",
                                      weighted = NULL,
                                      diag = TRUE,
                                      add.colnames = TRUE,
                                      add.rownames = TRUE)
sprintf("A total number of triangles : %s was found", sum(count_triangles(graphTot)))


## Total FFL

#using triad_census igraph function we will get the specific triads (see documentation : <https://igraph.org/r/doc/triad_census.html>)


triad_classes <- c(
  "A,B,C :: unconnected", 
  "A->B,C :: dual directed", 
  "A<->B,C :: dual mutual", 
  "A<-B->C :: fan out", 
  "A->B<-C :: fan in", 
  "A->B->C :: cascade",
  "A<->B<-C :: mutual in",
  "A<->B->C :: mutual out",
  "A->B<-C,A->C :: FFL",
  "A<-B<-C, A->C :: FBL",
  "A<->B<->C :: bi-mutual",
  "A<-B->C, A<->C :: regulated mutual",
  "A->B<-C, A<->C :: regulating mutual",
  "A->B->C, A<->C :: mutual cascade",
  "A->B<->C, A<->C :: semi clique",
  "A<->B<->C, A<->C :: clique")

triad_TotR <- triad_census(graphTot)
Tot_FFL <- triad_TotR[9]

sprintf("Number of total FFL: %s", Tot_FFL)


## CFFL1- Coherent type 1

#This time we will only take into account the positive regulatory effect subset graph


graphPR<-graph_from_adjacency_matrix(adj_PR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = TRUE)


triad_PR <- triad_census(graphPR)
CIFFL1 <- triad_PR[9]

sprintf("Number of coherent FFL type 1 : %s", CIFFL1)


## IFFL2 - Incoherent type 2

#This time we will only take into account the negative regulatory effect subset graph


graphNR<-graph_from_adjacency_matrix(adj_NR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = TRUE)


triad_NR <- triad_census(graphNR)
IIFFL2 <- triad_NR[9]

sprintf("Number of incoherent FFL type 2 : %s", IIFFL2)


#This time we will only take into account the unknown regulatory effect subset graph


graphXR<-graph_from_adjacency_matrix(adj_XR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = TRUE)

triad_XR <- triad_census(graphXR)
XIFFL1 <- triad_XR[9]

sprintf("Number of FFL with all-unknown regulation : %s", XIFFL1)


## Function to implement triad count in randomized replicas

# using 10 replicas for this random triad type frequency counts


FFL_freq<- function(adj_mat, typemat, replicas){
  
  df_tripclasses = data.frame(matrix(NA, nrow = 16, ncol = replicas))
  
  for (rep in 1:replicas){
    
    rnd_adj_mat <- randomize_network_qualitative(adj_mat, typemat)
    
    graph_rnd<-graph_from_adjacency_matrix(rnd_adj_mat, 
                                           mode="directed",
                                           weighted = NULL,
                                           diag = TRUE,
                                           add.colnames = FALSE,
                                           add.rownames = FALSE)
    
    
    count_tripclasses <- triad_census(graph_rnd)
    
    df_tripclasses[,rep]<- count_tripclasses
    
  }
  
  return(df_tripclasses)
}


#Get triad census for all types of regulations

#10 replicas of random adj matrix per reg. effect


n<-10
df_TotR<-FFL_freq(adj_TotR_convertedto1, "+", n)
df_PR<-FFL_freq(adj_PR, "+", n)
df_NR<-FFL_freq(adj_NR, "-", n)
df_XR<-FFL_freq(adj_XR, "?", n)

# means of random network triad type counts in a df

df_triadav = cbind(rowMeans(df_TotR), 
                   rowMeans(df_PR), 
                   rowMeans(df_NR),
                   rowMeans(df_XR))

# std of random network triad type counts in vector

std_FFL<- c(sd(df_TotR[9,]), 
            sd(df_PR[9,]), 
            sd(df_NR[9,]),
            sd(df_XR[9,]))

# annexing triad count data from random and wt networks into the same df

totaldf<-cbind(triad_TotR, df_triadav[, 1], 
			triad_PR, df_triadav[, 2], 
			triad_NR, df_triadav[, 3], 
			triad_XR, df_triadav[, 4])

colnames(totaldf)<- c("Tot_WT", "Tot_rnd", "PR_WT", "PR_rnd", "NR_WT", "NR_rnd", "XR_WT", "XR_rnd")
rownames(totaldf)<- triad_classes

df_triad_RND_WT<- as.data.frame(totaldf)


# writing frequencies in csv

write.csv(df_triad_RND_WT,"./triad_comp.csv")

df_triad_RND_WT


# calculation of p of observing the real wt FFL count in the random normal distributions of FFL frequency per regulatory effect
  

p_TotR  = pnorm(Tot_FFL, df_triadav[9,1], std_FFL[1], lower.tail = FALSE)
p_PR  = pnorm(CIFFL1, df_triadav[9,2], std_FFL[2], lower.tail = FALSE)
p_NR  = pnorm(IIFFL2, df_triadav[9,3], std_FFL[3], lower.tail = FALSE)


sign_rnd <- (rbind(cbind(Tot_FFL, df_triadav[9,1],std_FFL[1], p_TotR ), 
					cbind(CIFFL1, df_triadav[9,2],std_FFL[2], p_PR ), 
					cbind(IIFFL2, df_triadav[9,3],std_FFL[3], p_NR )))

colnames(sign_rnd)<- c("real_val", "rnd_mean", "rnd_std", "p(X>=o)")
rownames(sign_rnd)<- c("TotR", "PR", "NR")

sign_rnd<- as.data.frame(sign_rnd)

write.csv(sign_rnd, "./signif_FFL")

sign_rnd


# Boxplot of FFL frequency distributions per regulatory effect

font_add_google("Montserrat", "times")

## Automatically use showtext to render text
showtext_auto()
par(family="times")

wt_values<-c(Tot_FFL, CIFFL1, XIFFL1)

df_plot<- t(rbind(df_TotR[9,1:10], df_PR[9,1:10], df_NR[9,1:10]))
rnd_boxplot<-boxplot(df_plot, 
			names=c("TotR","PR", "NR"), 
			ylab="Number of FFL in 10 replicas", 
			xlab="Type of autorregulation", 
			col=c("#336699", "#FF9999", "#66CCCC"), 
			ylim = c(0,18))


text( 
  x=c(1:c(Tot_FFL, CIFFL1, IIFFL2)), 
  y=rnd_boxplot$stats[nrow(rnd_boxplot$stats),] + 1, 
  paste("* wt = ",c(Tot_FFL, CIFFL1, IIFFL2),sep="")  
)


# 2. AUTOREGULATION FREQUENCY

#------------------------------------------------------------------------
  
#  Testing with the randomized network

#Frequency generation per `adj_matrix`:** Sum through the `adj_matrix` diagonal.

#Random Replica frequency generation** with `auto_freq_generation` function


auto_freq_generation<- function(adj_mat, typemat, replicas){
  
  rnd_freq_vector<-vector()
  
  for (rep in 1:replicas){
    
    rnd_freq <- sum(diag(randomize_network_qualitative(adj_mat, typemat)))
    
    rnd_freq_vector<-c(rnd_freq_vector, rnd_freq)
  }
  
  return(rnd_freq_vector)
}


## Total autorregulation

#WT
#---------------------------------------------------------------------

#Taking into account unique autorregulatory motifs (per tf/gene unit, regardless of type of regulation)


real_TotR_unique <- sum(diag(adj_TotR_convertedto1))
print(c("WT unique Autorregulatory motif frequency: ",real_TotR_unique))


#Random for autorregulatory motifs taking into unique regulatory motifs regardless of types of regulation


n<-100

rnd_adj_TotR_unique <- auto_freq_generation(adj_TotR_convertedto1, "TOT", n)


# calculate acc random autorregulation frequency distribution function and
# check if real value fits into distribution taking into account mean and distribution - one tailed test 

mu = mean(rnd_adj_TotR_unique)
sigma = sd(rnd_adj_TotR_unique)
p  = pnorm(real_TotR_unique, mu, sigma, lower.tail = FALSE)

sprintf("The probability of random autorregulatory motif frequency taking real value %s or higher is %i",real_TotR_unique,p)
sprintf("Mean is %s", mu)
sprintf("Sd is %s", sigma)


#Taking into account that different types of regulation in autorregulatory 
#motifs would conform different autorregulatory motifs by themselves


real_TotR <- sum(diag(adj_TotR))
print(c("WT Autorregulatory motif with regulation type frequency: ",real_TotR))


#Random for autorregulatory TOTAL motifs taking into account types of regulation


n<-100

rndARmotif_TotR <- auto_freq_generation(adj_TotR, "TOT", n)


# calculate accumulated distribution function and check if real value fits 
#into distribution taking into account mean and distribution - one tailed test 

mu = mean(rndARmotif_TotR)
sigma = sd(rndARmotif_TotR)
p  = pnorm(real_TotR, mu, sigma, lower.tail = FALSE)

sprintf("The probability of random autorregulatory motif frequency taking real value %s or higher is %i",real_TotR,p)
sprintf("Mean is %s", mu)
sprintf("Sd is %s", sigma)


## Positive autorregulation


real_PR <- sum(diag(adj_PR))
print(c("WT Positive Autorregulatory motif frequency: ",real_PR))


#Random PR Replica


n<-100

rnd_PR <- auto_freq_generation(adj_PR, "+", n)


# calculate accumulated distribution function and check if real value fits 
#into distribution taking into account mean and distribution - one tailed test 

mu = mean(rnd_PR)
sigma = sd(rnd_PR)
p  = pnorm(real_PR, mu, sigma, lower.tail = FALSE)

sprintf("The probability of random autorregulatory motif frequency taking real value %s or higher is %s",real_PR,p)
sprintf("Mean is %s", mu)
sprintf("Sd is %s", sigma)


## Negative autorregulation


real_NR <- sum(diag(adj_NR))
print(c("WT Negative Autorregulatory motif frequency: ", real_NR))


# Random NR replica


n<-100

rnd_NR <- auto_freq_generation(adj_NR, "-", n)



# calculate accumulated distribution function and check if real value fits 
#into distribution taking into account mean and distribution - one tailed test 

mu = mean(rnd_NR)
sigma = sd(rnd_NR)
p  = pnorm(real_NR, mu, sigma, lower.tail = FALSE)

sprintf("The probability of random autorregulatory motif frequency taking real value %s or higher is %i",real_NR,p)
sprintf("Mean is %s", mu)
sprintf("Sd is %s", sigma)


## Unknown autorregulation XR


real_XR <- sum(diag(adj_XR))
print(c("WT Unknown Autorregulatory motif frequency: ",real_XR))


#Random XR Replica


n<-100


rnd_XR <- auto_freq_generation(adj_XR, "?", n)

# calculate accumulated distribution function and check if real value fits
# into distribution taking into account mean and distribution - one tailed test 

mu = mean(rnd_XR)
sigma = sd(rnd_XR)
p  = pnorm(real_XR, mu, sigma, lower.tail = FALSE)

sprintf("The probability of random autorregulatory motif frequency taking real value %s or higher is %s",real_XR,p)
sprintf("Mean is %s", mu)
sprintf("Sd is %s", sigma)


### Boxplot


font_add_google("Montserrat", "times")

## Automatically use showtext to render text
showtext_auto()
par(family="times")

rnd_boxplot<-boxplot(cbind(rnd_adj_TotR_unique, rnd_PR, rnd_NR), 
			names=c("TotR","PR", "NR"), 
			ylab="Number of unique autorregulatory motifs", 
			lab="Type of autorregulation", 
			col=c("#336699", "#FF9999", "#66CCCC"), 
			ylim=c(0,7))


wt_values<-c(real_TotR_unique, real_PR, real_NR)
text( 
  x=c(1:wt_values), 
  y=rnd_boxplot$stats[nrow(rnd_boxplot$stats),] + 0.5, 
  paste("* wt = ",wt_values,sep="")  
)

#-------------------------------------------------------------------------------

# 3. LONGEST / SHORTEST CASCADE

# we will use the igraph functions diameter and distance over the directed graphs built from the regulatory effect adjacency matrices

## Random longest and shortest cascade generation function


analysis_cascade_rndnet<- function(adj_mat, typemat, replicas){
  
  rnd_shortC_vector<-vector() # shortest cascade
  rnd_longC_vector<-vector() # longest cascade
  
  
  for (rep in 1:replicas){
    
    rnd_adj_mat <- randomize_network_qualitative(adj_mat, typemat)
    
    graph_rnd<-graph_from_adjacency_matrix(rnd_adj_mat, 
                                           mode="directed",
                                           weighted = NULL,
                                           diag = TRUE,
                                           add.colnames = FALSE,
                                           add.rownames = FALSE)
    
    
    rnd_min_dist<-distances(graph_rnd)
    rnd_min_dist[rnd_min_dist==0]<-max(rnd_min_dist)
    rnd_min_dist[rnd_min_dist==1]<-max(rnd_min_dist)
    rnd_min_dist<-min(rnd_min_dist)
    rnd_shortC_vector <- c(rnd_shortC_vector, rnd_min_dist)
    
    rnd_max_dist<-diameter(graph_rnd, directed = TRUE)
    rnd_longC_vector <- c(rnd_longC_vector, rnd_max_dist)
    
    
  }
  
  return(cbind(rnd_shortC_vector, rnd_longC_vector))
  
}


## Total regulation

#Real motifs in regulatory network


graphTot<-graph_from_adjacency_matrix(adj_TotR_convertedto1, 
                                      mode="directed",
                                      weighted = NULL,
                                      diag = TRUE,
                                      add.colnames = TRUE,
                                      add.rownames = NULL)

max_dist<-diameter(graphTot, directed = TRUE)

mat_dist<-distances(graphTot)
mat_dist[mat_dist==Inf]=0
mat_dist[mat_dist==0]=max(mat_dist)
mat_dist[mat_dist==1]=max(mat_dist)
sprintf("Shortest paths lengths from TOTAL adjacency matrix cascade : %s",  min(mat_dist))
sprintf("Longest paths lengths from TOTAL adjacency matrix cascade : %s", max_dist)

maxdist_TotR<-max_dist
mindist_TotR<-min(mat_dist)


# Random motifs in TOTAL regulatory networks length values X, and p(X>=o(l)) calculation

n<-10
rnd_TotR_lsC <- analysis_cascade_rndnet(adj_TotR_convertedto1, "+", n)

rnd_TotR_sC<-as.vector(rnd_TotR_lsC[,1])
rnd_TotR_lC<-as.vector(rnd_TotR_lsC[,2])

mu1 = mean(rnd_TotR_sC)
sigma1 = sd(rnd_TotR_sC)
p1  = pnorm(min(mat_dist), mu1, sigma1, lower.tail = FALSE)

mu2 = mean(rnd_TotR_lC)
sigma2 = sd(rnd_TotR_lC)
p2  = pnorm(max_dist, mu2, sigma2, lower.tail = FALSE)

sprintf("The probability of shortest random autorregulatory cascade length taking real value %s or higher is %s",min(mat_dist),p1)
sprintf("Mean shortest cascade : %s",mu1)
sprintf("std shortest cascade : %s",sigma1)

sprintf("The probability of shortest random autorregulatory cascade length taking real value %s or higher is %s",max_dist,p2)
sprintf("Mean longest cascade : %s",mu2)
sprintf("std longest cascade : %s",sigma2)


## Positive regulation

#Real motifs in regulatory network


graphPR<-graph_from_adjacency_matrix(adj_PR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = TRUE)

max_dist<-diameter(graphPR, directed = TRUE)
mat_dist<-distances(graphPR)
mat_dist[mat_dist==Inf]=0
mat_dist[mat_dist==0]=max(mat_dist)
mat_dist[mat_dist==1]=max(mat_dist)
sprintf("Shortest paths from Positive regulation adjacency matrix cascade : %s",  min(mat_dist))
sprintf("Longest paths from Positive regulation adjacency matrix cascade : %s", max_dist)

maxdist_PR<-max_dist
mindist_PR<-min(mat_dist)


#Random motifs in PR regulatory networks and p calculation


rnd_PR_lsC <- analysis_cascade_rndnet(adj_PR, "+", n)
rnd_PR_sC<-as.vector(rnd_PR_lsC[,1])
rnd_PR_lC<-as.vector(rnd_PR_lsC[,2])

mu1 = mean(rnd_PR_sC)
sigma1 = sd(rnd_PR_sC)
p1  = pnorm(min(mat_dist), mu1, sigma1, lower.tail = FALSE)

mu2 = mean(rnd_PR_lC)
sigma2 = sd(rnd_PR_lC)
p2  = pnorm(max_dist, mu2, sigma2, lower.tail = FALSE)

sprintf("The probability of shortest random positive cascade length taking real value %s or higher is %s",min(mat_dist),p1)
sprintf("Mean shortest cascade : %s",mu1)
sprintf("std shortest cascade : %s",sigma1)

sprintf("The probability of shortest random positive cascade length taking real value %s or higher is %s",max_dist,p2)
sprintf("Mean longest cascade : %s",mu2)
sprintf("std longest cascade : %s",sigma2)


## Negative regulation

#Real motifs in regulatory network


graphNR<-graph_from_adjacency_matrix(adj_NR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = TRUE)

max_dist<-diameter(graphNR, directed = TRUE)
mat_dist<-distances(graphNR)
mat_dist[mat_dist==Inf]=0
mat_dist[mat_dist==0]=max(mat_dist)
mat_dist[mat_dist==1]=max(mat_dist)
sprintf("Shortest paths from Positive regulation adjacency matrix cascade : %s",  min(mat_dist))
sprintf("Longest paths from Positive regulation adjacency matrix cascade : %s", max_dist)


maxdist_NR<-max_dist
mindist_NR<-min(mat_dist)


#Random motifs in NR regulatory networks and p calculation


n<-10
rnd_NR_lsC <- analysis_cascade_rndnet(adj_NR, "-", n)
rnd_NR_sC<-as.vector(rnd_NR_lsC[,1])
rnd_NR_lC<-as.vector(rnd_NR_lsC[,2])

mu1 = mean(rnd_NR_sC)
sigma1 = sd(rnd_NR_sC)
p1  = pnorm(min(mat_dist), mu1, sigma1, lower.tail = FALSE)

mu2 = mean(rnd_NR_lC)
sigma2 = sd(rnd_NR_lC)
p2  = pnorm(max_dist, mu2, sigma2, lower.tail = FALSE)

sprintf("The probability of shortest random negative cascade length taking real value %s or higher is %s",min(mat_dist),p1)
sprintf("Mean shortest cascade : %s",mu1)
sprintf("std shortest cascade : %s",sigma1)

sprintf("The probability of shortest random negative cascade length taking real value %s or higher is %s",max_dist,p2)
sprintf("Mean longest cascade : %s",mu2)
sprintf("std longest cascade : %s",sigma2)


## Unknown regulation

# Measure will not be representative as only 20 instances have unknown regulation


graphXR<-graph_from_adjacency_matrix(adj_XR, 
                                     mode="directed",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = TRUE,
                                     add.rownames = NULL)

max_dist<-diameter(graphXR, directed = TRUE)

mat_dist<-distances(graphXR)
mat_dist[mat_dist==Inf]=0
mat_dist[mat_dist==0]=max(mat_dist)

sprintf("Shortest paths from TOTAL adjacency matrix cascade : %s",  min(mat_dist))
sprintf("Longest paths from TOTAL adjacency matrix cascade : %s", max_dist)

maxdist_XR<-max_dist
mindist_XR<-min(mat_dist)

# RANDOM UNKNOWN SC LC LENGTH CALCULATION AND P CALCULATION

n<-10
rnd_XR_lsC <- analysis_cascade_rndnet(adj_XR, "?", n)
rnd_XR_sC<-as.vector(rnd_XR_lsC[,1])
rnd_XR_lC<-as.vector(rnd_XR_lsC[,2])

mu1 = mean(rnd_XR_sC)
sigma1 = sd(rnd_XR_sC)
p1  = pnorm(min(mat_dist), mu1, sigma1, lower.tail = FALSE)

mu2 = mean(rnd_XR_lC)
sigma2 = sd(rnd_XR_lC)
p2  = pnorm(max_dist, mu2, sigma2, lower.tail = FALSE)

sprintf("The probability of shortest random negative cascade length taking real value %s or higher is %s",min(mat_dist),p1)
sprintf("Mean shortest cascade : %s",mu1)
sprintf("std shortest cascade : %s",sigma1)

sprintf("The probability of shortest random negative cascade length taking real value %s or higher is %s",max_dist,p2)
sprintf("Mean longest cascade : %s",mu2)
sprintf("std longest cascade : %s",sigma2)


# Boxplot

##longest cascade


font_add_google("Montserrat", "times")

## Automatically use showtext to render text
showtext_auto()
par(family="times")

rnd_boxplot<-boxplot(cbind(rnd_TotR_lC, rnd_PR_lC, rnd_NR_lC, rnd_XR_lC), 
				names=c("TotR","PR", "NR", "XR"), 
				ylab="Length of longest cascade", 
				xlab="Type of autorregulation", 
				col=c("#336699", "#FF9999", "#66CCCC"), 
				ylim=c(0,80))

wt_values<-c(maxdist_TotR, maxdist_PR, maxdist_NR, maxdist_XR)
text( 
  x=c(1:wt_values), 
  y=rnd_boxplot$stats[nrow(rnd_boxplot$stats),] + 7, 
  paste("* wt = ",wt_values,sep="")  
)


  