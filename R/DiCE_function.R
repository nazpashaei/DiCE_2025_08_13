#' Differential Centrality-Ensemble Analysis Based on Gene Expression Profiles and Protein-Protein Interaction Network
#'
#' This function offers a comprehensive framework for identifying and prioritizing disease-related genes.
#'
#' @param data A list of 2 elements:
#'   \itemize{
#'     \item \code{data} - A data frame of gene expression data with gene symbols as column names and the last column as class labels (Tumor, Normal).
#'     \item \code{DE} - A data frame of DEGs with three columns: \code{Gene.symbol}, \code{adj.P.Val}, and \code{logFC}.
#'   }
#' @param regulation_status A string specifying the direction of regulation to filter genes. Acceptable values: "up", "down", or "both".
#' @param species A string specifying the species. Acceptable values: "human" (9606), "mouse" (10090), or "rat" (10116).
#' @return A data frame of key genes with their gene names, log fold change (logFC), and ensemble ranking.
#' @examples
#' \dontrun{
#'   # Example data
#'   data <- list(
#'     data = data.frame(gene1 = c(1, 2), gene2 = c(3, 4), class = c("Tumor", "Normal")),
#'     DE = data.frame(Gene.symbol = c("gene1", "gene2"), adj.P.Val = c(0.01, 0.02), logFC = c(2, -2))
#'   )
#'   KeyGenes <- DiCE_function(data)
#' }
#' @export


DiCE_function <- function(data,regulation_status,species,method,pval_threshold, log2fc_threshold,pval_type,case_label, control_label){

  library(dplyr)
  library(tibble)
  library(FSelectorRcpp)
  library(igraph)
  library(data.table)
  library(NetWeaver)
  library(praznik)
  # The data must be a list of 2 elements. The first element is a data frame of gene expression data.
  # Column names are gene symbols; the last column is the class label (Tumor, Normal).
  # The second list element is a list of DEGs named DE and is a data frame with three columns: Gene.symbol, adj.P.Val, and logFC.  data = data
  #===================================
  str(data$DE)
  process_gene_names <- function(Gene.symbol) {
    # Count the occurrences of '///' in the gene name
    slashes_count <- length(gregexpr("///", Gene.symbol)[[1]])
    if (slashes_count == 1) {
      # If only one '///', consider the second word after '///'
      Gene.symbol <- gsub(".*///(\\w+).*", "\\1", Gene.symbol)
    } else if (slashes_count > 1) {
      # If more than one '///', consider the second word between '////'
      Gene.symbol <- gsub(".*///(.+?)///.*", "\\1", Gene.symbol)
    }

    return(Gene.symbol)
  }

  # Apply the function to each gene name
  #pseudogenes and lncRNA's
  data$DE$Gene.symbol <- sapply(data$DE$Gene.symbol, process_gene_names)
  colnames(data$data)<-sapply(colnames(data$data), process_gene_names)
  data$DE <- subset(data$DE, !grepl("LOC", data$DE$Gene.symbol))
  data$DE <- subset(data$DE, !grepl("LINC", data$DE$Gene.symbol))
  #data$DE$Gene.symbol<-toupper(data$DE$Gene.symbol)
  data$DE$Gene.symbol<-(data$DE$Gene.symbol)


  #Phase I: Construction of a candidate gene pool by DEA with a loose cutoff
  #=======================================
  if (!(pval_type %in% c("adj.P.Val", "P.Value"))) {
  stop("Invalid pval_type. Use 'adj.P.Val' or 'P.Value'.")
}
pvals <- data$DE[[pval_type]]

if (regulation_status == "Up") {
  dee <- data$DE[(pvals <= pval_threshold) & (data$DE$logFC >= log2fc_threshold), ]
} else if (regulation_status == "Down") {
  dee <- data$DE[(pvals <= pval_threshold) & (data$DE$logFC <= -log2fc_threshold), ]
} else {
  dee <- data$DE[(pvals <= pval_threshold) & (abs(data$DE$logFC) >= log2fc_threshold), ]
}
   
  dee1=dee
  colnames(dee1)=c("gene_name","P.Value","adj.P.Val","logFC")

  #+++++++++++++++++++++++Phase II: Selection of the top discriminative genes from the candidate pool obtained in Phase I using the Information Gain (IG) filter approach.

  class=data$data[,ncol(data$data)]
  if (!all(c(case_label, control_label) %in% class)) {
  stop("Class labels provided not found in data")
}
  #colnames(data$data)<-toupper(colnames(data$data))
  data$data=data$data[,colnames(data$data)%in%dee1$gene_name]

  rownames(data$data)<-NULL
  d_mat=data$data
  #d_mat1=scale(d_mat);
  d_mat1=d_mat;
  d_mat2=as.data.frame(d_mat1);

  newMydata=cbind(d_mat2,class)

  weights <- information_gain(class~., newMydata,type = "infogain");#View(weights);#c("infogain", "gainratio", "symuncert")
  s=as.data.frame(weights$importance)
  row.names(s)=weights$attributes
  #----------------------------------------Mutual information scores
  #weights <- miScores(d_mat2,class);weights=data.frame(weights);s=weights;
  #row.names(s)=row.names(weights)
  
  vi <- s
  vi$max <- apply(vi, 1, max)
  VS=vi[order(-vi$max),]
  rn=rownames(VS);ddds=cbind(rn,VS,1:nrow(VS));ddds=ddds[,-3];
  colnames(ddds)=c("gene_name","weight","rank");#View(ddds)
  ddds1<- ddds[-which(ddds$weight==0),];#View(ddds1);
  oo=ddds[,-1];oo1=oo[,c(2,1)]; plot(oo1)

  x2<- ddds[-which(ddds$weight<mean(ddds$weight,na.rm = TRUE)),];#View(x2);dim(x2);
  m1=merge(x2,dee1,by="gene_name");m2=m1[order(m1$rank),];#View(m2)
  setdiff(x2$gene_name,m2$gene_name)
  #++++++++++++++++++++++PPInetwork
  get_species_id <- function(species) {
     species <- tolower(species)
     if (species %in% c("human", "homo sapiens")) {
        return(9606)
    } else if (species %in% c("mouse", "mus musculus")) {
       return(10090)
    } else if (species %in% c("rat", "rattus norvegicus")) {
       return(10116)
    } else {
       stop("Unsupported species. Please use 'human', 'mouse', or 'rat'.")
    }
   }
  
  library(STRINGdb)
  string_db <- STRINGdb$new(version="12", species=get_species_id(species),score_threshold=400, input_directory="");#protocol="http"
  p=m2[,-c(2,3)];p=as.data.frame(p)
  p_mapped <- string_db$map(p, "gene_name", takeFirst=TRUE, removeUnmappedRows=TRUE, quiet=FALSE)
  neighbors <- string_db$get_neighbors(p_mapped$STRING_id);str(neighbors)
  inter<-string_db$get_interactions(p_mapped$STRING_id);str(inter)

  #+++++++++++++++++++++Expression matrix based on STRING name
   table(duplicated(p_mapped$gene_name));print(p_mapped[duplicated(p_mapped$gene_name),])
  dup_genes <- p_mapped$gene_name[duplicated(p_mapped$gene_name) | duplicated(p_mapped$gene_name, fromLast = TRUE)]
  p_mapped_unique <- p_mapped[!(p_mapped$gene_name %in% dup_genes), ]
  nn<-inter;whol=nn[,-3];colnames(whol)=c("node1","node2");bl=whol
  vertex=c(whol$node1,whol$node2);vertex=unique(vertex);length(vertex)
  p_mapped_filtered <- p_mapped[p_mapped$gene_name %in% dup_genes & p_mapped$STRING_id %in% vertex, ]
  p_mapped=rbind(p_mapped_unique,p_mapped_filtered)
  p_mapped1 <- p_mapped[!duplicated(p_mapped$gene_name), ]  # remove duplicates

  table(is.na(p_mapped1$STRING_id))
  convert_ORF <- function(x) {
      parts <- strsplit(x, "ORF")[[1]]
      paste0(parts[1], "orf", parts[2])
  }
   p_mapped1$gene_name <- ifelse(substr(p_mapped1$gene_name, 1, 1) == "C" & (grepl("ORF", substr(p_mapped1$gene_name, 3, 5))|grepl("ORF", substr(p_mapped1$gene_name, 4, 6))),
                             sapply(p_mapped1$gene_name, convert_ORF),
                            p_mapped1$gene_name)
  
  setdiff(colnames(d_mat),p_mapped1$gene_name)
  s=d_mat[,colnames(d_mat)%in%p_mapped1$gene_name];dim(s)
  s1=as.data.frame(s)
  setdiff(colnames(s1),p_mapped1$gene_name)

  setnames(s1,old=p_mapped1$gene_name, new=p_mapped1$STRING_id,skip_absent=TRUE);#skip_absent=TRUE)

  expression1=cbind(s1,class);dim(expression1)
  table(duplicated(colnames(expression1)))

  #++++++++++++++++++++Phase III: Modification of PPIs network by calculating Pearson correlation coefficients (c.c.) of paired genes in terms of gene expression for a specific phenotype or under one condition and assigning (1−|𝑐.𝑐|) as the normalized distance between paired genes by the PPI network.

  nn<-inter;whol=nn[,-3];colnames(whol)=c("node1","node2");bl=whol
  vertex=c(whol$node1,whol$node2);vertex=unique(vertex);length(vertex)

  test1=expression1[expression1$class=="control_label",];
  test1=test1[,-ncol(test1)];
  name=colnames(test1);
  df=test1; 

  #mat1 <- abs(cor(df,method = "pearson"))
  mat1 <- abs(cor(df,method = "spearman"))
  com=intersect(vertex,colnames(test1));length(com);length(vertex);setdiff(vertex,com)
  #------------------------
  dif1=setdiff(vertex,com);
  vertex=vertex[!vertex %in% c(dif1)]
  #--------------------
  BB=mat1[vertex,vertex]
  bb=as.data.frame(BB)
  aa=stack(bb);
  aa$ind=as.character(aa$ind)
  rn=row.names(BB);lr=ncol(BB)
  ind2<-rep(rn, times =lr)
  aa=cbind(aa,ind2)
  aa=aa[order(match(aa[,3],bl[,1])),]
  concat=paste(bl[,2],bl[,1])

  sss=aa[paste(aa$ind,aa$ind2)%in%concat,]
  sss=cbind(sss[3],sss[2],sss[1])
  colnames(sss)=c("source", "target", "weight")
  g <- graph_from_data_frame(sss,directed = FALSE)
  is_weighted(g)
  g <- set_edge_attr(g, "weight", value=(1-(sss$weight))+0.001)
  E(g)$weight[is.na(E(g)$weight)] <- 1

  V(g); E(g)
  gsize(g)
  gorder(g)
  is_weighted(g)
  #Phase IV: Topological analysis of sample-specific weighted PPI networks from the perspectives of Betweenness and Eigenvector centrality, followed by differential analysis on each centrality measure of individual genes.
  Strength=sort(strength(g),decreasing = TRUE);rank=seq(1,length(Strength),1);Strength=cbind(Strength,rank);#View(Strength)
  centralities <- cbind(Eig=evcent(g)$vector,Betweenness=betweenness(g));#View(centralities)

  centralities.Normal <-centralities
  #***********************************
  nn<-inter;whol=nn[,-3];colnames(whol)=c("node1","node2");bl=whol
  vertex=c(whol$node1,whol$node2);vertex=unique(vertex);length(vertex)

  test1=expression1[expression1$class=="case_label",];
  test1=test1[,-ncol(test1)];
  name=colnames(test1);
  df=test1
  #mat1 <- abs(cor(df,method = "pearson"))
  mat1 <- abs(cor(df,method = "spearman"))
  com=intersect(vertex,colnames(test1));length(com);length(vertex);setdiff(vertex,com)
  #------------------------
  dif1=setdiff(vertex,com);
  vertex=vertex[!vertex %in% c(dif1)]
  #--------------------
  BB=mat1[vertex,vertex]
  bb=as.data.frame(BB)
  aa=stack(bb);
  aa$ind=as.character(aa$ind)
  rn=row.names(BB);lr=ncol(BB)
  ind2<-rep(rn, times =lr)
  aa=cbind(aa,ind2)
  aa=aa[order(match(aa[,3],bl[,1])),]
  concat=paste(bl[,2],bl[,1])

  sss=aa[paste(aa$ind,aa$ind2)%in%concat,]
  sss=cbind(sss[3],sss[2],sss[1])
  colnames(sss)=c("source", "target", "weight")
  g <- graph_from_data_frame(sss,directed = FALSE)
  is_weighted(g)
  g <- set_edge_attr(g, "weight", value=(1-(sss$weight))+0.001)
  E(g)$weight[is.na(E(g)$weight)] <- 1
  V(g); E(g)
  gsize(g)
  gorder(g)
  is_weighted(g)
  Strength=sort(strength(g),decreasing = TRUE);rank=seq(1,length(Strength),1);Strength=cbind(Strength,rank);#View(Strength)
  centralities <- cbind(Eig=evcent(g)$vector,Betweenness=betweenness(g));#View(centralities)

  centralities.Tumor <-centralities

  #+++++++++++++++++++++++++++++++++++++++++Ensemble ranking
  #Phase V: Ensemble ranking to produce a robust and trustworthy level of rating based on the ranks established in the previous stage.

  df.tmp <- cbind(centralities.Tumor,centralities.Normal);#View(df.tmp)
  #-----------------------------------------------
  df.eig <- df.tmp[,c(1,3)]
  df.betweenness <- df.tmp[,c(2,4)]
  rank=seq(1,nrow(df.eig),1)

  rr=0;result=0
  for(i in 1:nrow(df.eig)){
    result[i] <-df.eig[i,1]-df.eig[i,2]
  }
  df.eig=cbind(df.eig,result)

  df.eig <- df.eig[order(abs(df.eig[,3]),decreasing = TRUE),];df.eig=(cbind(df.eig,rank))

  rr=0;result=0
  for(i in 1:nrow(df.betweenness)){
    result[i] <-df.betweenness[i,1]-df.betweenness[i,2]
  }

  df.betweenness=cbind(df.betweenness,result)
  df.betweenness <- df.betweenness[order(abs(df.betweenness[,3]),decreasing = TRUE),];
  df.betweenness=(cbind(df.betweenness,rank))

  #-------------------
  common_row_names <- intersect(rownames(df.betweenness), rownames(df.eig))
  df.betweenness_common <- df.betweenness[common_row_names, ]
  df.eig_common <- df.eig[common_row_names, ]
  df.tmp22 <- cbind(df.betweenness_common, df.eig_common)
  str(df.tmp22)
  df.tmp22[is.na(df.tmp22)] <- 0

  enr=cbind(df.tmp22[,4],df.tmp22[,8])
  ensembleRanking=ensemble_rank(enr,method='ProductOfRank')
  df.tmp22=cbind(df.tmp22,ensembleRanking);
  #sort based on ranking
  df.tmp22=df.tmp22[order(df.tmp22[,9],decreasing = TRUE ),]
  df.tmp22=(cbind(df.tmp22,rank));df.tmp22=as.data.frame(df.tmp22)
  pg=p_mapped1
  #gene. prioritized
  merged_df <- merge(pg, df.tmp22, by.x = "STRING_id", by.y = "row.names", all.y = TRUE)
  df2 <- merged_df[order(merged_df$rank.2,decreasing=FALSE),];#View(df2)
  #write.csv(df2, file = "gene.prioritized.csv")


  #++++++++++++++++++++++++++++identify-KeyGenes
  #df2=df2[-which((df2$Betweenness==0)&(df2$Betweenness.1==0)),]

  x<-df2
  # Compute thresholds based on the chosen method
  if (method == "mean") {
    B <- mean(df2$Betweenness)
    B1 <- mean(df2$Betweenness.1)
    E <- mean(df2$Eig)
    E1 <- mean(df2$Eig.1)
  } else if (method == "top25") {
    B <- quantile(df2$Betweenness, 0.75)
    B1 <- quantile(df2$Betweenness.1, 0.75)
    E <- quantile(df2$Eig, 0.75)
    E1 <- quantile(df2$Eig.1, 0.75)
  } else if (method == "median") {
    B <- median(df2$Betweenness)
    B1 <- median(df2$Betweenness.1)
    E <- median(df2$Eig)
    E1 <- median(df2$Eig.1)
  } else {
    stop("Method must be one of: 'mean', 'top25', or 'median'")
  }
   
  x1<- x[-which((df2$Betweenness<B)&(df2$Betweenness.1<B1)),];dim(x1);x2<- x[-which((df2$Eig<E)&(df2$Eig.1<E1)),];dim(x2)

  common1=intersect(x1$gene_name,x2$gene_name)
  common1

  lc=df2[df2$gene_name%in%common1,];dim(lc)
  DiCE.genes=0; DiCE.genes=lc[,c(2,4,14)];colnames(DiCE.genes)<-c("gene_name", "logFC", "ensemble.Ranking");
  DiCE.genes$Final.rank <- 1:nrow(DiCE.genes)
cat("Phase I =", nrow(dee1), "\n",
    "Phase II =", nrow(m2), "\n",
    "Phase III =", nrow(df2), "\n",
    "DiCe-genes =", length(common1), "\n")
# Open a connection for writing
file_conn <- file("DiCE.output.csv", open = "wt")

# Write the summary information first
writeLines(c(
    paste("Phase I =", nrow(dee1)),
    paste("Phase II =", nrow(m2)),
    paste("Phase III =", nrow(df2)),
    paste("DiCE-genes =", length(common1)),
    ""  # Blank line before the table
), file_conn)

# Merge df2 columns (result, result.1, rank.2) into dee1 by gene_name
dee1_merged <- merge(dee1, df2[, c("gene_name", "result", "result.1", "rank.2")], 
                     by = "gene_name", all.x = TRUE)



# Create a new column indicating presence in each phase
dee1_merged$Phase <- apply(dee1_merged, 1, function(row) {
  phases <- c()
  if (row["gene_name"] %in% dee1$gene_name) phases <- c(phases, "I")
  if (row["gene_name"] %in% m2$gene_name) phases <- c(phases, "II")
  if (row["gene_name"] %in% df2$gene_name) phases <- c(phases, "III")
  if (row["gene_name"] %in% DiCE.genes$gene_name) phases <- c(phases, "DiCE")
  
  if (length(phases) == 0) return("-")
  return(tail(phases, 1))  # return only the last (most recent) matched phase
})

# Merge DiCE.genes column (Final.rank) into the merged data
dee1_merged <- merge(dee1_merged, DiCE.genes[, c("gene_name", "Final.rank")], 
                     by = "gene_name", all.x = TRUE)

dee1_merged[is.na(dee1_merged)] <- "-"


colnames(dee1_merged)<-c("Offical gene symbol from input file", "p-value","adj.P.Val","log2FC","DB","DE","Ensemble.rank","DiCE Phase","Final.rank");

# Append the data
write.table(dee1_merged, file = file_conn, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)

# Close the file connection
close(file_conn)
#------------------------------------------:) end 
  View(DiCE.genes)
  return(DiCE.genes)
}
#------------------------------------------:) end

