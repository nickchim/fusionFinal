require(readr)
require(stringr)
require(tibble)
#library(biomaRt)


source("getPFAMDomain.R")

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
input <- read.delim("star_fusion_combo.tsv")

#Cancer Gene and Fusion Prep
geneList <-read_tsv("lists/CancerGeneList.tsv")
fuseList <-read_tsv("lists/FusionList2.txt")

input$Cancerous_Gene <- F
input$Cancerous_Gene_Symbol <- NA
input$Cancerous_Fusion <- F
input$Cancerous_Fusion_Symbol <- NA
input$Fusion_Cancer_Type <- NA

geneRep <- 1
fuseRep <- 1
repAdjustedGene <- geneList[geneList$Count >= geneRep,]
repAdjustedFuse <- fuseList[fuseList$Count >= fuseRep,]

#Target Prep
input$Targetable <- FALSE
input$Targetable_gene <- NA
input$Drug_name <- NA
input$Drug_chembl_id <- NA

Target_genes <- read_tsv("lists/Target_List.tsv")

#Get PFAM Domain
input <- getPFAMDomain(input)


#Domain Prep
input$Cancer_Domains <- NA
domainList <-read_tsv("lists/Domain_List.tsv")

#Big Apply Function
#Annotate with Cancerous gene, Cancerous Fusion, Targetable info and Domain targetability
bigApply <- function(x) {
  #Cancer Gene and Fusion Checker
  gene1 <- x["H_Gene"]
  gene2 <- x["T_Gene"]
  fuse <- paste(gene1,gene2,sep="--")
  
  gene1In <- gene1 %in% repAdjustedGene$Gene
  gene2In <- gene2 %in% repAdjustedGene$Gene
  fuseIn <- fuse %in% repAdjustedFuse$Fusions
  
  x["Cancerous_Gene"] <- gene1In || gene2In
  x["Cancerous_Fusion"] <- fuseIn
  
  if (gene1In && gene2In) {
    x["Cancerous_Gene_Symbol"] <- paste(gene1,gene2,sep=",")
  }else if (gene1In) {
    x["Cancerous_Gene_Symbol"] <- gene1
  }else if (gene2In) {
    x["Cancerous_Gene_Symbol"] <- gene2
  }
  
  if (fuseIn) {
    x["Cancerous_Fusion_Symbol"] <- paste(gene1,gene2,sep="--")
    
    cancers <- unlist(str_split(fuseList[fuseList$Fusions == paste(gene1,gene2,sep="--"),"Cancer_Types"],","))
    indices <- which(fuseList[fuseList$Fusions == paste(gene1,gene2,sep="--"),-c(1,2,ncol(fuseList))] != 0) + 2
    for (i in 1:length(cancers)) {
      if (!is.na(cancers[i])) {
        if (i == 1) {
          x["Fusion_Cancer_Type"] <- paste(cancers[i],indices[i],sep=":")
        }else {
          x["Fusion_Cancer_Type"] <- paste(x["Fusion_Cancer_Type"],paste(cancers[i],indices[i],sep=":"),sep=",")
        }
      }
    }
  }
  
  #Targetable
  if (x["H_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["H_Gene"]
    
    x["Targetable_gene"] <- x["H_Gene"]
    x["Drug_name"] <- Target_genes$Drug_name[index]
    x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
  }
  
  if (x["T_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["T_Gene"]
    
    if (is.na(x["Targetable_gene"])) {
      x["Targetable_gene"] <- x["T_Gene"]
      x["Drug_name"] <- Target_genes$Drug_name[index]
      x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
    }else {
      x["Targetable_gene"] <- paste(x["Targetable_gene"],x["T_Gene"],sep=",")
      x["Drug_name"] <- paste(x["Drug_name"],Target_genes$Drug_name[index],sep=",")
      x["Drug_chembl_id"] <- paste(x["Drug_chembl_id"],Target_genes$Drug_chembl_id[index],sep=",")
    }
    
  }
  
  #Domain Checker
  doms <- c()
  domList <- str_split(unlist(str_split(x["H_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  doms <- c()
  domList <- str_split(unlist(str_split(x["T_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  if (!is.na(x["Cancer_Domains"])) {
    list <- str_split(x["Cancer_Domains"],",")[[1]]
    x["Cancer_Domains"] <- paste(levels(factor(list)),collapse=",")
  }
  return(x)
}

#Applying bigApply
input <- apply(input, MARGIN = 1, FUN = bigApply)

#Cleaning up table and adjusting order of columns
input <- t(input)
input <- as.tibble(input)
input <- input[,c("Fusion_Name","H_Gene","T_Gene","Left_Breakpoint_Chr","Left_Breakpoint_Pos","Left_Breakpoint_Str","Right_Breakpoint_Chr","Right_Breakpoint_Pos","Right_Breakpoint_Str","Star_Junction_Readcount","Star_Spanning_Fragcount","FC_Common_Mapping_Reads","FC_Spanning_Pairs","FC_Spanning_Unique_Reads","Fusion_Transcript_Sequence","Fusion_Protein_Sequence","Fusion_Type","Catcher","Cancerous_Gene","Cancerous_Gene_Symbol","Cancerous_Fusion","Cancerous_Fusion_Symbol","Fusion_Cancer_Type","Targetable","Targetable_gene","Drug_name","Drug_chembl_id","H_Gene_PFAM_All","H_Gene_PFAM_IN_FUSION","T_Gene_PFAM_All","T_Gene_PFAM_IN_FUSION","Cancer_Domains")]

input <- add_column(input,H_Gene_PFAM_NOT_IN_FUSION=NA,.before = "T_Gene_PFAM_All")
input <- add_column(input,T_Gene_PFAM_NOT_IN_FUSION=NA,.before = "Cancer_Domains")

#Annotate Domains with in fusion and not in fusion
domainFunc <- function(x) {
  allFuse <- unlist(str_split(unlist(str_split(x["H_Gene_PFAM_All"],":")),"; "))
  allFuse <- allFuse[seq(1,length(allFuse),by=2)]
  inFuse <- unlist(str_split(unlist(str_split(x["H_Gene_PFAM_IN_FUSION"],":")),"; "))
  inFuse <- inFuse[seq(1,length(inFuse),by=2)]

  notInFuse <- allFuse[which(!allFuse %in% inFuse)]
  if (length(notInFuse) > 0) {x["H_Gene_PFAM_NOT_IN_FUSION"] <- paste(notInFuse,collapse=",")}

  allFuse <- unlist(str_split(unlist(str_split(x["T_Gene_PFAM_All"],":")),"; "))
  allFuse <- allFuse[seq(1,length(allFuse),by=2)]
  inFuse <- unlist(str_split(unlist(str_split(x["T_Gene_PFAM_IN_FUSION"],":")),"; "))
  inFuse <- inFuse[seq(1,length(inFuse),by=2)]


  notInFuse <- allFuse[which(!allFuse %in% inFuse)]
  if (length(notInFuse) > 0) {x["T_Gene_PFAM_NOT_IN_FUSION"] <- paste(notInFuse,collapse=",")}

  return(x)
}
input <- apply(input,FUN=domainFunc,MARGIN = 1)
input <- as.data.frame(t(input))



#Tiering Prep
input$Tier <- 4

#Tierize func
tierize <- function(x,input) {
  cancerousGene <- as.logical(x["Cancerous_Gene"])
  cancerousFusion <- as.logical(x["Cancerous_Fusion"])
  spanning <- FALSE
  if (!is.na(x["FC_Spanning_Pairs"]) && !is.na(x["FC_Spanning_Unique_Reads"])) {
    spanning <- x["FC_Spanning_Pairs"] > 3 || x["FC_Spanning_Unique_Reads"] > 3
  }else if (!is.na(x["Star_Spanning_Fragcount"]) && !is.na(x["Star_Junction_Readcount"])) {
    spanning <- x["Star_Spanning_Fragcount"] > 3 || x["Star_Junction_Readcount"] > 3
  }

  star <- x["Catcher"] == "Star"
  fusion <- x["Catcher"] == "Fusion_Catcher"
  starAndFusion <- FALSE
  if (sum(x["Fusion_Name"] == input[,"Fusion_Name"]) > 1) {
    if (any(input[x["Fusion_Name"] == input[,"Fusion_Name"],"Catcher"] != x["Catcher"])) {
      starAndFusion <- TRUE
    }
  }
  targetable <- as.logical(x["Targetable"])


  if (cancerousGene & spanning & starAndFusion & targetable) {
    tier <- "1A"
  }else if (cancerousFusion & spanning & starAndFusion & targetable) {
    tier <- "1B"
  }else if (cancerousGene & spanning & targetable) {
    tier <- "2A"
  }else if (cancerousFusion & spanning & targetable) {
    tier <- "2B"
  }else if (cancerousGene & targetable) {
    tier <- "3A"
  }else if (cancerousFusion & targetable) {
    tier <- "3B"
  }else {
    tier <- "4"
  }

  x["Tier"] <- tier
  return(x)
}
input <- t(apply(input, FUN = tierize, MARGIN = 1,input))

write.table(input, "star_fusion_analyzed.tsv", sep="\t", row.names = FALSE)
