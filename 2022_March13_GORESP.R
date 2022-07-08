

####  GO enrichment function
####  curr_exp is name of experiment
####  mat = experimental matrix
####  coln = coln of mat with fitness defect scores for enrichment
####  sig = significance threshold for GO enrichment of fitness defect scores
####  bp_path = path of GOSET RDS file
####  minGeneSetSize = minimum size of geneSet to include in enrichment
####  maxSetSize = (not an option) maximum size of geneSet to include in enrichment -- set to 300


# RETURNS a dataframe of enrichment results, sorted by increasing FDR value. The columns are:
#         term = name of gene set
#         querySetFraction = the fraction of the query set that overlaps with the term set
#         geneSetFraction = the fraction of the term set that overlaps with the query set
#         foldEnrichment = the fold enrichment of the query set with the term genes
#         P = P value estimating the significance with which the query set is enriched with the term genes
#         FDR = FDR value estimating the significance of enrichment
#         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#	  maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
####
####  query fraction is length of overlap divided by number of genes scored significant in screen
####  geneSet is length of genes in the geneSet
####  overlap is length of intersect
####  bgRate round rate is nGenes /total number of rows in unigene set (matrix)
####  foldenrichment query fraction /bgRate
####  foldenrichment * bgRate = query fraction



runGORESP = function (fdrThresh = 0.2, curr_exp = colnames(mat)[coln], mat,coln,sig = 1,
                      bp_path = "2021_Decenber30_GO_BP.RDS",go_path = "2021_Decenber30_GOID_GOBP_SGD.txt",
                      bp_input = NULL, go_input = NULL,minGeneSetSize = 5,maxGeneSetSize=300){
#   NONSPECIFIC.TERMS <- list(mf=c("MOLECULAR_FUNCTION", "BINDING", "CATALYTIC ACTIVITY"),
#                           cc=c("CELL", "CELL CORTEX PART", "CELL DIVISION SITE PART", "CELL FRACTION", "CELL PART", "CELL PERIPHERY", "CELL PROJECTION PART", "CELL WALL PART", "CELLULAR_COMPONENT", "CHROMOSOMAL PART", "CYTOPLASMIC PART", "CYTOPLASMIC VESICLE PART", "CYTOSKELETAL PART", "CYTOSOLIC PART", "ENDOPLASMIC RETICULUM PART", "ENDOSOMAL PART", "EXTERNAL ENCAPSULATING STRUCTURE", "EXTERNAL ENCAPSULATING STRUCTURE PART", "EXTRINSIC TO MEMBRANE", "GOLGI APPARATUS PART", "INSOLUBLE FRACTION", "INTEGRAL TO MEMBRANE", "INTEGRAL TO MEMBRANE OF MEMBRANE FRACTION", "INTRACELLULAR", "INTRACELLULAR ORGANELLE", "INTRACELLULAR ORGANELLE LUMEN", "INTRACELLULAR ORGANELLE PART", "INTRACELLULAR PART", "INTRACELLULAR MEMBRANE-BOUNDED ORGANELLE", "INTRACELLULAR NON-MEMBRANE-BOUNDED ORGANELLE", "INTRINSIC TO MEMBRANE", "MEMBRANE", "MEMBRANE-BOUNDED ORGANELLE", "MEMBRANE-ENCLOSED LUMEN", "MEMBRANE FRACTION", "MEMBRANE PART", "MICROBODY PART", "MICROTUBULE ORGANIZING CENTER PART", "MITOCHONDRIAL MEMBRANE PART", "MITOCHONDRIAL PART", "NON-MEMBRANE-BOUNDED ORGANELLE", "NUCLEAR CHROMOSOME PART", "NUCLEAR MEMBRANE PART", "NUCLEAR PART", "NUCLEOLAR PART", "NUCLEOPLASM PART", "ORGANELLE", "ORGANELLE INNER MEMBRANE", "ORGANELLE LUMEN", "ORGANELLE MEMBRANE", "ORGANELLE OUTER MEMBRANE", "ORGANELLE PART", "ORGANELLE SUBCOMPARTMENT", "PERIPHERAL TO MEMBRANE OF MEMBRANE FRACTION", "PEROXISOMAL PART", "PLASMA MEMBRANE ENRICHED FRACTION", "PLASMA MEMBRANE PART", "VACUOLAR PART", "VESICULAR FRACTION"),
#                           bp=c("POSITIVE REGULATION OF MACROMOLECULE METABOLIC PROCESS", "REGULATION OF CELLULAR COMPONENT ORGANIZATION", "POSITIVE REGULATION OF METABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR METABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR PROCESS", "REGULATION OF CELLULAR PROCESS", "CELLULAR NITROGEN COMPOUND BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF NITROGEN COMPOUND METABOLIC PROCESS", "REGULATION OF CATALYTIC ACTIVITY", "POSITIVE REGULATION OF CATALYTIC ACTIVITY", "REGULATION OF MOLECULAR FUNCTION", "POSITIVE REGULATION OF CELLULAR COMPONENT ORGANIZATION", "REGULATION OF ORGANELLE ORGANIZATION", "POSITIVE REGULATION OF CATABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR CATABOLIC PROCESS", "POSITIVE REGULATION OF MOLECULAR FUNCTION", "REGULATION OF CATABOLIC PROCESS", "REGULATION OF CELLULAR CATABOLIC PROCESS", "CELLULAR RESPONSE TO CHEMICAL STIMULUS", "CELLULAR RESPONSE TO ORGANIC SUBSTANCE", "POSITIVE REGULATION OF BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF CELLULAR BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF MACROMOLECULE BIOSYNTHETIC PROCESS", "CELLULAR CARBOHYDRATE METABOLIC PROCESS", "REGULATION OF CELLULAR PROTEIN METABOLIC PROCESS", "REGULATION OF PROTEIN METABOLIC PROCESS", "NEGATIVE REGULATION OF BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF CELLULAR BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF CELLULAR MACROMOLECULE BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF MACROMOLECULE METABOLIC PROCESS", "RESPONSE TO EXTERNAL STIMULUS", "RESPONSE TO EXTRACELLULAR STIMULUS", "CELLULAR HOMEOSTASIS", "HOMEOSTATIC PROCESS", "REGULATION OF HOMEOSTATIC PROCESS", "ORGANIC SUBSTANCE TRANSPORT", "CELLULAR NITROGEN COMPOUND CATABOLIC PROCESS", "ORGANIC ACID BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF ORGANELLE ORGANIZATION", "ORGANELLE FISSION", "NEGATIVE REGULATION OF CELLULAR COMPONENT ORGANIZATION", "NEGATIVE REGULATION OF NITROGEN COMPOUND METABOLIC PROCESS", "CELLULAR DEVELOPMENTAL PROCESS", "MAINTENANCE OF LOCATION IN CELL","REGULATION OF DEVELOPMENTAL PROCESS","SMALL MOLECULE CATABOLIC PROCESS","ORGANIC ACID TRANSPORT","CARBOXYLIC ACID TRANSPORT", "CELLULAR RESPONSE TO EXTERNAL STIMULUS","NEGATIVE REGULATION OF RESPONSE TO STIMULUS","RESPONSE TO ENDOGENOUS STIMULUS","CELLULAR RESPONSE TO ENDOGENOUS STIMULUS","REGULATION OF LIGASE ACTIVITY", "CELLULAR COMPONENT MACROMOLECULE BIOSYNTHETIC PROCESS","REGULATION OF CELLULAR KETONE METABOLIC PROCESS", "POSITIVE REGULATION OF ORGANELLE ORGANIZATION", "RIBONUCLEOPROTEIN COMPLEX BIOGENESIS", "PROTEIN COMPLEX SUBUNIT ORGANIZATION", "PROTEIN COMPLEX BIOGENESIS", "PROTEIN COMPLEX ASSEMBLY", "CELLULAR PROTEIN COMPLEX ASSEMBLY", "RIBONUCLEOPROTEIN COMPLEX SUBUNIT ORGANIZATION", "RIBONUCLEOPROTEIN COMPLEX ASSEMBLY", "REGULATION OF PROTEIN COMPLEX ASSEMBLY", "PROTEIN COMPLEX DISASSEMBLY", "RIBONUCLEOPROTEIN COMPLEX LOCALIZATION", "RIBONUCLEOPROTEIN COMPLEX EXPORT FROM NUCLEUS", "CELLULAR PROTEIN COMPLEX DISASSEMBLY", "REGULATION OF PROTEIN COMPLEX DISASSEMBLY", "PROTEIN COMPLEX LOCALIZATION", "POSITIVE REGULATION OF PROTEIN COMPLEX ASSEMBLY", "CELLULAR PROTEIN COMPLEX LOCALIZATION", "NEGATIVE REGULATION OF PROTEIN COMPLEX DISASSEMBLY", "NEGATIVE REGULATION OF PROTEIN COMPLEX ASSEMBLY", "SMALL NUCLEOLAR RIBONUCLEOPROTEIN COMPLEX ASSEMBLY", "RIBONUCLEOPROTEIN COMPLEX DISASSEMBLY", "CHAPERONE-MEDIATED PROTEIN COMPLEX ASSEMBLY", "POSITIVE REGULATION OF PROTEIN COMPLEX DISASSEMBLY", "NEGATIVE REGULATION OF MACROMOLECULE BIOSYNTHETIC PROCESS", "CELLULAR COMPONENT MOVEMENT", "CELLULAR COMPONENT DISASSEMBLY", "REGULATION OF CELLULAR COMPONENT SIZE", "CELLULAR COMPONENT MAINTENANCE", "REGULATION OF CELLULAR COMPONENT BIOGENESIS", "CELLULAR COMPONENT DISASSEMBLY AT CELLULAR LEVEL", "CELLULAR COMPONENT MAINTENANCE AT CELLULAR LEVEL", "NEGATIVE REGULATION OF CELLULAR METABOLIC PROCESS", "RESPONSE TO ORGANIC SUBSTANCE", "CELLULAR CHEMICAL HOMEOSTASIS", "CHEMICAL HOMEOSTASIS", "REGULATION OF RESPONSE TO STIMULUS", "POSITIVE REGULATION OF RESPONSE TO STIMULUS"),
#                           bp.lenient=c("POSITIVE REGULATION OF MACROMOLECULE METABOLIC PROCESS", "REGULATION OF CELLULAR COMPONENT ORGANIZATION", "POSITIVE REGULATION OF METABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR METABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR PROCESS", "REGULATION OF CELLULAR PROCESS", "REGULATION OF CATALYTIC ACTIVITY", "POSITIVE REGULATION OF CATALYTIC ACTIVITY", "REGULATION OF MOLECULAR FUNCTION", "POSITIVE REGULATION OF CELLULAR COMPONENT ORGANIZATION", "REGULATION OF ORGANELLE ORGANIZATION", "POSITIVE REGULATION OF CATABOLIC PROCESS", "POSITIVE REGULATION OF CELLULAR CATABOLIC PROCESS", "POSITIVE REGULATION OF MOLECULAR FUNCTION", "REGULATION OF CATABOLIC PROCESS", "REGULATION OF CELLULAR CATABOLIC PROCESS", "CELLULAR RESPONSE TO CHEMICAL STIMULUS", "CELLULAR RESPONSE TO ORGANIC SUBSTANCE", "POSITIVE REGULATION OF BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF CELLULAR BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF MACROMOLECULE BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF CELLULAR BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF CELLULAR MACROMOLECULE BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF MACROMOLECULE METABOLIC PROCESS", "RESPONSE TO EXTERNAL STIMULUS", "RESPONSE TO EXTRACELLULAR STIMULUS", "CELLULAR HOMEOSTASIS", "HOMEOSTATIC PROCESS", "REGULATION OF HOMEOSTATIC PROCESS", "ORGANIC SUBSTANCE TRANSPORT", "ORGANIC ACID BIOSYNTHETIC PROCESS", "NEGATIVE REGULATION OF ORGANELLE ORGANIZATION", "ORGANELLE FISSION", "NEGATIVE REGULATION OF CELLULAR COMPONENT ORGANIZATION", "CELLULAR DEVELOPMENTAL PROCESS", "MAINTENANCE OF LOCATION IN CELL","REGULATION OF DEVELOPMENTAL PROCESS","SMALL MOLECULE CATABOLIC PROCESS","ORGANIC ACID TRANSPORT", "CELLULAR RESPONSE TO EXTERNAL STIMULUS","NEGATIVE REGULATION OF RESPONSE TO STIMULUS","RESPONSE TO ENDOGENOUS STIMULUS","CELLULAR RESPONSE TO ENDOGENOUS STIMULUS","REGULATION OF LIGASE ACTIVITY", "CELLULAR COMPONENT MACROMOLECULE BIOSYNTHETIC PROCESS", "POSITIVE REGULATION OF ORGANELLE ORGANIZATION", "NEGATIVE REGULATION OF MACROMOLECULE BIOSYNTHETIC PROCESS", "CELLULAR COMPONENT MOVEMENT", "CELLULAR COMPONENT DISASSEMBLY", "REGULATION OF CELLULAR COMPONENT SIZE", "CELLULAR COMPONENT MAINTENANCE", "REGULATION OF CELLULAR COMPONENT BIOGENESIS", "CELLULAR COMPONENT DISASSEMBLY AT CELLULAR LEVEL", "CELLULAR COMPONENT MAINTENANCE AT CELLULAR LEVEL", "NEGATIVE REGULATION OF CELLULAR METABOLIC PROCESS", "RESPONSE TO ORGANIC SUBSTANCE", "CELLULAR CHEMICAL HOMEOSTASIS", "CHEMICAL HOMEOSTASIS", "REGULATION OF RESPONSE TO STIMULUS", "POSITIVE REGULATION OF RESPONSE TO STIMULUS"),
#                           complexes=c("GOLGI APPARATUS","CELL CORTEX","CELL WALL","CELLULAR BUD","CHROMOSOME","CYTOPLASM","CYTOPLASMIC MEMBRANE-BOUNDED VESICLE","CYTOSKELETON","ENDOMEMBRANE SYSTEM","ENDOPLASMIC RETICULUM","MEMBRANE FRACTION","MEMBRANE","MICROTUBULE ORGANIZING CENTER","MITOCHONDRIAL ENVELOPE","MITOCHONDRION","HETEROGENEOUS NUCLEAR RIBONUCLEOPROTEIN COMPLEX","NUCLEOLUS","NUCLEUS","PEROXISOME","PLASMA MEMBRANE","SITE OF POLARIZED GROWTH","VACUOLE","POLAR MICROTUBULE","SMALL NUCLEAR RIBONUCLEOPROTEIN COMPLEX","SMALL NUCLEOLAR RIBONUCLEOPROTEIN COMPLEX","TRANSCRIPTION FACTOR COMPLEX","CDC73-PAF1 COMPLEX","SIGNAL RECOGNITION PARTICLE", "ARP2-3 PROTEIN COMPLEX", "CCR4-NOT NOT2-NOT5 SUBCOMPLEX", "CDC48-UFD1-NPL4 COMPLEX", "EKC-KEOPS PROTEIN COMPLEX", "HDA COMPLEX", "HRD1 UBIQUITIN LIGASE ERAD-L COMPLEX", "MRNA CAPPING ENZYME COMPLEX", "NRD1-NAB3-SEN1 TERMINATION COMPLEX", "RIC1-RGP1 COMPLEX", "SNRNP U2", "SNRNP U6", "RNA POLYMERASE III COMPLEX"))
# ##############
library(gplots)
##############
mycolors = c(
  "darkorange1",
  "dodgerblue",
  "darkgreen",
  "navy",
  "mediumpurple"  ,
  "royalblue3",
  "darkolivegreen4",
  "firebrick",
  "cyan4",
  "hotpink3",
  "plum4",
  "blue",
  "magenta4",
  "skyblue3",
  "green4",
  "red3",
  "steelblue3",
  "tomato",
  "purple4",
  "goldenrod3",
  "steelblue",
  "darkred",
  "lightpink3",
  "darkorchid",
  "lightblue3",
  "dimgrey",
  "chocolate1",
  "seagreen3",
  "darkkhaki",
  "darksalmon"
)
##############
CLUST.COL <- c("#FF00CC","#33CCFF", "#33CC00", "#9900FF", "#FF9900", "#FFFF00", "#FFCCFF", "#FF0000", "#006600", "#009999", "#CCCC00", "#993300", "#CC99CC", "#6699CC","#CCCCFF", "#FFCC99", "#9966FF", "#CC6600", "#CCFFFF", "#99CC00", "#FF99FF", "#0066FF", "#66FFCC", "#99CCFF", "#9999CC", "#CC9900", "#CC33FF", "#006699", "#F5DF16", "#B5185E", "#99FF00", "#00FFFF", "#990000", "#CC0000", "#33CCCC", "#CC6666", "#996600", "#9999FF", "#3366FF")
rc=col2hex(mycolors)
prunedCol <- "#BEBEBE"
CLUST.COL = c(CLUST.COL,rc)




##### required functions
# computes the number of unique pairs given the number of items to consider
# maxVal - the maximum number of items
# RETURNS a 2-column matrix where each row contains a different pair, specified with item indices
getUniquePairs = function (maxVal)
{
  firstI <- rep(1:(maxVal - 1), (maxVal - 1):1)
  secondI <- sapply(2:maxVal, function(x) {
    x:maxVal
  })
  cbind(firstI, unlist(secondI))
}



##############

# generates a dataframe of gene sets that are *not* significantly enriched, yet they contain query genes,
# i.e. genes in chemical-genetic interactions
# queryGeneSets - named list of queryGeneSets, where each list element is a vector of query genes
#           - the names are filenames (typically) identifying different experiments
# enrichMat - dataframe with enrichment stats for gene sets (one per row), with the following columns:
#             filename, term, geneSetFraction, FDR, overlapGenes, maxOverlapGeneScore
#           - see documentation for the output of hyperG() for descriptions of these columns
#           - rows with the same value, x, in the filename column specify enrichment results for
#             the set of query genes in queryGeneSets with name=x
# scoreMat - score matrix/dataframe; row names are gene IDs and column names are filenames
#          - each column contains a different set of scores
# termsToExclude - vector of terms (i.e. names of gene sets) to exclude from the results; can be NULL
# fdrThresh - FDR threshold; only show gene sets that do not pass this significance threshold
# RETURNS a dataframe of gene sets that are not significantly enriched (one per row),
#         sorted by decreasing maxOverlapGeneScore value. The dataframe includes these columns:
#         filename = filename identifying the query gene set
#         term = gene set name
#         geneSetFraction = the fraction of the term set that overlaps with the query set
#         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#             the scores of the genes are shown in parentheses
#	  maxOverlapGeneScore = the maximum score of the overlapGenes
#         unenrichedGenes = a |-separated list of genes in the overlap of the query set and the term set
#             that *also* do not belong to any significantly enriched term set;
#
# genesNotInEnrichedTerm = function (queryGeneSets, enrichMat, scoreMat, termsToExclude ,
#   fdrThresh = 0.1)
# {
#   scoreMat <- as.matrix(scoreMat)
#   enrichMat <- enrichMat[!(enrichMat$term %in% termsToExclude),
#     , drop = F]
#   lens <- sapply(queryGeneSets, length)
#   queryGeneSets <- queryGeneSets[lens > 0]
#   oGenes <- strsplit(enrichMat$overlapGenes, "\\|")
#   oGenes <- lapply(oGenes, function(genes) {
#     genes <- strsplit(genes, "\\(")
#     sapply(genes, function(vec) {
#       vec[1]
#     })
#   })
#   rowI <- split(1:nrow(enrichMat), enrichMat$filename)
#   enrichI <- match(names(queryGeneSets), names(rowI))
#   extraGenes <- queryGeneSets[is.na(enrichI)]
#   queryGeneSets <- queryGeneSets[!is.na(enrichI)]
#   rowI <- rowI[enrichI[!is.na(enrichI)]]
#   tmp <- lapply(1:length(queryGeneSets), function(expI) {
#     setdiff(queryGeneSets[[expI]], unlist(oGenes[rowI[[expI]]]))
#   })
#   names(tmp) <- names(queryGeneSets)
#   extraGenes <- c(extraGenes, tmp)
#   lens <- sapply(extraGenes, length)
#   extraGenes <- extraGenes[lens > 0]
#   if (length(extraGenes) > 0) {
#     lens <- lens[lens > 0]
#     extraGenes <- data.frame(filename = rep(names(extraGenes),
#       lens), gene = unlist(extraGenes), stringsAsFactors = F)
#     i <- match(extraGenes$gene, rownames(scoreMat))
#     i <- cbind(i, match(extraGenes$filename, colnames(scoreMat)))
#     extraGenes$score <- round(scoreMat[i], 2)
#     extraGenes <- extraGenes[order(extraGenes$score, decreasing = T),
#       ]
#     i <- split(1:nrow(extraGenes), extraGenes$filename)
#     extraGenes <- lapply(i, function(curRow) {
#       tmp <- paste(extraGenes$gene[curRow], "(", extraGenes$score[curRow],
#         ")", sep = "")
#       c(extraGenes$score[curRow[1]], paste(tmp, collapse = "|"))
#     })
#   }
#   tmp <- lapply(1:length(queryGeneSets), function(expI) {
#     curRow <- rowI[[expI]]
#     sigI <- curRow[enrichMat$FDR[curRow] <= fdrThresh]
#     unenrichedGenes <- setdiff(queryGeneSets[[expI]], unlist(oGenes[sigI]))
#     curRow <- setdiff(curRow, sigI)
#     if (length(curRow) == 0) {
#       return(list(rowI = NULL, unenrichedGenes = NULL))
#     }
#     unenrichedGenes <- lapply(oGenes[curRow], function(genes) {
#       intersect(unenrichedGenes, genes)
#     })
#     lens <- sapply(unenrichedGenes, length)
#     unenrichedGenes <- unenrichedGenes[lens > 0]
#     curRow <- curRow[lens > 0]
#     if (length(curRow) == 0) {
#       return(list(rowI = NULL, unenrichedGenes = NULL))
#     }
#     expI <- match(enrichMat$filename[curRow[1]], colnames(scoreMat))
#     unenrichedGenes <- lapply(unenrichedGenes, function(curGenes) {
#       geneI <- match(curGenes, rownames(scoreMat))
#       geneStr <- scoreMat[geneI, expI]
#       names(geneStr) <- curGenes
#       geneStr <- round(sort(geneStr, decreasing = T), 2)
#       geneStr <- paste(names(geneStr), "(", geneStr, ")",
#         sep = "")
#       paste(geneStr, collapse = "|")
#     })
#     list(rowI = curRow, unenrichedGenes = unenrichedGenes)
#   })
#   unenrichedMat <- enrichMat[unlist(lapply(tmp, function(ob) {
#     ob$rowI
#   })), ]
#   unenrichedMat$unenrichedGenes <- unlist(lapply(tmp, function(ob) {
#     ob$unenrichedGenes
#   }))
#   if (length(extraGenes) > 0) {
#     unenrichedMat <- unenrichedMat[c(rep(1, length(extraGenes)),
#       1:nrow(unenrichedMat)), ]
#     toDoI <- 1:length(extraGenes)
#     unenrichedMat$filename[toDoI] <- names(extraGenes)
#     unenrichedMat$term[toDoI] <- "OTHER"
#     unenrichedMat$overlapGenes[toDoI] <- sapply(extraGenes,
#       function(vec) {
#         vec[2]
#       })
#     unenrichedMat$maxOverlapGeneScore[toDoI] <- as.numeric(sapply(extraGenes,
#       function(vec) {
#         vec[1]
#       }))
#     unenrichedMat$unenrichedGenes[toDoI] <- unenrichedMat$overlapGenes[toDoI]
#     if (!is.null(unenrichedMat$pruneOutcome)) {
#       unenrichedMat$pruneOutcome[toDoI] <- "OTHER"
#     }
#     naCol <- setdiff(colnames(unenrichedMat), c("filename",
#       "term", "overlapGenes", "maxOverlapGeneScore", "unenrichedGenes"))
#     colI <- match(naCol, colnames(unenrichedMat))
#     unenrichedMat[toDoI, colI] <- NA
#   }
#   rownames(unenrichedMat) <- NULL
#   unenrichedMat <- unenrichedMat[order(unenrichedMat$geneSetFraction,
#     decreasing = T), ]
#   unenrichedMat[order(unenrichedMat$maxOverlapGeneScore, decreasing = T),
#     ]
# }

# computes enrichment using the hypergeometric test, and uses the resulting P values with
# the Benjamini Hochberg method to estimate FDR values
# querySet - character vector of genes in query set
# geneSets - named list of gene sets to test for significant overlap w/ the query set
# scoreMat - dataframe of gene scores
#          - first column = scores, gene column
#          - can be NULL
# uni - character vector of genes in the universe (i.e. background set)
#     - if NULL, must specify uniSize
# uniSize - the # of genes in the universe
# minSetSize, maxSetSize - min/max # of genes in geneSets (after restricting to the gene universe)
# RETURNS a dataframe of enrichment results, sorted by increasing FDR value. The columns are:
#         term = name of gene set
#         querySetFraction = the fraction of the query set that overlaps with the term set
#         geneSetFraction = the fraction of the term set that overlaps with the query set
#         foldEnrichment = the fold enrichment of the query set with the term genes
#         P = P value estimating the significance with which the query set is enriched with the term genes
#         FDR = FDR value estimating the significance of enrichment
#         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#	  maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
hyperG = function (querySet, geneSets, uni, scoreMat, minSetSize = minGeneSetSize,
  maxSetSize = maxGeneSetSize, uniSize = NA)
{
  if (!is.null(uni)) {
    geneSets <- lapply(geneSets, intersect, uni)
    lens <- sapply(geneSets, length)
    geneSets <- geneSets[lens >= minSetSize & lens <= maxSetSize]
    uniSize <- length(uni)
  }
  if (!is.null(scoreMat)) {
    scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),
      ]
    if (!is.null(uni)) {
      i <- match(uni, scoreMat$gene)
      scoreMat <- scoreMat[sort(i[!is.na(i)]), ]
    }
    scoreMat$score <- round(scoreMat$score, 2)
  }
  enrichInfo <- sapply(geneSets, function(geneSet) {
    overlapSet <- intersect(querySet, geneSet)
    pVal <- phyper(length(overlapSet) - 1, length(geneSet),
      uniSize - length(geneSet), length(querySet), lower.tail = F)
    if (length(overlapSet) > 0) {
      overlapSet <- sort(overlapSet)
    }
    overlapSize <- length(overlapSet)
    if (is.null(scoreMat)) {
      maxScore <- NA
    }
    else {
      i <- sort(match(overlapSet, scoreMat$gene))
      maxScore <- scoreMat$score[i[1]]
      overlapSet <- paste(scoreMat$gene[i], "(", scoreMat$score[i],
        ")", sep = "")
    }
    overlapSet <- paste(overlapSet, collapse = "|")
    bgRate <- length(geneSet)/uniSize
    foldEnrich <- overlapSize/length(querySet)/bgRate
    c(overlapSet, overlapSize/length(geneSet), foldEnrich,
      pVal, maxScore, overlapSize/length(querySet))
  })
  enrichInfo <- t(enrichInfo)
  enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[,
    6]), geneSetFraction = as.numeric(enrichInfo[, 2]), foldEnrichment = as.numeric(enrichInfo[,
      3]), P = as.numeric(enrichInfo[, 4]), FDR = p.adjust(as.numeric(enrichInfo[,
        4]), method = "BH"), overlapGenes = enrichInfo[, 1],
    maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = F)
  rownames(enrichCol) <- NULL
  enrichCol = enrichCol[order(enrichCol$FDR), ]
}
#


#######hiphop:::overlapCoeff
#######the overlap of genesets for all combinations
#######If set X is a subset of Y or the converse then the overlap coefficient is equal to 1.

# compute the overlap coefficient given a pair of (gene) sets
# gsPairList - a list of two sets (each set is a vector of IDs)
# RETURNS the overlap coefficient

overlapCoeff = function (gsPairList)
{
  length(intersect(gsPairList[[1]], gsPairList[[2]]))/min(length(gsPairList[[1]]),
    length(gsPairList[[2]]))
}

####### given enrichInfo after hyperG, computes edgeMat for making enrichment map
# generates an enrichment map in xgmml format using hypergeometric test statistics
# enrichInfo - dataframe with enrichment stats for gene sets (one per row), with the following columns:
#              term, geneSetFraction, querySetFraction, FDR, overlapGenes, maxOverlapGeneScore
#            - see documentation for the output of hyperG() for descriptions of these columns
# geneSets - named list of gene sets tested for significant overlap w/ the query set,
#              restricted to genes in the universe
# outFile - the output xgmml file will be saved to this location
# fdrThresh - FDR threshold; only show gene sets that pass this significance threshold
# overlapThresh - an edge between a pair of enriched gene sets will only be shown if the overlap coefficient
#              is >= overlapThresh
# nonEnrichInfo - dataframe with info on gene sets that are *not* significantly enriched (one per row)
#             with the following columns:
#             term, overlapGenes, maxOverlapGeneScore, geneSetFraction, unenrichedGenes
#            - see documentation for the output of genesNotInEnrichedTerm() for descriptions of these columns
#            - can be NULL
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreName - score label to use in top overlap gene barplots
# plotForEachEnrichedTerm - if TRUE, a top overlap gene barplot will be created for each enriched term;
#              if FALSE, a barplot will be created for each enriched term cluster
# goTable - a dataframe with the following columns describing GO terms:
#         - "term" (GO term), "id" (GOID)
#         - if provided (i.e. not NULL), the GO ID numbers of the enriched GO terms will be saved in
#           the output xgmml file as a node attribute called "GOID".
#         - the GOID allows for an easy link to the GO website page for the associated GO term

###############
#### query set is genes with significant fitness defect, uni is all the genes in the data matrix, the union
clusterEnrich= function(enrichInfo, geneSets, outFile, fdrThresh = 0.1, overlapThresh = 0.5,
                                  nonEnrichInfo = NULL, barModGenes = NULL, scoreName = "Fitness defect score",
                                  plotForEachEnrichedTerm = F,  go_path = NULL, go_input = NULL){




  go_file = file.path(go_path)

  if(!is.null(go_input)) {goTable = go_input} else if(!is.null(go_path)){
goTable = read.delim(go_file,stringsAsFactors = F,check.names = F)}
   #i
  nodeSizeRange <- c(10, 40)
  prunedCol <- "#BEBEBE"
  labelWidth <- 20
  edgeWidthRange <- c(1, 5)
  overlapCoeffRange <- c(overlapThresh, 1)

  if (!is.null(nonEnrichInfo)) {
    nonEnrichInfo$maxOverlapGeneScore <- round(nonEnrichInfo$maxOverlapGeneScore,
                                               2)
    nonEnrichInfo$geneSetFraction <- round(nonEnrichInfo$geneSetFraction *100, 1)
    if (is.null(nonEnrichInfo$nGenes)) {
      lens <- sapply(geneSets, length)
      i <- match(nonEnrichInfo$term, names(lens))
      nonEnrichInfo$nGenes <- lens[i]
    }
    tmp <- strsplit(nonEnrichInfo$overlapGenes, "\\|")
    w <- which(is.na(tmp))
    if(length(w)>0) tmp = tmp[-w]
    if (is.null(nonEnrichInfo$unenrichedGenes)) {
      nonEnrichInfo$overlapGenes <- sapply(tmp, paste,
                                           collapse = "| ")
    }
    else {
      unEnriched <- strsplit(nonEnrichInfo$unenrichedGenes,
                             "\\|")
      tmp.mod <- sapply(1:length(tmp), function(termI) {
        vec <- tmp[[termI]]
        i <- match(unEnriched[[termI]], tmp[[termI]])
        vec[i] <- paste("<b>", vec[i], "</b>", sep = "")
        paste(vec, collapse = "| ")
      })

      nonEnrichInfo$overlapGenes <- tmp.mod
    }

    if (is.null(enrichInfo)) {

      return()
    }
  }

  enrichInfo <- enrichInfo[enrichInfo$FDR <= fdrThresh, , drop = F]

  enrich = enrichInfo

  if(!is.null(go_path)) {
    toDoI <- match(enrichInfo$term, goTable$term)
    enrichInfo$GOID = goTable$GOID[toDoI]
  }

  if (nrow(enrichInfo) == 0) {

    return()
  }
  enrichInfo$formattedLabel <- sapply(enrichInfo$term, function(curLabel) {
    curLabel <- strwrap(curLabel, labelWidth)
    paste(curLabel, collapse = "\n")
  })
  i <- match(enrichInfo$term, names(geneSets))
  if (any(is.na(i))) {
    stop("Could not find gene sets for ", sum(is.na(i)),
         " enriched terms.")
  }
  geneSets <- geneSets[i]
  if (is.null(enrichInfo$nGenes)) {
    enrichInfo$nGenes <- sapply(geneSets, length)
  }
  tmpSize <- -log10(enrichInfo$FDR)
  maxVal <- max(tmpSize[!is.infinite(tmpSize)])
  tmpSize[is.infinite(tmpSize)] <- maxVal + 2
  gsSizeRange <- range(tmpSize)
  if (gsSizeRange[1] == gsSizeRange[2]) {
    gsSizeRange[1] <- -log10(fdrThresh)
    gsSizeRange[2] <- gsSizeRange[2] + 1
  }
  tmpSize <- (tmpSize - gsSizeRange[1])/(gsSizeRange[2] - gsSizeRange[1])
  tmpSize <- nodeSizeRange[1] + tmpSize * (nodeSizeRange[2] -nodeSizeRange[1])
  enrichInfo$size <- round(tmpSize, 2)
  if (nrow(enrichInfo) == 1) {
    enrichInfo$cluster <- CLUST.COL[1]
    edgeMat <- NULL
  }



#########################   EDGE MAT   ########################
  else {
    pairI <- getUniquePairs(length(geneSets))
    distVal <- apply(pairI, 1, function(onePair) {
      overlapCoeff(geneSets[onePair])
    })
    distVal[distVal < overlapThresh] <- 0
    edgeMat <- data.frame(nodeA = pairI[, 1], nodeB = pairI[,2], coeff = distVal)
    enrichInfo$cluster <- prunedCol
    if (is.null(enrichInfo$pruneOutcome)) {
      termI <- 1:nrow(enrichInfo)
    }
    else {
      termI <- which(enrichInfo$pruneOutcome == enrichInfo$term)
    }
    if (length(termI) == 1) {
      enrichInfo$cluster[termI] <- CLUST.COL[1]
    }
    else {
      i <- which((edgeMat$nodeA %in% termI) & (edgeMat$nodeB %in%
                                                 termI))
      enrichInfo$id = termI

      g=igraph::graph_from_data_frame(edgeMat[which(edgeMat$coeff!=0),],directed = F,vertices = enrichInfo$id)
      adj = igraph::as_adjacency_matrix(g)
      clusters = igraph::clusters(g)
      clusters = split(names(clusters$membership),clusters$membership)


      clusters <- lapply(clusters, as.numeric)

      lens <- sapply(clusters, length)
      clusters <- data.frame(id = unlist(clusters), cluster = CLUST.COL[rep(1:length(clusters),
                                                                            lens)], stringsAsFactors = F)
      enrichInfo$cluster[clusters$id] <- clusters$cluster
    }
    edgeMat <- edgeMat[edgeMat$coeff > 0, , drop = F]
    if (nrow(edgeMat) > 0) {
      edgeMat$size <- (edgeMat$coeff - overlapCoeffRange[1])/(overlapCoeffRange[2] -
                                                                overlapCoeffRange[1])
      edgeMat$size <- edgeWidthRange[1] + edgeMat$size *
        (edgeWidthRange[2] - edgeWidthRange[1])
      edgeMat$coeff <- round(edgeMat$coeff, 2)
      edgeMat$size <- round(edgeMat$size, 2)
    }
    else {
      edgeMat <- NULL
    }
  }
  otherI <- order(enrichInfo$cluster)
  otherI <- otherI[order(enrichInfo$FDR[otherI])]
  termI <- which(enrichInfo$cluster[otherI] != prunedCol)
  if (length(termI) < length(otherI)) {
    otherI <- c(otherI[termI], otherI[-termI])
  }
  enrichInfo$id <- 1:nrow(enrichInfo)
  enrichInfo <- enrichInfo[otherI, , drop = F]
  enrichInfo$geneSetFraction <- round(enrichInfo$geneSetFraction *
                                        100, 1)
  enrichInfo$querySetFraction <- round(enrichInfo$querySetFraction *
                                        100, 1)


  #if (is.null(edgeMat)) print("edgeMat is NULL")


  if (!is.null(edgeMat)) {
    nam = c("source","target","label","overlapCoeff","width")
    orig = c("nodeA","nodeB","label","coeff","size")

    src = names(geneSets)[edgeMat$nodeA]

    trg = names(geneSets)[edgeMat$nodeB]
    edgeMat$label = paste(src,"(overlap)",trg)
    m = match(names(edgeMat),orig)
    names(edgeMat) = nam[m]
  }


  output = list(enrichInfo = enrichInfo,edgeMat = edgeMat)

  return(output)
}
###############

compSCORE <- function(mat,coln, sig = 1){

  df = data.frame(score = mat[,coln],stringsAsFactors = F)
  library(dplyr)
  df$gene = rownames(mat)
  rownames(df) = df$gene
  df$index=0
  wdf = which(df$score >= sig)
  df$index[wdf]=1
  df = df[,c('index','score','gene')]
  df = df %>% arrange(desc(score))
  df
}

  ###### bulk of code

  score = compSCORE(mat,coln,sig = sig)

  fdrThresh = as.numeric(fdrThresh)

  bp_file = file.path(bp_path)
  scoreMat = score


  queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index == 1)]))
  uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))

  if(!is.null(bp_input))  {bp = bp_input} else {bp <- readRDS(bp_file)}

  #uniGenes.mn <- unique(intersect(scoreMat$gene,unlist(bp,use.names=F)))

  #### intersect of the geneSets with the backgroundSet, filtering for size


  enrichMat.mn <- hyperG(querySet = queryGenes.mn, geneSets = bp,
                           uni = uniGenes.mn, scoreMat = score, minSetSize = minGeneSetSize,
                           maxSetSize = maxGeneSetSize, uniSize = NA)
  queryGeneSets = list()
  queryGeneSets[[curr_exp]] = queryGenes.mn
  enrichMat.mn$filename <- curr_exp
  enrichMat_Ordered = enrichMat.mn[with(enrichMat.mn, order(FDR,
                                                            -foldEnrichment)), ]
  scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),]
  scoreMat <- scoreMat[match(uniGenes.mn, scoreMat$gene), "score",drop = F]
  rownames(scoreMat) <-   uniGenes.mn
  colnames(scoreMat) <- curr_exp
#
#   nonEnrichMat.mn <- genesNotInEnrichedTerm(queryGeneSets,
#                                               enrichMat.mn, scoreMat, NONSPECIFIC.TERMS$bp, fdrThresh)
#
  maxSetSize = maxGeneSetSize
  #### intersect of the geneSets with the backgroundSet, filtering for size
  bp <- lapply(bp, intersect, uniGenes.mn)
  lens <- sapply(bp, length)
  bp <- bp[lens >= minGeneSetSize & lens <= maxGeneSetSize]

  q = clusterEnrich(enrichInfo = enrichMat.mn, geneSets = bp,
                      fdrThresh = fdrThresh, overlapThresh = 0.5,
                      nonEnrichInfo = NULL, barModGenes = NULL,
                      scoreName = "score", plotForEachEnrichedTerm = F,go_path = go_path)

  edgeMat = q$edgeMat
  enrichInfo = q$enrichInfo
  library(dplyr)
  if(!is.null(enrichInfo)) {

    enrichInfo = enrichInfo %>% arrange(FDR)
    enrichInfo$nOverlap = round(enrichInfo$geneSetFraction*enrichInfo$nGenes/100)
    enrichInfo$nQuery = round((enrichInfo$geneSetFraction*enrichInfo$nGenes/100)/(enrichInfo$querySetFraction/100))

    w = which(names(enrichInfo) %in% c("querySetFraction", "geneSetFraction" ,
                                   "foldEnrichment" , "P" , "FDR" ))
    enrichInfo[,c("querySetFraction","geneSetFraction", "foldEnrichment")] =
      round(enrichInfo[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2)

    enrichInfo[,c("P", "FDR")] =
      signif(enrichInfo[,c("P", "FDR")],digits = 3)
    nam = c(
      "filename",
      "GOID",
      "term",
      "nGenes"  ,
      "nQuery",
      "nOverlap"     ,
      "querySetFraction",
      "geneSetFraction" ,
      "foldEnrichment"  ,
      "P"            ,
      "FDR"      ,
      "overlapGenes",
      "maxOverlapGeneScore" ,
      "cluster"      ,
      "id" ,
      "size"    ,
      "formattedLabel"
    )
      enrichInfo = enrichInfo[,nam]


  }else { print("no GO enrichment!") }

  return(list(enrichInfo = enrichInfo , edgeMat = q$edgeMat))
}
####################################
####################################

####################################
####################################


