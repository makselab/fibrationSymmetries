get.shuffled.edges <- function(oldEdges) {
  graph <- igraph::graph_from_edgelist(oldEdges)
  outDegrees <- igraph::degree(graph = graph, mode = "out", loops = T, normalized = FALSE)
  outDegrees <- outDegrees[sample(igraph::vcount(graph))]
  inDegrees <- igraph::degree(graph = graph, mode = "in", loops = T, normalized = FALSE)
  inDegrees <- inDegrees[sample(igraph::vcount(graph))]
  
  newGraph <- igraph::sample_degseq(out.deg = outDegrees, in.deg = inDegrees, method = "simple")
  newEdges <- as.data.frame(igraph::as_edgelist(newGraph), stringsAsFactors = F)
  newEdges[] <- apply(newEdges, 2, as.character)
  
  return(newEdges)
}

get.new.edges <- function(raw_edges, method = c("degreeSequence")) {
  if(ncol(raw_edges) == 2) {
    newEdges <- get.shuffled.edges(raw_edges)
  } else {
    newEdges = data.frame(V1 = character(),
                          V2 = character(),
                          V3 = character(),
                          stringsAsFactors = FALSE)
    
    uniqueWeights = unique(raw_edges$V3)
    for(i in 1:length(uniqueWeights)) {
      toAdd = get.shuffled.edges(oldEdges = as.matrix(raw_edges[raw_edges$V3 == uniqueWeights[i], 1:2]))
      toAdd$V3 = uniqueWeights[i]
      newEdges = rbind(newEdges, toAdd)
    }
  }
  return(newEdges)
}

zscore.to.pvalue <- function(x) {
  return(1 - (pracma::erf(x) + 1) / 2)
}

#' Find fibration building blocks and their p-values
#'
#' This function finds fibration building blocks Specify only one file or raw_edges.
#'
#' @param raw_edges List of edges
#' @param file File from which edgelist can be read (specify only one file or raw_edges)
#' @param sep Field separator character used in read.table. If missing is set at " ".
#' @param header A logical value indicating whether "file" contains the names of the variables as its first line. If missing is set at F.
#' @param sampleSize Size of the sample on which to find p-values. Default is 10000.
#' @param mode Which p-values need to be found: nl (by fiber numbers), class (by class) or by block name.
#' @return The count of building blocks in the network, their occurence in the random network and corresponding Z-Scores and p-values.
#' @export
get.building.block.pvalues <- function(raw_edges = NA, file = NA, sep = " ", header = F, sampleSize = 10000, mode = c("nl", "class", "blockName")) {
  if(mode != "nl") {
    stop("Only nl mode is supported currently")
  }
  raw_edges = get.raw.edges(raw_edges = raw_edges, file = file, sep = sep, header = header)
  weighted = as.logical(ncol(raw_edges) - 2)
  buildingBlocks = get.building.blocks(raw_edges = raw_edges)
  buildingBlocks = buildingBlocks %>%
    dplyr::group_by(nl) %>%
    dplyr::summarise(Count = n())
  
  modelSummary = data.frame(nl = character(),
                            Count = character(),
                            Trial = character(),
                            stringsAsFactors = FALSE)
  for(i in 1:sampleSize) {
    syntheticBlocks = get.building.blocks(raw_edges = get.new.edges(raw_edges))
    
    syntheticBlocks =
      syntheticBlocks %>%
      dplyr::group_by(nl) %>%
      dplyr::summarise(Count = n())
    
    syntheticBlocks$Trial = i
    
    modelSummary = rbind(modelSummary, syntheticBlocks)
  }
  
  modelSummary = tidyr::spread(modelSummary, key = nl, value = Count)
  modelSummary[is.na(modelSummary)] = 0
  
  modelSummary =
    data.frame(nl = colnames(modelSummary[, -1]),
               MeanRandom = apply(modelSummary[, -1], 2, mean),
               SDRandom = apply(modelSummary[, -1], 2, sd))
  
  buildingBlocks = merge(buildingBlocks, modelSummary, by = "nl", all.x = T)
  buildingBlocks[is.na(buildingBlocks)] = 0
  
  buildingBlocks$ZScore = round((buildingBlocks$Count - buildingBlocks$MeanRandom) / buildingBlocks$SDRandom, digits = 3)
  buildingBlocks$PValue = round(zscore.to.pvalue(buildingBlocks$ZScore), digits = 3)
  
  return(buildingBlocks)
}