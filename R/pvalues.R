progress.bar <- function (x, max = 100, startTime, nowTime) {
  timePerPercent = difftime(nowTime, startTime, units ="secs") / x
  percent <- x / max * 100
  cat(sprintf("\r[%-50s] %d%% Est time remaining: %d s",
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent), floor((max - x) * timePerPercent)
              )
      )
  if (x == max)
    cat('\n')
}

get.shuffled.edges.degreeSequence <- function(raw_edges) {
  graph <- igraph::graph_from_edgelist(raw_edges)
  outDegrees <- igraph::degree(graph = graph, mode = "out", loops = T, normalized = FALSE)
  outDegrees <- outDegrees[sample(igraph::vcount(graph))]
  inDegrees <- igraph::degree(graph = graph, mode = "in", loops = T, normalized = FALSE)
  inDegrees <- inDegrees[sample(igraph::vcount(graph))]
  
  newGraph <- igraph::sample_degseq(out.deg = outDegrees, in.deg = inDegrees, method = "simple")
  newEdges <- as.data.frame(igraph::as_edgelist(newGraph), stringsAsFactors = F)
  newEdges[] <- apply(newEdges, 2, as.character)
  
  return(newEdges)
}

get.shuffled.edges.erdosrenyi <- function(raw_edges) {
  graph <- igraph::sample_gnm(n = length(unique(c(raw_edges[, 1], raw_edges[, 2]))), m = nrow(raw_edges), directed = T, loops = T)
  newEdges <- as.data.frame(igraph::as_edgelist(graph), stringsAsFactors = F)
  newEdges[] <- apply(newEdges, 2, as.character)
  return(newEdges)
}

get.shuffled.edges <- function(raw_edges, method = "degreeSequence") {
  if(method == "degreeSequence") {
    return(get.shuffled.edges.degreeSequence(raw_edges))
  }
  if(method == "erdosrenyi") {
    return(get.shuffled.edges.erdosrenyi(raw_edges))
  }
}

get.new.edges <- function(raw_edges, method = c("degreeSequence", "erdosrenyi")) {
  if(length(method) > 1) {
    stop("Specify randomization method to run this function")
  }
  if(ncol(raw_edges) == 2) {
    newEdges <- get.shuffled.edges(raw_edges, method = method)
  } else {
    newEdges = data.frame(V1 = character(),
                          V2 = character(),
                          V3 = character(),
                          stringsAsFactors = FALSE)
    
    uniqueWeights = unique(raw_edges$V3)
    for(i in 1:length(uniqueWeights)) {
      toAdd = get.shuffled.edges(raw_edges = as.matrix(raw_edges[raw_edges$V3 == uniqueWeights[i], 1:2]), method = method)
      toAdd$V3 = uniqueWeights[i]
      newEdges = rbind(newEdges, toAdd)
    }
  }
  return(newEdges)
}

zscore.to.pvalue <- function(x) {
  return(1 - (pracma::erf(x) + 1) / 2)
}

#' Get p-values corresponding to fiber building blocks.
#'
#' get.building.block.pvalues() function finds fiber building block p-values.
#' That is, the probability (estimated on the sample of the size "sampleSize") of having the number of building blocks
#' in the graph (specified using "raw_edges", "file", "sep" and "header" variables)
#' with the certain fiber numbers or belonging to the certain class or having the certain block name (depending on the "mode" variable)
#' in the random network (with the same degree sequence or created by Erdos-Renyi depending on "method" variable).
#'
#' @param raw_edges 2 or 3 column data frame specifying the list of edges (specify only one file or raw_edges)
#' @param file Path to the file with the edgelist. Make sure to specify "sep" (if different from " ") and "header" (if different from FALSE) to be passed to the read.table function.
#' @param sep To be used with the "file" variable. Defines the field separator character to be used in the read.table() function. Is set with " " by default.
#' @param header To be used with the "file" variable. A logical value indicating whether "file" contains the names of the variables as its first line. Is set as FALSE by default.
#' @param sampleSize Size of the sample on which to find p-values. Default is 10000.
#' @param mode Which p-values need to be found: nl, class or block name.
#' @param method Randomization method: "degreeSequence" or "erdosrenyi".
#' @param progressBar A logical value indicating whether to show the progress bar.
#' @return The count of building blocks in the network, their occurence in the random network and corresponding Z-Scores and p-values.
#' @export
get.building.block.pvalues <- function(raw_edges = NA, file = NA, sep = " ", header = F, sampleSize = 10000, mode = c("nl", "Class", "BlockName"),
                                       method = c("degreeSequence", "erdosrenyi"), progressBar = T) {
  if(length(method) > 1) {
    stop("Specify randomization method to run this function")
  }
  raw_edges = get.raw.edges(raw_edges = raw_edges, file = file, sep = sep, header = header)
  weighted = as.logical(ncol(raw_edges) - 2)
  buildingBlocks = get.building.blocks(raw_edges = raw_edges, progressBar = F)
  buildingBlocks = buildingBlocks %>%
    dplyr::group_by_at(mode) %>%
    dplyr::summarise(Count = n())
  colnames(buildingBlocks)[1] = "Mode"
  
  modelSummary = data.frame(Mode = character(),
                            Count = character(),
                            Trial = character(),
                            stringsAsFactors = FALSE)
  
  start.time = Sys.time()
  for(i in 1:sampleSize) {
    if(progressBar == T) {
      progress.bar(i, max = sampleSize, startTime = start.time, nowTime = Sys.time())
    }
    syntheticBlocks = get.building.blocks(raw_edges = get.new.edges(raw_edges, method = method), progressBar = F)
    
    syntheticBlocks =
      syntheticBlocks %>%
      dplyr::group_by_at(mode) %>%
      dplyr::summarise(Count = n())
    colnames(syntheticBlocks)[1] = "Mode"
    
    syntheticBlocks$Trial = i
    
    modelSummary = rbind(modelSummary, syntheticBlocks)
  }
  
  modelSummary = tidyr::spread(modelSummary, key = Mode, value = Count)
  modelSummary[is.na(modelSummary)] = 0
  
  modelSummary =
    data.frame(Mode = colnames(modelSummary[, -1]),
               MeanRandom = apply(modelSummary[, -1], 2, mean),
               SDRandom = apply(modelSummary[, -1], 2, sd))
  
  buildingBlocks = merge(buildingBlocks, modelSummary, by = "Mode", all.x = T)
  buildingBlocks[is.na(buildingBlocks)] = 0
  
  buildingBlocks$ZScore = round((buildingBlocks$Count - buildingBlocks$MeanRandom) / buildingBlocks$SDRandom, digits = 3)
  buildingBlocks$PValue = round(zscore.to.pvalue(buildingBlocks$ZScore), digits = 3)
  
  colnames(buildingBlocks)[1] = mode
  
  return(buildingBlocks)
}