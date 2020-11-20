get.raw.edges <- function(raw_edges, file, sep, header) {
  if(!is.na(file) & first(!is.na(raw_edges))) {
    stop("Both raw_edges and file are specified, specify only one")
  }

  if(!is.na(file)) {
    if(is.na(raw_edges)) {
      raw_edges <- read.table(file = file, sep = sep, header = header, stringsAsFactors = F, quote = "")
      if(ncol(raw_edges) < 2) {stop("File contains less than 2 columns, check that \"sep\" value is right")}
    } else {
      stop("No file or edgelist specified, please see manual for usage")
    }
  } else {
    raw_edges[] <- apply(raw_edges, 2, as.character)
  }

  return(raw_edges)
}

#' Find balanced coloring using Kamei algorithm
#'
#' This function finds a balanced coloring. Specify only one file or raw_edges.
#'
#' @param raw_edges List of edges
#' @param file File from which edgelist can be read (specify only one file or raw_edges)
#' @param sep Field separator character used in read.table. If missing is set at " ".
#' @param header A logical value indicating whether "file" contains the names of the variables as its first line. If missing is set at F
#' @param directed A logical value indicating whether the network is directed or not. If missing is set at F
#' @param look.for.no.input.nodes A logical value indicating whether colors of nodes with no inputs need to be fixed. If missing is set at T
#' @return Data frame with 2 columns: Name (for the node name) and Color (for the node color id)
#' @export
get.balanced.coloring.Kamei <- function(raw_edges = NA, file = NA, sep = " ", header = F, directed = F, look.for.no.input.nodes = T) {
  ####################################
  # Preprocessing of nodes and edges #
  ####################################
  raw_edges = get.raw.edges(raw_edges = raw_edges, file = file, sep = sep, header = header)
  weighted = as.logical(ncol(raw_edges) - 2)

  graph = igraph::graph_from_edgelist(as.matrix(raw_edges[, 1:2]), directed = directed)

  integer.nodes = igraph::as_data_frame(graph, what = "vertices")
  if(nrow(integer.nodes) > 30000) {
    stop(paste0("Number of nodes (", nrow(integer.nodes), ") is too big to run this version of the code due to possible memory overflow"))
  }
  integer.nodes$Id = 1:nrow(integer.nodes)
  colnames(integer.nodes)[1] = "Name"
  integer.nodes = integer.nodes[, c("Id", "Name")]

  integer.edges = data.frame(igraph::as_edgelist(graph, names = F), stringsAsFactors = F)
  colnames(integer.edges) = c("Source", "Target")

  integer.weights = NULL
  if(weighted == T) {
    integer.weights = unique(raw_edges[, 3])
    integer.edges$Weight = as.integer(factor(raw_edges[, 3], levels = integer.weights))
  }
  ####################################
  ####### End of preprocessing #######
  ####################################

  integer.nodes$Color = 1

  no.input.nodes = NULL
  if(directed == T & look.for.no.input.nodes == T) {
    initial.colors <- (1:nrow(integer.nodes))[!1:nrow(integer.nodes) %in% unique(dplyr::filter(integer.edges, Source != Target)[, 2])]
    no.input.nodes <- (1:nrow(integer.nodes))[!1:nrow(integer.nodes) %in% unique(integer.edges[, 2])]
    integer.nodes$Color[initial.colors] <- 2:(length(initial.colors) + 1)
    integer.nodes$Fixed = 0
    if(length(no.input.nodes) != 0) {
      integer.nodes$Fixed[no.input.nodes] <- 1
    }
  }

  integer.nodes$Color <- getBalancedColoring(integer.nodes$Color, integer.nodes$Fixed, as.matrix(integer.edges), directed, weighted, length(integer.weights))
  return(dplyr::select(integer.nodes, -c(Id, Fixed)))
}

#' Find the input set color vector given the coloring
#'
#' This function finds the input set color vector given the coloring. Note, can be undirected
#'
#' @param nodes dataframe with 2 columns: Name, Color.
#' @param raw_edges List of edges
#' @param file File from which edgelist can be read (specify only one file or raw_edges)
#' @param sep Field separator character used in read.table. If missing is set at " ".
#' @param header A logical value indicating whether "file" contains the names of the variables as its first line. If missing is set at F
#' @param directed A logical value indicating whether the network is directed or not. If missing is set at F
#' @return
#' @export
get.input.set.color.vector <- function(nodes, raw_edges = NA, file = NA, sep = " ", header = F, directed = F) {
  ####################################
  # Preprocessing of nodes and edges #
  ####################################
  raw_edges = get.raw.edges(raw_edges = raw_edges, file = file, sep = sep, header = header)
  weighted = as.logical(ncol(raw_edges) - 2)

  graph = igraph::graph_from_edgelist(as.matrix(raw_edges[, 1:2]), directed = directed)

  integer.nodes = igraph::as_data_frame(graph, what = "vertices")
  if(nrow(integer.nodes) > 30000) {
    stop(paste0("Number of nodes (", nrow(integer.nodes), ") is too big to run this version of the code due to possible memory overflow"))
  }
  integer.nodes$Id = 1:nrow(integer.nodes)
  colnames(integer.nodes)[1] = "Name"
  integer.nodes = integer.nodes[, c("Id", "Name")]

  integer.edges = data.frame(igraph::as_edgelist(graph, names = F), stringsAsFactors = F)
  colnames(integer.edges) = c("Source", "Target")

  integer.weights = NULL
  if(weighted == T) {
    integer.weights = unique(raw_edges[, 3])
    integer.edges$Weight = as.integer(factor(raw_edges[, 3], levels = integer.weights))
  }
  ####################################
  ####### End of preprocessing #######
  ####################################

  integer.nodes <- merge(x = integer.nodes, y = nodes, by = "Name", sort = F)

  nodes
  vectors <-
    calculateVectors(nodeColors = integer.nodes$Color,
                     edges = as.matrix(integer.edges),
                     directed = directed,
                     weighted = weighted,
                     numberOfWeights = length(integer.weights),
                     fromR = T)
  vectors <- matrix(unlist(vectors), ncol = nrow(nodes))
  vectors <- data.frame(vectors, stringsAsFactors = F)
  colnames(vectors) <- integer.nodes$Name

  vectors <- vectors[, nodes$Name]
  return(vectors)
}
