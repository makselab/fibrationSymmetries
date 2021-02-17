get.raw.edges <- function(raw_edges, file, sep, header) {
  if(!is.na(file) & dplyr::first(!is.na(raw_edges))) {
    stop("Both raw_edges and file are specified, specify only one")
  }

  if(!is.na(file)) {
    if(is.na(raw_edges)) {
      raw_edges <- read.table(file = file, sep = sep, header = header, stringsAsFactors = F, quote = "")
      if(ncol(raw_edges) < 2) {stop("File contains less than 2 columns, check that \"sep\" value is right")}
      if(header == T) {
        colnames(raw_edges) = paste0("V", 1:ncol(raw_edges))
      }
    } else {
      stop("No file or edgelist specified, please see manual for usage")
    }
  } else {
    raw_edges[] <- apply(raw_edges, 2, as.character)
  }
  
  return(raw_edges)
}

#' Find the minimal balanced coloring of the graph using the algorithm introduced by Kamei and Cock.
#'
#' The get.balanced.coloring.Kamei() function finds a minimal balanced coloring based on the algorithm developed in (1).
#' This function supports directed and undirected cases as well as weighted and unweighted cases.
#' Directionality is defined by the variable "directed" and weightedness is defined by the number of columns in the input data (2 for undirected and 3 for directed).
#' Input graph is specified using "raw_edges", "file", "sep" and "header" variables.
#' (1) Kamei H, Cock PJ. A. Computation of balanced equivalence relations and their lattice for a coupled cell network. SIAM J Appl Dyn Syst. 2013;12,352-382.
#'
#' @param raw_edges 2 or 3 column data frame specifying the list of edges (specify only one file or raw_edges)
#' @param file Path to the file with the edgelist. Make sure to specify "sep" (if different from " ") and "header" (if different from FALSE) to be passed to the read.table function.
#' @param sep To be used with the "file" variable. Defines the field separator character to be used in the read.table() function. Is set with " " by default.
#' @param header To be used with the "file" variable. A logical value indicating whether "file" contains the names of the variables as its first line. Is set as FALSE by default.
#' @param directed A logical value indicating whether the network is directed or not. Is set as FALSE by default.
#' @param look.for.no.input.nodes Only valid for the directed case.
#' Coloring is balanced if all the nodes of the same color have same ISCVs (Input Set Color Vectors as defined in SI VI of Leifer I, Morone F, Reis SDS, Andrade JS, Sigman M, Makse HA. Circuits with broken fibration symmetries perform core logic computations in biological networks).
#' Obvious problem with this definition is that all the nodes that have no inputs have the same ISCVs and therefore have to be the same color.
#' This may cause a cascade effect where nodes that receive from nodes with no inputs will also have the color same to each other.
#' There is no reason for the nodes that are assigned the same color as a result of this problem to be synchronized, because the information they receive is different.
#' To fight this, the algorithm will force all the nodes with no inputs to have different colors at every iteration.
#' To turn this functionality off set "look.for.no.input.nodes" as FALSE. Make sure you understand what you are doing.
#' @return Data frame with 2 columns: Name (for the node name) and Color (for the node color id corresponding to the minimal balanced coloring).
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
  integer.nodes$Color = 1
  integer.nodes$Fixed = 0

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

  no.input.nodes = NULL
  if(directed == T & look.for.no.input.nodes == T) {
    initial.colors <- (1:nrow(integer.nodes))[!1:nrow(integer.nodes) %in% unique(dplyr::filter(integer.edges, Source != Target)[, 2])]
    no.input.nodes <- (1:nrow(integer.nodes))[!1:nrow(integer.nodes) %in% unique(integer.edges[, 2])]
    integer.nodes$Color[initial.colors] <- 2:(length(initial.colors) + 1)
    if(length(no.input.nodes) != 0) {
      integer.nodes$Fixed[no.input.nodes] <- 1
    }
  }

  integer.nodes$Color <- getBalancedColoring(integer.nodes$Color, integer.nodes$Fixed, as.matrix(integer.edges), directed, weighted, length(integer.weights))
  return(dplyr::select(integer.nodes, -c(Id, Fixed)))
}

#' Find the ISCV (Input Set Color Vector) of all nodes in the graph given the coloring.
#'
#' The get.input.set.color.vector() function finds ISCVs (Input Set Color Vectors as defined in SI VI of (1)) of all nodes in the graph given the coloring of the graph.
#' Coloring is specified by "nodes" variable formatted as the output of the get.balanced.coloring.Kamei() function.
#' Input graph is specified using "raw_edges", "file", "sep" and "header" variables.
#' (1) Leifer I, Morone F, Reis SDS, Andrade JS, Sigman M, Makse HA. Circuits with broken fibration symmetries perform core logic computations in biological networks. PLoS Comput Biol 2020;16(6):e1007776.
#'
#' @param nodes Variable is formatted as the output of the get.balanced.coloring.Kamei() function. Dataframe with 2 columns: Name, Color.
#' @param raw_edges 2 or 3 column data frame specifying the list of edges (specify only one file or raw_edges)
#' @param file Path to the file with the edgelist. Make sure to specify "sep" (if different from " ") and "header" (if different from FALSE) to be passed to the read.table function.
#' @param sep To be used with the "file" variable. Defines the field separator character to be used in the read.table() function. Is set with " " by default.
#' @param header To be used with the "file" variable. A logical value indicating whether "file" contains the names of the variables as its first line. Is set as FALSE by default.
#' @param directed A logical value indicating whether the network is directed or not. Is set as FALSE by default.
#' @return A dataframe each column of which is an ISCV of each node of the graph.
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
