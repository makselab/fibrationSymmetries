#' Find balanced coloring using Kamei algorithm
#'
#' This function find a balanced coloring
#'
#' @param raw_edges list of edges
#' @return Node table
#' @export
get.balanced.coloring.Kamei <- function(raw_edges = NA, file = NA, sep = "\t", header = F, directed = F, weighted = F, look.for.no.input.nodes = T) {
  Rcpp::sourceCpp("kamei_coloring.cpp")
  raw_edges[] <- apply(raw_edges, 2, as.character)
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

  if(weighted == T) {
    integer.weights = unique(raw_edges[, 3])
    integer.edges$Weight = as.integer(factor(raw_edges[, 3], levels = integer.weights))
  }

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

  if(weighted == T) {
    integer.nodes$Color <- getBalancedColoring(integer.nodes$Color, integer.nodes$Fixed, as.matrix(integer.edges), directed, weighted, length(integer.weights))
  } else {
    integer.nodes$Color <- getBalancedColoring(integer.nodes$Color, integer.nodes$Fixed, as.matrix(integer.edges), directed, weighted, 0)
  }
  return(dplyr::select(integer.nodes, -c(Id, Fixed)))
}
