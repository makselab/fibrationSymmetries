#' @importFrom magrittr "%>%"

isOnlyMainFiber <- function(block, fiberId) {
  regulatorIds <- block[block$FiberId != fiberId, "FiberId"]
  return(anyDuplicated(regulatorIds) == 0)
}

isFiberSendingToRegulators <- function(edges) {
  feedback <- edges %>%
    dplyr::filter(SourceType == "Fiber") %>%
    dplyr::filter(TargetType == "Regulator")
  return(nrow(feedback) != 0)
}

areAllNodesFromBlockInFiber <- function(block) {
  return(nrow(block[block$NodeType == "Regulator", ]) == 0)
}

isSizeOfInputSetOne <- function(block, edges) {
  edges <- edges[!duplicated(edges[, 1:2]), ]
  fiberNode <- block %>%
    dplyr::filter(NodeType == "Fiber") %>%
    dplyr::select(Node) %>%
    dplyr::summarise(Node = dplyr::first(Node))
  fiberNode <- as.character(fiberNode)
  return(length(edges[edges$Target == fiberNode, "Source"]) == 1)
}

doFibersSendToFibers <- function(edges){
  fiberFiberInteractions <- edges %>%
    dplyr::filter(SourceType == "Fiber") %>%
    dplyr::filter(TargetType == "Fiber")
  return(nrow(fiberFiberInteractions) > 0)
}

fiberHasOneInput <- function(edges) {
  NumberOfInputs <- edges %>%
    dplyr::filter(TargetType == "Fiber") %>%
    dplyr::group_by(Target) %>%
    dplyr::summarise(NumberOfInputs = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(NumberOfInputs = dplyr::first(NumberOfInputs))
  return(NumberOfInputs$NumberOfInputs[1] == 1)
}

isOneNodeFromFiberRegulator <- function(edges) {
  numberOfFiberInputs <- edges %>%
    dplyr::filter(SourceType == "Fiber") %>%
    dplyr::filter(TargetType == "Fiber") %>%
    dplyr::group_by(Target) %>%
    dplyr::summarise(NumberOfInputs = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(NumberOfInputs = dplyr::first(NumberOfInputs))
  return(numberOfFiberInputs$NumberOfInputs == 1)
}

chainCondition <- function(block, edges) {
  # chain condition: only one node has no output, only one node has 2 outputs, all other nodes have 1 input and 1 output
  for(i in 1:nrow(block)) {
    numberOfOutputs <- edges %>%
      dplyr::filter(Source == block$Node[i]) %>%
      dplyr::summarise(NumberOfOutputs = dplyr::n())

    numberOfInputs <- edges %>%
      dplyr::filter(Target == block$Node[i]) %>%
      dplyr::summarise(NumberOfInputs = dplyr::n())

    block$NumberOfOutputs[i] <- numberOfOutputs$NumberOfOutputs[1]
    block$NumberOfInputs[i] <- numberOfInputs$NumberOfInputs[1]
  }
  loopbackNode <- edges %>%
    dplyr::filter(Source == Target)
  loopbackNode <- loopbackNode$Source[1]
  if(is.na(loopbackNode)) {return(F)}
  
  loopTwoOutputs <- block[block$Node == loopbackNode, "NumberOfOutputs"] == 2
  oneNoOutputNode <- nrow(block[block$NumberOfOutputs == 0, ]) == 1
  otherNodesOneOne <- (nrow(block[(block$NumberOfOutputs != 0) & (block$NumberOfOutputs != 2), ]) == nrow(block) - 2)

  return(loopTwoOutputs & oneNoOutputNode & otherNodesOneOne)
}

allSameInput <- function(edges) {
  sources <- edges %>%
    dplyr::group_by(Source) %>%
    dplyr::summarise(N = dplyr::n())
  return(nrow(sources) == 1)
}

classify.block <- function(block, edges, fiberId) {
  Class = ""
  BlockName = ""
  nl = ""
  # this big structure of ifs is hard to understand, but it is drawn in block diagram in file blockdiagram.xml
  if(isFiberSendingToRegulators(edges)) {
    if(doFibersSendToFibers(edges)) {
      Class <- "Feedback Fiber"
      BlockName <- "Feedback Fiber"
      nl <- "Fibonacci"
    } else {
      Class <- "Fibonacci n = 1"
      BlockName <- "Fibonacci n = 1"
      nl <- "n = 1, l = 0"
    }
  } else {
    if(!isOnlyMainFiber(block, fiberId)) {
      Class <- "Multi-layered Fiber"
      BlockName <- "Multi-layered Fiber"
      nl <- "Multi-layered Fiber"
    } else {
      if(areAllNodesFromBlockInFiber(block)) {
        if(isSizeOfInputSetOne(block, edges)) {
          if(chainCondition(block, edges)) {
            Class <- "Chain"
            K <- nrow(block)
            BlockName <- paste(K, "CF", sep = "-")
            nl <- "n = 1, l = 0"
          } else {
            if(allSameInput(edges)) {
              Class <- "Synchronized Star Fiber"
              K <- nrow(block) - 1
              BlockName <- paste(K, "SSF", sep = "-")
              nl <- "n = 1, l = 0"
            } else {
              Class <- "Chain-Star"
              BlockName <- "Chain-Star"
              nl <- "n = 1, l = 0"
            }
          }
        } else {
          Class <- "n > 1"
          BlockName <- "n > 1"
          nl <- "n > 1"
        }
      } else {
        if(doFibersSendToFibers(edges)) {
          if(isOneNodeFromFiberRegulator(edges)) {
            L <- nrow(block[block$NodeType == "Regulator", ])
            K <- nrow(block[block$NodeType == "Fiber", ])
            Class <- "Feed-Forward Fiber"
            BlockName <- paste(L, K, "FFF", sep = "-")
            nl <- paste("n = 1, l = ", L, sep = "")
          } else {
            Class <- "Unclassified"
            BlockName <- "Unclassified"
            nl <- "Unclassified"
          }
        } else {
          if(fiberHasOneInput(edges)) {
            Class <- "Unsynchronized Star Fiber"
            K <- nrow(block) - 1
            BlockName <- paste(K, "USF", sep = "-")
            nl <- "n = 0, l = 1"
          } else {
            L <- nrow(block[block$NodeType == "Regulator", ])
            K <- nrow(block[block$NodeType == "Fiber", ])
            Class <- "FAN Fiber"
            BlockName <- paste(L, K, "FAN", sep = "-")
            nl <- paste("n = 0, l = ", L, sep = "")
          }
        }
      }
    }
  }
  return(c(Class, BlockName, nl))
}

#' Get fiber building blocks classification and/or figures and raw data.
#'
#' get.building.blocks() is used to get fiber building block classification and/or figures and raw data.
#' Returns the list of fiber building blocks as defined in (1) that correspond to the network specified using "raw_edges", "file", "sep" and "header" variables.
#' Classification is returned by default. Classification is obtained following the block diagram specified in "blockClassification.xml" file in the package folder.
#' To get additional output specify "outputFolder" variable and "csv", "png" or "pdf" variables.
#' (1) Morone F, Leifer I, Makse HA. Fibration symmetries uncover the building blocks of biological networks. Proc Natl Acad Sci USA. 2020;117(15):83068314.
#'
#' @param raw_edges 2 or 3 column data frame specifying the list of edges (specify only one file or raw_edges)
#' @param file Path to the file with the edgelist. Make sure to specify "sep" (if different from " ") and "header" (if different from FALSE) to be passed to the read.table function.
#' @param sep To be used with the "file" variable. Defines the field separator character to be used in the read.table() function. Is set with " " by default.
#' @param header To be used with the "file" variable. A logical value indicating whether "file" contains the names of the variables as its first line. Is set as FALSE by default.
#' @param outputFolder Path to output csv, png or pdf files in. Will be created if doesn't exist.
#' @param csv A logical value indicating whether building blocks node and edge csv files need to be written to the outputFolder.
#' @param png A logical value indicating whether illustrations of the building blocks need to be outputted as png files in the outputFolder.
#' @param pdf A logical value indicating whether the summary of the building blocks including figures, classification and node names needs to be put in the pdf file "structures.pdf" in the output folder.
#' @param progressBar A logical value indicating whether to show the progress bar.
#' @return A dataframe that includes:
#' FiberId (the color id corresponding to the output of the get.balanced.coloring.Kamei() function),
#' Nodes (list of nodes in the fiber building block),
#' Fiber (list of nodes in the fiber),
#' Regulators (list of regulator nodes),
#' Class (Classification corresponding to the block name),
#' BlockName (Block name),
#' nl (|n,l> classification of the block).
#' @export
get.building.blocks <- function(raw_edges = NA, file = NA, sep = " ", header = F, outputFolder = NA, csv = F, png = F, pdf = F, progressBar = T) {
  if((csv != F | png != F | pdf != F) & is.na(outputFolder)) {
    stop("Specify outputFolder to get csv, png or pdf files")
  }
  if(!is.na(outputFolder)) {
    if(!dir.exists(outputFolder)) {
      dir.create(outputFolder)
    }
  }
  balancedColoring = get.balanced.coloring.Kamei(raw_edges = raw_edges, file = file, sep = sep, header = header, directed = T)
  raw_edges = get.raw.edges(raw_edges = raw_edges, file = file, sep = sep, header = header)
  weighted = as.logical(ncol(raw_edges) - 2)

  nonTrivialColors = as.matrix(
    balancedColoring %>%
      dplyr::group_by(Color) %>%
      dplyr::summarise(FiberSize = dplyr::n()) %>%
      dplyr::filter(FiberSize > 1) %>%
      dplyr::select(Color)
  )

  blocks = as.data.frame(matrix(rep("", 7 * length(nonTrivialColors)), ncol = 7), stringsAsFactors = F)
  colnames(blocks) = c("FiberId", "Nodes", "Fiber", "Regulators", "Class", "BlockName", "nl")

  graph = igraph::graph_from_edgelist(as.matrix(raw_edges[, 1:2]))
  graph = igraph::set_edge_attr(graph, "Weight", value = raw_edges$V3)

  if(weighted & (png | pdf)) {
    raw_edges$Color = dplyr::group_indices(raw_edges, V3)
    numberOfColors = max(raw_edges$Color)
    if(numberOfColors < 9 & numberOfColors > 2) {
      edgeColors = RColorBrewer::brewer.pal(numberOfColors, "Set1")
    } else {
      edgeColors = grDevices::rainbow(numberOfColors)
    }
    raw_edges$Color = edgeColors[raw_edges$Color]

    graph = igraph::set_edge_attr(graph, "Color", value = raw_edges$Color)

    legendColors = raw_edges[!duplicated(raw_edges$Color), ]
    legendColors$Color = factor(legendColors$Color, levels = edgeColors)
    legendColors = dplyr::arrange(legendColors, Color)$V3
  }

  if(pdf == T) {
    cat(" ", sep = " ", file = paste0(outputFolder, "/structures.Rmd"))
  }

  start.time = Sys.time()
  for(i in 1:nrow(blocks)) {
    if(progressBar == T) {
      progress.bar(i, max = nrow(blocks), startTime = start.time, nowTime = Sys.time())
    }
    fiberId = nonTrivialColors[i]
    fiberGenes = balancedColoring$Name[balancedColoring$Color == fiberId]

    ####################################
    #### Find fibers and regulators ####
    ####################################
    blockGenes = fiberGenes
    colorsToAdd = fiberId

    while(length(colorsToAdd) > 0) {
      currentColor = colorsToAdd[1]
      colorsToAdd = colorsToAdd[-1]

      newRegulatorGenes = unique(raw_edges$V1[raw_edges$V2 %in% dplyr::filter(balancedColoring, Color == currentColor)$Name])
      newRegulatorGenes = newRegulatorGenes[!newRegulatorGenes %in% blockGenes]
      blockGenes = c(blockGenes, newRegulatorGenes)
      if(length(newRegulatorGenes) == 0) {next}

      newColors = balancedColoring$Color[balancedColoring$Name %in% newRegulatorGenes]
      newColors = newColors[duplicated(newColors)]

      if(length(newColors) == 0) {next}

      colorsToAdd = c(colorsToAdd, newColors)
    }

    blocks$FiberId[i] = fiberId
    blocks$Nodes[i] = paste(blockGenes, collapse = ", ")
    blocks$Fiber[i] = paste(fiberGenes, collapse = ", ")
    blocks$Regulators[i] = paste(blockGenes[!blockGenes %in% fiberGenes], collapse = ", ")
    ####################################
    ## Fibers and regulators are found #
    ####################################

    ####################################
    ##### Classify building blocks #####
    ####################################
    block = data.frame(Node = fiberGenes, FiberId = fiberId, NodeType = "Fiber", stringsAsFactors = F)
    if(length(blockGenes) != length(fiberGenes)) {
      block = rbind(block,
                    data.frame(Node = blockGenes[!blockGenes %in% fiberGenes], FiberId = -1, NodeType = "Regulator", stringsAsFactors = F)
      )
    }
    for(j in 1:nrow(block)) {
      if(block$FiberId[j] != -1) {next}
      block$FiberId[j] = balancedColoring$Color[balancedColoring$Name == block$Node[j]]
    }

    subgraph = igraph::induced_subgraph(graph, blockGenes, impl = "auto")

    blockEdges = as.data.frame(igraph::as_edgelist(subgraph), stringsAsFactors = F)
    colnames(blockEdges) = c("Source", "Target")
    if(weighted) {
      blockEdges$Weight = igraph::E(subgraph)$Weight
    }

    for(j in 1:nrow(block)) {
      blockEdges$SourceType[blockEdges$Source == block$Node[j]] = block$NodeType[j]
      blockEdges$TargetType[blockEdges$Target == block$Node[j]] = block$NodeType[j]
    }

    blockClassification = classify.block(block = block, edges = blockEdges, fiberId)

    blocks$Class[i] = blockClassification[1]
    blocks$BlockName[i] = blockClassification[2]
    blocks$nl[i] = blockClassification[3]
    ####################################
    #### Building blocks classified ####
    ####################################
    if(csv) {
      outputBlock = block
      colnames(outputBlock)[1] = "Label"
      outputBlock$Id = outputBlock$Label
      outputBlock <- outputBlock[, c("Id", "Label", "FiberId")]
      write.csv(outputBlock, file = paste0(outputFolder, "/", fiberId, "_nodes.csv"), quote = F, row.names = F)
      outputEdges = blockEdges
      colnames(outputEdges)[4:5] = c("Type", "Color")
      outputEdges$Type = "directed"
      outputEdges$Color = igraph::E(subgraph)$Color
      write.csv(outputEdges, file = paste0(outputFolder, "/", fiberId, "_edges.csv"), quote = F, row.names = F)
    }

    if(png | pdf) {
      block$Node = factor(block$Node, igraph::as_data_frame(subgraph, what = "vertices")$name)
      block = dplyr::arrange(block, Node)

      igraph::V(subgraph)$color = dplyr::group_indices(block, FiberId)

      png(filename = paste0(outputFolder, "/", fiberId, ".png"), width = 640, height = 360)
      oldMargins <- par("mar")
      par(mar = c(0, 0, 0, 0))
      if(!weighted) {
        plot(subgraph)
      } else {
        plot(subgraph, edge.color = igraph::E(subgraph)$Color)
        legend(x = 1.2, y = 1.1, legend = legendColors,
               col = edgeColors, lty = 1, lwd = 3, cex = 1,
               text.font = 4, bg = 'white')
      }
      par(mar = oldMargins)
      dev.off()

      if(pdf == T) {
        if(grepl("n =", blocks$nl[i])) {
          cat(paste0("|", blocks$nl[i], "> \\"), sep = "\n",
              file = paste0(outputFolder, "/structures.Rmd"),
              append = T)
        } else {
          cat(paste0(blocks$nl[i], " \\"), sep = "\n",
              file = paste0(outputFolder, "/structures.Rmd"),
              append = T)
        }
        cat(paste0("Fiber nodes: ", blocks$Fiber[i], " \\"), sep = "\n",
            file = paste0(outputFolder, "/structures.Rmd"),
            append = T)
        cat(paste0("Regulator nodes: ", blocks$Regulators[i], " \\"), sep = "\n",
            file = paste0(outputFolder, "/structures.Rmd"),
            append = T)
        cat(paste0("![Alt text](", outputFolder, "/", fiberId, ".png) \\"), sep = "\n",
            file = paste0(outputFolder, "/structures.Rmd"),
            append = T)
        cat("\\newpage", sep = "\n",
            file = paste0(outputFolder, "/structures.Rmd"),
            append = T)
      }
    }
  }

  if(pdf == T) {
    rmarkdown::render(input = paste0(outputFolder, "/structures.Rmd"), output_format = rmarkdown::pdf_document())
    file.remove(paste0(outputFolder, "/structures.Rmd"))
    if(png == F) {
      for(fiberId in blocks$FiberId)
        file.remove(paste0(outputFolder, "/", fiberId, ".png"))
    }
  }

  return(blocks)
}
