# Functions to work with fibration symmetries.

## Installation

```r
# Install the package from the GitHub directly in RStudio by running:
devtools::install_github("ianleifer/fibrationSymmetries")
```

## Usage
Package features 4 functions:

1) get.balanced.coloring.Kamei() finds the minimal balanced coloring of the graph using the algorithm introduced by Kamei and Cock.

Example usage:
nodes = get.balanced.coloring.Kamei(file = "Network.txt", directed = T)

nodes 

get.input.set.color.vector()	Find the ISCV (Input Set Color Vector) of all nodes in the graph given the coloring.

get.building.blocks()	Get fiber building blocks classification and/or figures and raw data.

```xml
<myxml>
   <someElement />  
</myxml>
```

get.building.block.pvalues()	Get p-values corresponding to fiber building blocks.

Package can be installed directly in RStudio by running
devtools::install_github("ianleifer/fibrationSymmetries")