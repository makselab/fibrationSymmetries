#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

std::vector <std::vector<int>> calculateVectors(std::vector<int> nodeColors,
		      std::vector< std::vector<int> > &connections,
		      bool directed,
		      bool weighted,
		      int numberOfWeights) {
	int numberOfColors = *max_element(nodeColors.begin(), nodeColors.end()) + 1;
	int numberOfNodes = nodeColors.size();

	/* Create a 2D vector array to store the matrix of input color vectors.
	First dimension is the node id and second is color id.
	Array needs to be of size numberOfNodes x (numberOfColors * numberOfWeights) because
	each node has an input set color vector (see SI VI in "Circuits with broken fibration symmetries
	perform core logic computations in biological networks. Ian Leifer, Flaviano Morone,
	Saulo D. S. Reis, Jose´ S. Andrade Jr., Mariano Sigman, Hernan A. Makse" for the definition
	of the ISCV) that counts how many nodes of each color send inputs of the given weight
	to the given node. */
	std::vector< std::vector<int> > vectors(numberOfNodes);
	for(int i = 0; i < numberOfNodes; i++) {
		vectors[i].resize(numberOfColors * (weighted?numberOfWeights:1));
	}

	for(int i = 0; i < connections.size(); i++) {
		if(directed == false) {
			int pos = 0;
			if(numberOfWeights == 0) {
				pos = nodeColors[connections[i][1]];
			} else {
				pos = nodeColors[connections[i][1]] * numberOfWeights + connections[i][2];
			}
			vectors[connections[i][0]][pos]++;
		}
		int pos = 0;
		if(numberOfWeights == 0) {
			pos = nodeColors[connections[i][0]];
		} else {
			pos = nodeColors[connections[i][0]] * numberOfWeights + connections[i][2];
		}
		vectors[connections[i][1]][pos]++;
	}
	return(vectors);
}

// returns new number of colors
int classifyNodes(std::vector< std::vector<int> > vectors, std::vector<int> &nodeColors, bool directed, std::vector <int> fixedColorNodes) {
	// first let`s find how many unique types of vectors are out there
	std::vector< std::vector<int> > vectorTypes;
	vectorTypes.push_back(vectors[0]);

	for(int i = 1; i < nodeColors.size(); i++) {
		bool add = 1;
		for(int j = 0; j < vectorTypes.size(); j++) {
			if(vectors[i] == vectorTypes[j]) {add = 0;}
		}
		if(add == 1) {vectorTypes.push_back(vectors[i]);}
	}

	// now let`s recolor node types for next step
	for(int i = 0; i < nodeColors.size(); i++) {
		for(int j = 0; j < vectorTypes.size(); j++) {
			if(vectors[i] == vectorTypes[j]) {nodeColors[i] = j + fixedColorNodes.size();}
		}
	}

	if(directed != 0) {
		for(int i = 0; i < fixedColorNodes.size(); i++) {
			nodeColors[fixedColorNodes[i]] = i;
		}
	}
	return(vectorTypes.size() + fixedColorNodes.size());
}

// [[Rcpp::export]]
std::vector<int> getBalancedColoring(std::vector<int> nodeColors, std::vector<bool> colorFixed, IntegerMatrix edges, bool directed, bool weighted, int numberOfWeights) {
	// Safety
	if(weighted == false) {numberOfWeights = 0;}

	int numberOfColors = *max_element(nodeColors.begin(), nodeColors.end());
	int numberOfNodes = nodeColors.size();
	for(int i = 0; i < nodeColors.size(); i++) {
		nodeColors[i] = nodeColors[i] - 1;
	}

	// Fill in nodes with the fixed color
	std::vector <int> fixedColorNodes;
	for(int i = 0; i < colorFixed.size(); i++) {
		if(colorFixed[i] == true) {
			fixedColorNodes.push_back(i);
		}
	}

	// create connections (that is edges)
	std::vector <std::vector<int>> connections(edges.nrow());
	for(int i = 0; i < connections.size(); i++) {
		connections[i].resize(weighted?3:2);
		for(int j = 0; j < (weighted?3:2); j++) {
			connections[i][j] = edges(i, j) - 1;
		}
	}

	while(1) {
		std::vector< std::vector<int> > vectors(numberOfNodes);
		vectors = calculateVectors(nodeColors, connections, directed, weighted, numberOfWeights);

		int nOC = classifyNodes(vectors, nodeColors, directed, fixedColorNodes);

		if(nOC == numberOfColors) {break;}
		else {numberOfColors = nOC;}
	}
	for(int i = 0; i < nodeColors.size(); i++) {
		nodeColors[i]++;
	}
	return(nodeColors);
}
