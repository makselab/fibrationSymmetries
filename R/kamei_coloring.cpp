#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

std::vector <std::vector<int>> calculateVectors(std::vector<int> nodeColors,
                      std::vector< std::vector<int> > vectors,
		      std::vector< std::vector<int> > &connections,
		      bool directed,
		      bool weighted,
		      int numberOfWeights) {
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
int classifyNodes(std::vector< std::vector<int> > vectors, std::vector<int> &nodeColors, bool directed, std::vector <int> noInputNodes) {
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
			if(vectors[i] == vectorTypes[j]) {nodeColors[i] = j + noInputNodes.size();}
		}
	}

	if(directed != 0) {
		for(int i = 0; i < noInputNodes.size(); i++) {
			nodeColors[noInputNodes[i]] = i;
		}
	}
	return(vectorTypes.size() + noInputNodes.size());
}

// [[Rcpp::export]]
std::vector<int> getBalancedColoring(std::vector<int> nodeColors, std::vector<bool> colorFixed, IntegerMatrix edges, bool directed, bool weighted, int numberOfWeights) {
	if(weighted == false) {numberOfWeights = 0;}
	int numberOfColors = *max_element(nodeColors.begin(), nodeColors.end());
	std::vector <int> noInputNodes;
	for(int i = 0; i < colorFixed.size(); i++) {
		if(colorFixed[i] == true) {
			noInputNodes.push_back(i);
		}
	}

	int numberOfNodes = nodeColors.size();
	
	for(int i = 0; i < nodeColors.size(); i++) {
		nodeColors[i] = nodeColors[i] - 1;
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
		// create 2D vector array to store all vectors belonging to each node
		/* Explanation why array is of size numberOfNodes x (numberOfColors * numberOfWeights). There are two ways how to do it.
		Either it can be done as a 3D array and then we will need two realisations for weighted and non-weighted design.
		Or vectors themselves can be formed in a bit weird way, but we will classify nodes comparing vectors not worrying about their structure.
		It improves readability and simpliness only paying with the strange enumeration of array */
		std::vector< std::vector<int> > vectors(numberOfNodes);
		for(int i = 0; i < numberOfNodes; i++) {
			vectors[i].resize(numberOfColors * (weighted?numberOfWeights:1));
		}

		vectors = calculateVectors(nodeColors, vectors, connections, directed, weighted, numberOfWeights);
		int nOC = classifyNodes(vectors, nodeColors, directed, noInputNodes);

		if(nOC == numberOfColors) {break;}
		else {numberOfColors = nOC;}
	}
	for(int i = 0; i < nodeColors.size(); i++) {
		nodeColors[i]++;
	}
	return(nodeColors);
}
