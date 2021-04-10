#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <fstream>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include <bits/stdc++.h>


class Graph
{

	unsigned int nVertices;
	unsigned int nEdges;
	unsigned int numActiveVertices;
	std::vector<boost::dynamic_bitset<> > adjMat;
	std::vector<std::vector<double> > weightedMat;
	std::vector<int> degrees;
	boost::dynamic_bitset<> active;
	std::vector<std::string> labels;
	std::vector<std::unordered_set<unsigned int> > adjList;

public:
	Graph();
	void readGraph(std::ifstream& );
	void readWeighted(std::ifstream&);
	double density(void);
	int maxDegree(void);
	int minDegree(void);
	boost::dynamic_bitset<> neighbors(int);
	boost::dynamic_bitset<> closedNeighbors(int);
	std::unordered_set<unsigned int> neighborsList(int);
	std::vector<int> wNeighbors(int);
	void addEdge(int,int);
	void deleteEdge(int,int);
	void deleteVertex(int);
	void deleteLowerDegreeVertices(int);
	void printGraph(void);
	void printWeighted(void);
	~Graph();

	inline int numVertices(void){
		return nVertices;
	}

	inline int numEdges(void){
		return nEdges;
	}

	inline bool edgeExists(int v1,int v2){
		return adjMat[v1][v2];
	}

	inline bool isActive(int i){
		return active[i];
	}

	inline void inactive(int i){
		active[i]=false;
	}

	inline bool hasWeight(int v1, int v2){
		if(edgeWeight(v1,v2)!=0)
			return true;
		else
			return false;
	}

	inline double edgeWeight(int v1, int v2){
		return weightedMat[v1][v2];
	}

	inline bool vertexExists(int v){
		return active[v];
	}

	inline int degree(int v){
		return degrees[v];
	}

	inline std::string label(int v){
		return labels[v];
	}
	inline int indexFromLabel(std::string v){
		for (int i = 0; i < labels.size(); ++i)
		{
			if (labels[i].compare(v)==0)
			{
				return i;
			}
		}
	}
	
};

#endif