
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include "Graph.h"

using namespace std;


Graph::Graph(void){
	//cout<<"Creating graph object"<<endl;
	return;
}

void Graph::readWeighted(ifstream& infile){
	int v1,v2;
	double weight;
	string line,temp;
	stringstream ss;
	int ncols=0;

	getline(infile,line);

	getline(infile,line);
	ss.clear();
	ss<<line;
	while(ss>>temp){
		ncols++;
	}

	infile.clear();
	infile.seekg(0,ios::beg);


	infile>>nVertices>>nEdges;
	numActiveVertices=nVertices;
	std::vector<double> bvector(nVertices,0);
	for (int i = 0; i < nVertices; i++)
	{
		weightedMat.push_back(bvector);
		degrees.push_back(0);
		active.push_back(true);
	}


	if(ncols<3)
	{
		while(infile>>v1>>v2){
		weightedMat[v1][v2]=1;
		weightedMat[v2][v1]=1;
		degrees[v1]++;
		degrees[v2]++;
	}
	}
	else{
		while(infile>>v1>>v2>>weight){
		weightedMat[v1][v2]=weight;
		weightedMat[v2][v1]=weight;
		degrees[v1]++;
		degrees[v2]++;
		}
	}

	return;
}

void Graph::readGraph(ifstream& infile){
	string v1,v2;
	int k=0;
	std::map<string, int> labelMap;

	infile>>nVertices>>nEdges;
	numActiveVertices=nVertices;
	adjMat.resize(nVertices);
	degrees.resize(nVertices);
	active.resize(nVertices);
	boost::dynamic_bitset<> bvector(nVertices);
	adjList.resize(nVertices);
	for (int i = 0; i < nVertices; i++)
	{
		//adjMat.push_back(bvector);
		adjMat[i]=bvector;
		degrees[i]=0;
		active[i]=true;
	}

	while(infile>>v1>>v2){

		if (labelMap.find(v1)==labelMap.end())
		{
			labels.push_back(v1);
			labelMap.insert(make_pair(v1,k));
			k++;
		}
		if (labelMap.find(v2)==labelMap.end())
		{
			labels.push_back(v2);
			labelMap.insert(make_pair(v2,k));
			k++;
		}

		if(!adjMat[labelMap[v1]][labelMap[v2]] && labelMap[v1]!=labelMap[v2])
		{
			adjMat[labelMap[v1]][labelMap[v2]]=true;
			adjMat[labelMap[v2]][labelMap[v1]]=true;
			adjList[labelMap[v1]].insert(labelMap[v2]);
			adjList[labelMap[v2]].insert(labelMap[v1]);
			degrees[labelMap[v1]]++;
			degrees[labelMap[v2]]++;
		}
	}

	if (k>nVertices)
	{
		cerr<<"Bad file format. Vertices don\'t match labels.\n";
		exit(1);
	}

	return;
}


double Graph::density(void){
	double V=nVertices;
	double E=nEdges;
	return (E/((V*(V-1))/2.0));
}

int Graph::maxDegree(void){
	int maxDegree=0;
	for (int i = 0; i < nVertices; i++)
	{
		if(degrees[i]>maxDegree)
			maxDegree=degrees[i];
	}
	return maxDegree;
}

int Graph::minDegree(void){
	int minDegree;
	minDegree=nVertices-1;
	for (int i = 0; i < nVertices; i++)
	{

		if(degrees[i]<minDegree)
			minDegree=degrees[i];
	}
	return minDegree;
}

/*std::vector<int> Graph::neighbors(int v){
	std::vector<int> array;
	for (int i = 0; i < nVertices; i++)
	{
		if(edgeExists(v,i))
			array.push_back(i);
	}
	return array;
}*/

boost::dynamic_bitset<> Graph::neighbors(int v){
	return adjMat[v];
}

std::unordered_set<unsigned int> Graph::neighborsList(int v)
{
	return adjList[v];
}

/*std::vector<int> Graph::closedNeighbors(int v){
	std::vector<int> array;
	array.push_back(v);
	for (int i = 0; i < nVertices; i++)
	{
		if(edgeExists(v,i))
			if (active[i])
			array.push_back(i);
	}
	return array;
}*/

boost::dynamic_bitset<> Graph::closedNeighbors(int v){
	boost::dynamic_bitset<> cneigh=neighbors(v);
	cneigh.flip(v);
	return cneigh;
}


std::vector<int> Graph::wNeighbors(int v){
	std::vector<int> array;
	for (int i = 0; i < nVertices; i++)
	{
		if(hasWeight(v,i))
			if (active[i])
			array.push_back(i);
	}
	return array;
}

void Graph::addEdge(int v1,int v2){
	adjMat[v1][v2]=true;
	adjMat[v2][v1]=true;
	degrees[v1]++;
	degrees[v2]++;
	nEdges++;
	return;
}

void Graph::deleteEdge(int v1,int v2){
	adjMat[v1][v2]=false;
	adjMat[v2][v1]=false;
	degrees[v1]--;
	degrees[v2]--;
	nEdges--;
	return;
}

void Graph::deleteVertex(int v){

	boost::dynamic_bitset<> array;
	array=neighbors(v);
	active[v]=false;
	for (boost::dynamic_bitset<>::size_type i = 0; i < array.size(); i++)
	{
		deleteEdge(v,array[i]);
	}

	return;
}

void Graph::deleteLowerDegreeVertices(int k){
	for (int i = 0; i < nVertices; i++)
	{
		if(vertexExists(i))
			if (degree(i)<k)
				deleteVertex(i);
	}
}

void Graph::printGraph(void){
	for (int i = 0; i < nVertices; i++)
	{
		for (boost::dynamic_bitset<>::size_type  j = 0; j < nVertices; j++)
		{
			cout<< adjMat[i][j]<<" ";
		}
		cout<< endl;
	}
	return;
}

void Graph::printWeighted(void){
	for (int i = 0; i < nVertices; i++)
	{
		for (int j = 0; j < nVertices; j++)
		{
			cout<< weightedMat[i][j]<<" ";
		}
		cout<< endl;
	}
	return;
}

Graph::~Graph(void){
	//cout<<"deleting graph"<<endl;
	return;
}
