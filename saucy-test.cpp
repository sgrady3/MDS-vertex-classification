//g++ -o saucy-test saucy-test.cpp Graph.cpp -I /home/stephen/code/saucy-3.0/
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "Graph.h"	//my graph object

#ifdef __cplusplus
extern "C"
{
#endif

#include "saucy.h"

#ifdef __cplusplus
}
#endif



using namespace std;

char *marks;	//global variable to mark a vertex when found in automorphism
int **automorphs;	//global variable to mark when two vertices share an automorphism

//saucyio.c function to fix adjacencies
static int
init_fixadj1(int n, int *adj)
{
	int val, sum, i;

	/* Translate adj values to real locations */
	val = adj[0]; sum = 0; adj[0] = 0;
	for (i = 1; i < n; ++i) {
		sum += val;
		val = adj[i];
		adj[i] = sum;
	}
	return sum + val;
}

//saucyio.c function to fix edges
static void
init_fixadj2(int n, int e, int *adj)
{
	int i;

	/* Translate again-broken sizes to adj values */
	for (i = n-1; i > 0; --i) {
		adj[i] = adj[i-1];
	}
	adj[0] = 0;
	adj[n] = e;
}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


//function to store or print automorphism when found. Copied and modified from on_automorphism in saucyio.c
int on_automorphism(int n, const int *gamma, int k, int *support, void *arg)
{


	qsort(support,k,sizeof(int),compare);
	int standin_k, j;
	for (int i = 0; i < k; ++i)
	{
		standin_k=support[i];
		if (marks[standin_k]) continue;
		marks[standin_k]=1;
		//cout<<"("<<standin_k;
		for (j = gamma[standin_k]; j != standin_k; j = gamma[j]) {
			marks[j] = 1;
			//cout<<" "<<j;
			automorphs[standin_k][j]=1;
		}
		//cout<<")";
	}
	//cout<<"\n";
	for (int i = 0; i < k; ++i)
	{
		marks[support[i]]=0;
	}

	return 1;
}

// main function to interface with saucy
std::vector<int> get_automorphisms(int n, int e, std::vector<boost::dynamic_bitset<> > & adjMat)
{
	cout<<"getting automorphisms\n";
	struct saucy_stats stats; //internal saucy stats
	int *aout, *eout, *ain, *ein, *colors; //pointers to saucy 
	//int n=adjMat.size();
	//int e = g.numEdges();
	struct saucy *s;	//internal saucy struct
	struct saucy_graph sg; //internal saucy graph struct
	struct saucy_graph* sg_ptr=&sg; //pointer to saucy graph
	size_t index;
	aout = (int*)calloc((n+1), sizeof(int));   //allocate memory for saucy adjacency
	eout = (int*)malloc(2 * e * sizeof(int));  //allocate memory for saucy edes
	colors = (int*)malloc(n * sizeof(int));    //allocate memory for saucy colors (only one color in our case for each vertex)
	s= saucy_alloc(n);	//internal suacy function to allocate saucy memory
	sg->n=n;            //set saucy graph number of vertices to n
	sg->e=e;    		//set saucy graph number of edges to e
	sg->adj=aout;		//set saucy adjacency to allocated aout
	sg->edg=eout;		//set saucy edges to eout

	ain = aout + (0);	//used for digraphs
	ein = eout + (0);	//used for digraphs
	

	for (int i = 0; i < n; ++i)
  	{ 
  		colors[i]=0;							//set color for vertex i to 0
	   	index=adjMat[i].find_first();
	   	while(index!=boost::dynamic_bitset<>::npos) 
	   	{
	   		if (index>i)
	   		{
	   			++aout[i];					//add adjacency for i
	   			++ain[index];				//add adjacency for index
	   		}
	   		index=adjMat[i].find_next(index);
	   	}
  	}
  	init_fixadj1(n,aout);					//function copied from saucyio.c to fix adjacencies

  	for (int i = 0; i < n; ++i)
  	{ 
	   	index=adjMat[i].find_first();
	   	while(index!=boost::dynamic_bitset<>::npos) 
	   	{
	   		if (index>i)
	   		{
	   			eout[aout[i]++] = index;	//set edges between i and index
				ein[ain[index]++] = i;		//set edges between index and i
	   		}
	   		index=adjMat[i].find_next(index);
	   	}
  	}
  	init_fixadj2(n, 2 * e, aout);			//function copied from saucyio.c to fix edges

  	marks =(char*)calloc(n, sizeof(char));	//global variable function to mark if vertex has been added to automorphism when we find an automrphism. Used in on_automorphism
  	automorphs =(int**)calloc(n,sizeof(int**));	//global variable bit array used to store if two vertices are in an automorphism
  	for (int i = 0; i < n; ++i)
  	{
  		automorphs[i]=(int*)calloc(n,sizeof(int));	//allocate memory for automorphs bit array
  	}

  	saucy_search(s,sg,0,colors,on_automorphism,0,&stats); //internal saucy function to generate automorphisms
  	std::vector<int> auto_morphisms;
	//cout<<"auto_morphisms size: "<<auto_morphisms.size()<<"\n";
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if(automorphs[i][j])
			{
				//cout<<i<<" "<<j<<"\n";
				auto_morphisms.push_back(i);	//create vector of automorphimisms where every i%2==0 and i+1 are in an automorphhism
				auto_morphisms.push_back(j);
			}
		}
	}

	saucy_free(s);	//free up saucy
	//cout<<"freed saucy\n";
	free(marks);	//free up marks
	//cout<<"freed marks\n";	
	free (automorphs);	//free up automorphs bitarray
	//cout<<"freed automorphs\n";
	//cout<<"automorphisms size: "<<auto_morphisms.size()<<"\n";
  	//free(sg);
  	return auto_morphisms;
}


// depth first search recursive fucntion
void DFS(int v, std::vector<bool>& visited,std::vector<int>& orbit,std::vector<std::vector<int> >&adjList)
{
	visited[v]=true;
	orbit.push_back(v);

	for (int i = 0; i < adjList[v].size(); ++i)
	{
		if (!visited[adjList[v][i]])
		{
			DFS(adjList[v][i],visited,orbit,adjList);
		}
	}
}


//use connected components to put together automorphsims into orbits
std::vector<std::vector<int> > get_orbits(std::vector<int>& auto_morphisms, int n)
{
	std::vector<int> list;
	std::vector<int> orbit;
	std::vector<std::vector<int> > orbits;
	std::vector<std::vector<int> > adjList;
	std::vector<bool> visited(n);
	for (int i = 0; i < n; ++i)
	{
		adjList.push_back(list);
	}

	for (int i = 0; i < auto_morphisms.size(); i=i+2)
	{
		adjList[auto_morphisms[i]].push_back(auto_morphisms[i+1]);
		adjList[auto_morphisms[i+1]].push_back(auto_morphisms[i]);
	}

	//cout<<"doing the dfs\n";
	for (int i = 0; i < n; ++i)
	{
		if (!visited[i])
		{
			orbit.clear();
			DFS(i,visited,orbit,adjList);
			orbits.push_back(orbit);
		}
	}
	/*
	for (int i = 0; i < orbits.size(); ++i)
	{
		for (int j = 0; j < orbits[i].size(); ++j)
		{
			cout<<orbits[i][j]<<" ";
		}
		cout<<"\n";
	}
	*/

	return orbits;
}

//get function to get automorphisms then to get orbits from automorphisms
std::vector<std::vector<int> > find_auto_get_orbits(Graph& g, std::vector<boost::dynamic_bitset<> > & adjMat)
{
	//cout<<"getting automorphisms\n";
	std::vector<int> auto_morphisms=get_automorphisms(g.numVertices(),g.numEdges(),adjMat);
	cout<<"getting orbits\n";
	std::vector<std::vector<int> > orbits=get_orbits(auto_morphisms,g.numVertices());
	return orbits;
}


int main(int argc, char const *argv[])
{
	Graph g; //my graph object
	ifstream graphFile(argv[1]); //get graph file from command line
	g.readGraph(graphFile); //read in graph from file
	graphFile.close(); //close graph file
	const int num=g.numVertices(); //get number of vertices
	std::vector<boost::dynamic_bitset<> > adjMat(num); //set up bit adjacency matrix
	for (int i = 0; i < num; ++i)
	{
		adjMat[i]=g.closedNeighbors(i);
	}
	//cout<<"number of vertices: "<<g.numVertices()<<"\n";
	//cout<<"getting automorphisms\n";
	//std::vector<int> auto_morphisms=get_automorphisms(g,adjMat);
	//cout<<"got automorphisms\n";
	//cout<<"number of vertices: "<<g.numVertices()<<"\n";
	//cout<<"returned the automorphisms\n";
	//cout<<"automorphisms size: "<<auto_morphisms.size()<<"\n";
	/*
	for (int i = 0; i < auto_morphisms.size(); i=i+2)
	{
		//cout<<"let's try this\n";
		cout<<auto_morphisms[i]<<" "<<auto_morphisms[i+1]<<"\n";
	}
	*/
	//cout<<"getting orbits\n";
	//std::vector<std::vector<int> > orbits=get_orbits(auto_morphisms,num);
	std::vector<std::vector<int> > orbits=find_auto_get_orbits(g,adjMat); //get orbits from saucy
	/*
	for (int i = 0; i < orbits.size(); ++i)
	{
		if (orbits[i].size()>1)
		{
			for (int j = 0; j < orbits[i].size(); ++j)
			{
				cout<<orbits[i][j]<<" ";
			}
			cout<<"\n";
		}
	}
	*/
	cout<<"got orbits!\n";
	return 0;
}