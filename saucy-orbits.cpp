
#include "saucy-orbits.h"

#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

char *marks=0;
int **automorphs=0;

/* return value >1 indicates error */
static int
dupe_check(int n, int *adj, int *edg)
{
	int i, j, self_loop_ctr;
	int *dupe_tmp = (int*)calloc(n, sizeof(int));
	if (!dupe_tmp) {
		cout<<"can't allocate memory"<<"\n";
		free(dupe_tmp);
		return 2;
	}

	/* Check outgoing edges of each vertex for duplicate endpoints */
	for (i = 0; i < n; ++i) {
		self_loop_ctr = 0;
		for (j = adj[i] ; j < adj[i+1] ; j++) {
			/* Self-loops lead to two entries of edg[j]==i,
			 * so we check for those and only worry if we see
			 * 3 hits of edg[j]==i (which means 2 self-loops).
			 */
			if (edg[j] == i) {
				++self_loop_ctr;
				if (self_loop_ctr > 2) {
					cout<<"self loop duplicate\n";
					cout<<"duplicate edge in input\n";
					free(dupe_tmp);
					return 1;
				}
			}
			/* If we have recorded this vertex as connected to i,
			 * we have a dupe.
			 * Using i+1 because we used calloc above, and 0 is
			 * a valid vertex index.
			 */
			else if (dupe_tmp[edg[j]] == i+1) {
				cout<<"The other duplicate one\n";
				cout<<"duplicate edge in input\n";
				free(dupe_tmp);
				return 1;
			}
			dupe_tmp[edg[j]] = i+1;
		}
	}

	free(dupe_tmp);
	return 0;
}

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


std::vector<int> get_automorphisms(int n, int e, std::vector<boost::dynamic_bitset<> > & adjMat)
{
	cout<<"getting automorphisms\n";
	struct saucy_stats stats;
	int *aout=0; 
	int *eout=0;
	int *ain=0; 
	int *ein=0;
	int *colors=0;
	cout<<"allocated first ints\n";
	//int n=adjMat.size();
	//int e = g.numEdges();
	struct saucy *s=0;
	struct saucy_graph sg;
	struct saucy_graph* sg_ptr=&sg;
	cout<<"allocated structs\n";
	size_t index;
	aout = (int*)calloc((n+1), sizeof(int));
	eout = (int*)malloc(2 * e * sizeof(int));
	colors = (int*)malloc(n * sizeof(int));

	s= saucy_alloc(n);
	cout<<"allocated saucy\n";
	sg.n=n;
	cout<<"not the n\n";
	sg.e=e;
	cout<<"not the e\n";
	sg.adj=aout;
	cout<<"not the adj\n";
	sg.edg=eout;
	cout<<"set up the saucy graph\n";
	ain = aout + (0);
	ein = eout + (0);
	

	for (int i = 0; i < n; ++i)
  	{ 
  		colors[i]=0;
	   	index=adjMat[i].find_first();
	   	while(index!=boost::dynamic_bitset<>::npos) 
	   	{
	   		if (index>i)
	   		{
	   			++aout[i];
	   			++ain[index];
	   		}
	   		index=adjMat[i].find_next(index);
	   	}
  	}
  	init_fixadj1(n,aout);

  	for (int i = 0; i < n; ++i)
  	{ 
	   	index=adjMat[i].find_first();
	   	while(index!=boost::dynamic_bitset<>::npos) 
	   	{
	   		if (index>i)
	   		{
	   			eout[aout[i]++] = index;
				ein[ain[index]++] = i;
	   		}
	   		index=adjMat[i].find_next(index);
	   	}
  	}
  	init_fixadj2(n, 2 * e, aout);

  	if (dupe_check(n, aout, eout))
  	{
  		cout<<"there are duplicate edges\n";
  		exit(1);
  	}

  	marks =(char*)calloc(n, sizeof(char));
  	automorphs =(int**)calloc(n,sizeof(int**));
  	for (int i = 0; i < n; ++i)
  	{
  		automorphs[i]=(int*)calloc(n,sizeof(int));
  	}
  	cout<<"doing the saucy search\n";
  	saucy_search(s,sg_ptr,0,colors,on_automorphism,0,&stats);
  	std::vector<int> auto_morphisms;
	//cout<<"auto_morphisms size: "<<auto_morphisms.size()<<"\n";
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if(automorphs[i][j])
			{
				//cout<<i<<" "<<j<<"\n";
				auto_morphisms.push_back(i);
				auto_morphisms.push_back(j);
			}
		}
	}

	saucy_free(s);
	//cout<<"freed saucy\n";
	free(marks);
	//cout<<"freed marks\n";
	free (automorphs);
	//cout<<"freed automorphs\n";
	//cout<<"automorphisms size: "<<auto_morphisms.size()<<"\n";
  	//free(sg);
  	return auto_morphisms;
}


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

std::vector<std::vector<int> > find_auto_get_orbits(Graph& g, std::vector<boost::dynamic_bitset<> > & adjMat)
{
	cout<<"find auto get orbits\n";
	std::vector<int> auto_morphisms=get_automorphisms(g.numVertices(),g.numEdges(),adjMat);
	cout<<"getting orbits\n";
	std::vector<std::vector<int> > orbits=get_orbits(auto_morphisms,g.numVertices());
	return orbits;
}