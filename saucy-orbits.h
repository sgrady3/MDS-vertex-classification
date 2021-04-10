#include "Graph.h"	//my graph object

#ifdef __cplusplus
extern "C"
{
#endif

#include "saucy.h"

#ifdef __cplusplus
}
#endif


static int init_fixadj1(int, int*);

static void init_fixadj2(int, int, int*);

int compare (const void*, const void*);

int on_automorphism(int, const int*, int, int*, void*);

std::vector<int> get_automorphisms(int, int, std::vector<boost::dynamic_bitset<> >&);

void DFS(int, std::vector<bool>&,std::vector<int>&,std::vector<std::vector<int> >&);

std::vector<std::vector<int> > get_orbits(std::vector<int>&, int);

std::vector<std::vector<int> > find_auto_get_orbits(Graph& , std::vector<boost::dynamic_bitset<> >&);


