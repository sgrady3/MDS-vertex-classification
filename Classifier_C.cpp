// g++ -m64 -std=c++11 -o essential-DS Graph.cpp essential-DS.cpp computeOrbits.cpp /home/stephen/code/nauty26r11/nauty.a -I /home/stephen/code/nauty26r11 -I /home/stephen/code/gurobi752/linux64/include/ -L /home/stephen/code/gurobi752/linux64/lib/ -lgurobi_c++ -lgurobi75 -lm


#include <stdlib.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <chrono>
#include "gurobi_c++.h"
#include "Graph.h"
//#include "computeOrbits.h"
#include "saucy-orbits.h"


using namespace std;

struct Classifications
{
   int MDSsize=0;
   int byOrbit=0;
   int byTrueOrbit=0;
   int byessentialrule=0;
   int bySubset=0;
   int intswap=0;
   int byredundantrule=0;
   int MDStest=0;
   boost::dynamic_bitset<> MDS;                         //bitset for the initial MDS
   boost::dynamic_bitset<> essential;                   //bitset of essential vertices
   boost::dynamic_bitset<> intermittent;                //bitset of intermittent vertices
   boost::dynamic_bitset<> redundant;                   //bitset of redundant vertices
   boost::dynamic_bitset<> done;                        //bitset of vertices that have been classified
   boost::dynamic_bitset<> orbitComputed;               //bitset
   std::vector<std::vector<int> > orbits;
   std::vector<std::vector<int> > propersubset;
   std::vector<std::vector<int> > propersuperset;
   std::vector<std::vector<int> > subset;
   std::vector<GRBVar> variables;

};

bool maxSort(int i, int j) {return (i>j);}

void printVec(std::vector<int>& v)
{
    const int numV=v.size();
    for (int i = 0; i < numV; ++i)
    {
        cout<<v[i]<<" ";
    }
    cout<<"\n";
}

void printBitsSet(boost::dynamic_bitset<>& v)
{
    const int numV=v.size();
    for (int i = 0; i < numV; ++i)
    {
        if (v[i])
        {
            cout<<i<<" ";
        }
    }
    cout<<"\n";
}


//determine orbits with call to nauty
std::vector<std::vector<int> > createOrbits(std::vector<boost::dynamic_bitset<> >& bitsets, int num_edges)
{
    std::vector<std::vector<int> > orbitSets;
    //std::vector<std::vector<int> > denseorbitSets;
    cout<<"starting orbits\n";
    //orbitSets=sparse_orbits(bitsets,num_edges);
    //denseorbitSets=orbits(bitsets);

    int numV=bitsets.size();
    for (int i = 0; i < numV; ++i)
    {
        if (orbitSets[i].size()>1)
        {
            for (int j = 1; j < orbitSets[i].size(); ++j)
            {
                orbitSets[orbitSets[i][j]]=orbitSets[i];
            }
        }
    }
    cout<<orbitSets.size()<<"\n";
    return orbitSets;
}

//determine orbits with saucy
std::vector<std::vector<int> > createOrbits_saucy(Graph& g, std::vector<boost::dynamic_bitset<> >& bitsets)
{
    std::vector<std::vector<int> > orbitSets;
    orbitSets=find_auto_get_orbits(g,bitsets);
    std::vector<std::vector<int> > orbits(g.numVertices());
    for (int i = 0; i < orbitSets.size(); ++i)
    {
        for (int j = 0; j < orbitSets[i].size(); ++j)
        {
            //cout<<nauty_orbits[i][j]<<" ";
            orbits[orbitSets[i][j]]=orbitSets[i];
        }
        //cout<<"\n";
    }
    return orbits;
}

//create a Gurobi model
void createModel(GRBModel& model, std::vector<boost::dynamic_bitset<> >& bitsets, Classifications& classes, Graph& g)
{

    int numV=bitsets.size();
    for (int i = 0; i < numV; ++i)
    {
        //classes.subset[i].resize(numV);
        bitsets[i]=g.closedNeighbors(i);
        classes.variables[i]=model.addVar(0,1,0,GRB_BINARY,to_string(i));
    }

    GRBLinExpr expr, objexpr;
    objexpr=0;
    size_t index,index2;
    for (int i = 0; i < numV; ++i)
    {
        objexpr+=classes.variables[i];
        expr=0;
        index=bitsets[i].find_first();
        while(index!=boost::dynamic_bitset<>::npos) 
        {
            if (i!=index && bitsets[index].is_subset_of(bitsets[i]))
            {
                //cout<<index<<" is subset of "<<i<<"\n";
                classes.subset[i].push_back(index);
            }
            if (i!=index && bitsets[index].is_proper_subset_of(bitsets[i]))
            {
                //cout<<index<<" is proper subset of "<<i<<"\n";
                classes.propersubset[i].push_back(index);
                classes.propersuperset[index].push_back(i);
            }
            expr+=classes.variables[index];
            index=bitsets[i].find_next(index);
        }
        model.addConstr(expr>=1,to_string(i));

    }
    model.setObjective(objexpr,GRB_MINIMIZE);
    model.set(GRB_IntParam_OutputFlag,0);
    return;
}

//classify all unclassified vertices that share an orbit with v
int classifyOrbit(Classifications& classes, int v, int classification)
{
    int num_classified=0;
    //cout<<"classifying "<<v<<" in orbits code as "<<classification<<"\n";
    //cout<<classes.orbits.size()<<"\n";
    if (classes.orbits[v].size()>1)
    {
        //cout<<"the first if statement checked out\n";
        //cout<<"orbits size "<<classes.orbits[v].size()<<"\n";
        for (int i = 0; i < classes.orbits[v].size(); ++i)
        {
            //cout<<i<<" "<<classes.orbits[v][i]<<"\n";
            if (!classes.done[classes.orbits[v][i]])
            {
                if(classification==2)
                {
                    classes.essential[classes.orbits[v][i]]=true;
                    classes.done[classes.orbits[v][i]]=true;
                    classes.orbitComputed[classes.orbits[v][i]]=true; 
                    classes.byOrbit++;
                    num_classified++;
                }
                else if (classification==1)
                {
                    classes.intermittent[classes.orbits[v][i]]=true;
                    classes.done[classes.orbits[v][i]]=true;
                    classes.orbitComputed[classes.orbits[v][i]]=true;
                    classes.byOrbit++;
                    num_classified++;
                }
                else if (classification==0)
                {
                    //cout<<"setting redundant\n";
                    classes.redundant[classes.orbits[v][i]]=true;
                    //cout<<"setting done\n";
                    classes.done[classes.orbits[v][i]]=true;
                    //cout<<"setting orbitComputed\n";
                    classes.orbitComputed[classes.orbits[v][i]]=true;
                    classes.byOrbit++;
                    num_classified++;
                    //cout<<"done\n";
                }
                else{
                    cout<<"error in classifyOrbit, improper classification given\n";
                }
            }
        }           
    }
    return num_classified;
}


//Rule 1 in paper essential subset rule if vertex is subset of essential vertex then it is redunant
void essentialSubset(Classifications& classes, int v)
{
    int vertex;
    for (int i = 0; i < classes.propersubset[v].size(); ++i)
    {
        vertex=classes.propersubset[v][i];
        if(!classes.done[vertex])
        {
            classes.redundant[vertex]=true;
            classes.bySubset++;
            classes.done[vertex]=true;
            classes.bySubset+=classifyOrbit(classes,vertex,0);
            classes.variables[vertex].set(GRB_DoubleAttr_UB,0.0);
        }
    }
    return;
}

//two pendant rule
void twoPendantRule(std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g)
{
    std::vector<int> gates2pendants;
    std::vector<int> pendantEssential;
    int numpendant;
    for (int i = 0; i < bitsets.size(); ++i)
    {
        numpendant=0;
        for (int j = 0; j < bitsets.size(); ++j)
        {
            if (bitsets[i][j] && i!=j)
            {
                if(bitsets[j].count()==2)
                {
                    numpendant++;
                }

            }
        }
        if(numpendant>1)
            pendantEssential.push_back(i);
    }
    cout<<"Number pendantEssential: "<<pendantEssential.size()<<"\n";
    return;
}

//Rule 4 in paper implementation of essential rule
bool essentialRule2(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g, int v)
{
    //bool isEssential;
    //boost::dynamic_bitset<> subgraphvertices(g.numVertices()),nonsubgraphvertices(g.numVertices()),coveredvertices(g.numVertices()),uncovered(g.numVertices());

    //std::vector<int> subgraph;

    std::vector<int> neighborhoods;
    int numV=g.numVertices();
    for (int i = 0; i < classes.propersubset[v].size()-1; ++i)
    {
        for (int j = i+1; j < classes.propersubset[v].size(); ++j)                                          //check all pairs of proper subset vertices of v
        {
            if ((bitsets[classes.propersubset[v][i]] & bitsets[classes.propersubset[v][j]]).count()==1)     //check if only common neighbor of i and j is v
            {
                //cout<<"two proper subsets\n";
                neighborhoods.clear();
                for (int k = 0; k < numV; ++k)
                {
                    if (k!=v)
                    {
                        if (bitsets[classes.propersubset[v][i]][k])
                        {
                            neighborhoods.push_back(k);
                        }
                        if (bitsets[classes.propersubset[v][j]][k])
                        {
                            neighborhoods.push_back(k);
                        }
                    }
                }
                if (includes(classes.propersubset[v].begin(), classes.propersubset[v].end(),neighborhoods.begin(), neighborhoods.end()))    //check if all common neighborhoods are proper subsets of v
                {

                    return true;
                }
            }   
        }
    }
    return false;
}



//Rule 3 in paper perform intermittent swap on unclassiied MDS vertices
void intermittentSwap(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g)
{
    int numV=g.numVertices();
    size_t index,index2;
    boost::dynamic_bitset<> private_neighbors;
    index=classes.MDS.find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        private_neighbors=bitsets[index];
    //deterime if remaining MDS vertices are essential or intermittent     
        index2=classes.MDS.find_first();
        while(index2!=boost::dynamic_bitset<>::npos)                            //determine private neighbors of vertex in MDS
        {                                       
            if (index!=index2)
            {
                private_neighbors=private_neighbors-bitsets[index2];
            }

            index2=classes.MDS.find_next(index2);
        }
        for (int i = 0; i < numV; ++i)
        {
            if (i!=index && private_neighbors.is_subset_of(bitsets[i]))         //check if private neighbors are subset of another vertex if so both vertices are intermittent
            {
                if (!classes.done[index])
                {
                    classes.done[index]=true;
                    classes.intswap++;
                    classes.intermittent[index]=true;
                    //cout<<"classifying by orbits\n";
                    classes.intswap+=classifyOrbit(classes,index,1);
                    //cout<<"done classifying by orbits\n";
                }
                if (!classes.done[i])
                {
                    classes.done[i]=true;
                    classes.intswap++;
                    classes.intermittent[i]=true;
                    //cout<<"classifying by orbits\n";
                    classes.intswap+=classifyOrbit(classes,i,1);
                }
            }
        }

        index=classes.MDS.find_next(index);
    }
    
    return;
}

//call gurobi and solve model
int solve_model(GRBModel& model)        
{
    model.optimize();
    return model.get(GRB_IntAttr_Status);
}


//for unclassified MDS vertices solve ILP-exclude to determine essential or intermittent
void essential_or_intermittent(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, GRBModel& model, std::vector<float>& MDS_compute_times)
{
    std::chrono::time_point<std::chrono::system_clock> start, end; 
    int num_classified;
    int status;
    size_t index;
    index=classes.MDS.find_first();
    while(index!=boost::dynamic_bitset<>::npos)
    {
        if (!classes.done[index])
        {
            classes.variables[index].set(GRB_DoubleAttr_UB,0.0);

            start=std::chrono::system_clock::now();
            status=solve_model(model);

            end=std::chrono::system_clock::now();
            std::chrono::duration<float> elapsed_seconds = end - start;
            MDS_compute_times.push_back(elapsed_seconds.count());
            classes.MDStest++;
            classes.variables[index].set(GRB_DoubleAttr_UB,1.0);

            if (status==2 && classes.MDSsize==model.get(GRB_DoubleAttr_ObjVal))
            {   
                classes.intermittent[index]=true;
                classes.done[index]=true;
                num_classified=classifyOrbit(classes,index,1);
                classes.byTrueOrbit+=num_classified;
            }
            else
            {
                classes.essential[index]=true;
                classes.variables[index].set(GRB_DoubleAttr_LB,1.0);
                classes.done[index]=true;
                num_classified=classifyOrbit(classes,index,2);
                classes.byTrueOrbit+=num_classified;
                essentialSubset(classes,index);
            }
        }
        index=classes.MDS.find_next(index);
    }

    return;
}

//Rule 2 in paper determine if all vertices in a closed neighborhood are adjacent to an essential vertex
void redundantRule(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets)
{
    int numV=bitsets.size();
    
    bool allAdjacentToEssential;
    size_t index;
    for (int i = 0; i < numV; ++i)
    {
        if (!classes.done[i]&&!classes.MDS[i])
        {
            allAdjacentToEssential=true;
            index=bitsets[i].find_first();
            while(index!=boost::dynamic_bitset<>::npos) 
            {
                if ((bitsets[index]&classes.essential).count()==0)
                {
                    allAdjacentToEssential=false;
                    break;
                }
                index=bitsets[i].find_next(index);
            }
            if (allAdjacentToEssential)
            {
                classes.redundant[i]=true;
                classes.byredundantrule++;
                classes.done[i]=true;
                classes.byredundantrule+=classifyOrbit(classes,i,0);
                classes.variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    return;
}

//for unclassified vertices solve ILP-include to determine redundant or intermittent
void redundant_or_intermittent(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, GRBModel& model,std::vector<float>& MDS_compute_times)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int status;
    int numV=bitsets.size();
    for (int i = 0; i < numV; ++i)
    {
        if (!classes.done[i]&&!classes.MDS[i])
        {
            classes.variables[i].set(GRB_DoubleAttr_LB,1.0);
            start=std::chrono::system_clock::now();
            status=solve_model(model);
            end=std::chrono::system_clock::now();
            std::chrono::duration<float> elapsed_seconds = end - start;
            MDS_compute_times.push_back(elapsed_seconds.count());
            classes.MDStest++;
            classes.variables[i].set(GRB_DoubleAttr_LB,0.0);
            if (status==2 && classes.MDSsize==model.get(GRB_DoubleAttr_ObjVal))
            {
                classes.intermittent[i]=true;
                classes.done[i]=true;
                classes.byTrueOrbit+=classifyOrbit(classes,i,1);
            }
            else
            {
                classes.redundant[i]=true;
                classes.done[i]=true;
                classes.byTrueOrbit+=classifyOrbit(classes,i,0);   
                classes.variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    return;
}






int main(int argc, char const *argv[])
{
	ifstream graphFile(argv[1]);                                       //file with edge list of graph
	ofstream outputFile(argv[2]);

	Graph g;                                                           //create graph object
    cout<<"reading graph\n";
	g.readGraph(graphFile);                                            //read in graph from file
	graphFile.close();
    cout<<g.numVertices()<<" vertices and "<<g.numEdges()<<" edges\n";
	GRBEnv* env=NULL;                                                  //initialize Gurobi
	env = new GRBEnv();                                                //set up Gurobi enviornment

    GRBModel model = GRBModel(*env);                                   //initialize Gurobi model object

    const int numV=g.numVertices();

    std::vector<boost::dynamic_bitset<> > bitsets(numV);               //adjacency matrix
    
    std::vector<GRBConstr> constraints(numV);                          //ILP contraints
    struct Classifications classes;                                    //struct to hold all MDS and class information
    classes.MDS.resize(numV);                                          //initial MDS as bitset
    classes.essential.resize(numV);                                    //essential vertices as bitset
    classes.intermittent.resize(numV);                                 //intermittent vertices as bitset
    classes.redundant.resize(numV);                                    //redundant vertices as bitset
    classes.done.resize(numV);                                         //classified vertices as bitset
    classes.propersubset.resize(numV);                                 //bitset for vertices that are a proper subset
    classes.propersuperset.resize(numV);                               //bitset for vertices that are a proper superset
    classes.subset.resize(numV);                                       //bitset for verticest that are a subset
    classes.variables.resize(numV);                                    //gurobi variable information
    classes.orbitComputed.resize(numV);                                //vertices classified by orbits

    std::vector<std::vector<int> > superset(numV);                     //vertices that are a superset
    std::vector<int> numSuperset(numV,0);                              //number of supersets per vertex

    int numEss,numInt,numRed;
    numEss=numInt=numRed=0;

    createModel(model,bitsets,classes,g);                              //create Gurobi ILP model

    std::chrono::time_point<std::chrono::system_clock> start, end;     //timing variables for various portions of code
    start=std::chrono::system_clock::now();                            
    //classes.orbits=createOrbits(bitsets,g.numEdges());                 //generate orbits using nauty
    classes.orbits=createOrbits_saucy(g,bitsets);                      //generate orbits using saucy
    end=std::chrono::system_clock::now();
    std::chrono::duration<float> elapsed_seconds = end - start;
    cout<<"Orbits time: "<<elapsed_seconds.count() << "s\n";
    //twoPendantRule(bitsets,g);
    int byessrule2=0;
    bool ess;
    std::vector<std::vector<int> > subcomps;

    std::vector<float> MDS_compute_times;
    
    start=std::chrono::system_clock::now();
    for (int i = 0; i < numV; ++i)
    {
        if (!classes.done[i] && classes.propersubset[i].size()>0 && classes.propersuperset[i].size()==0)
        {
            if(essentialRule2(classes,bitsets,g,i))                    //determine if vertex is essential (was called 2 because second iteration of essential rule)
            {
                //cout<<i<<" is essential by rule 2.0\n";
                classes.essential[i]=true;                             //set to essential
                classes.variables[i].set(GRB_DoubleAttr_LB,1.0);       //set vertex's varialbe to 1 so in all solutions
                classes.byessentialrule++;
                classes.done[i]=true;
                essentialSubset(classes,i);                            //perform rule 1
                classes.byessentialrule+=classifyOrbit(classes,i,2);
                byessrule2++;
            }   
        }
    }
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout<<"Essential rule 2 time: "<<elapsed_seconds.count() << "s\n";

    start=std::chrono::system_clock::now();
    model.optimize();                                                  //compute initial MDS
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    MDS_compute_times.push_back(elapsed_seconds.count());
    classes.MDStest++;
    //cout<<MDStest<<" out of "<<numV<<" done\n";

    classes.MDSsize=model.get(GRB_DoubleAttr_ObjVal);
    cout<<"MDSsize: "<<classes.MDSsize<<"\n";

    for (int i = 0; i < numV; ++i)
    {
    	if (classes.variables[i].get(GRB_DoubleAttr_X)==1)            //get MDS vertices
    	{
    		classes.MDS[i]=true;
    	}
    }
    cout<<"intermittentSwap rule\n";
    start=std::chrono::system_clock::now();
    intermittentSwap(classes,bitsets,g);                              //perform rule 3
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout<<"Intermittent swap time: "<<elapsed_seconds.count() << "s\n";
    //cout<<classes.intswap<<"\n";
    cout<<"essential_or_intermittent\n";
    essential_or_intermittent(classes,bitsets,model,MDS_compute_times); //deterime if remaining MDS vertices are essential or intermittent
    cout<<"redundantRule\n";
    start=std::chrono::system_clock::now();
    redundantRule(classes,bitsets);                                     //perform Rule 2
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout<<"Redundant rule time: "<<elapsed_seconds.count() << "s\n";
    cout<<"redundant_or_intermittent\n";
    redundant_or_intermittent(classes,bitsets,model,MDS_compute_times); //determine if remaining vertices are redundant or intermittent


    //print out runtime metrics and number of each class
    float total_MDS_time, avg_MDS_time;
    total_MDS_time=0;
    float max_time=MDS_compute_times[0];
    float min_time=MDS_compute_times[0];
    for (int i = 0; i < MDS_compute_times.size(); ++i)
    {
        if (MDS_compute_times[i]>max_time) max_time=MDS_compute_times[i];
        if (MDS_compute_times[i]<min_time) min_time=MDS_compute_times[i];
        total_MDS_time+=MDS_compute_times[i];
    }
    cout<<"Total MDS time: "<<total_MDS_time<<"\n";
    cout<<"Average MDS time: "<<total_MDS_time/MDS_compute_times.size()<<"\n";
    cout<<"Max MDS compute time: "<<max_time<<"\n";
    cout<<"Min MDS compute time: "<<min_time<<"\n\n";



    cout<<"Number Essential: ";//<<numEss<<"\n";
    cout<<classes.essential.count()<<"\n";
    cout<<"Number Redundant: ";//<<numRed<<"\n";
    cout<<classes.redundant.count()<<"\n";
    cout<<"Number Intermittent: ";//<<numInt<<"\n";
    cout<<classes.intermittent.count()<<"\n";
    cout<<"Number by true orbit: "<<classes.byTrueOrbit<<"\n";
    cout<<"Number by orbit: "<<classes.byOrbit<<"\n";
    cout<<"Number by essential rule: "<<classes.byessentialrule<<"\n";
    cout<<"Number by subset: "<<classes.bySubset<<"\n";
    cout<<"Number by intermittent swap: "<<classes.intswap<<"\n";
    cout<<"Number by redundant rule: "<<classes.byredundantrule<<"\n";
    cout<<"Number of MDS tests: "<<classes.MDStest<<"\n";
    cout<<"\n |V|: "<<numV<<" |C| "<<classes.essential.count()+classes.redundant.count()+classes.intermittent.count()<<"\n";

    return 0;
}
