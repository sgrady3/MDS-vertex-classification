// g++ -m64 -std=c++11 -o essential-DS-nacher Graph.cpp essential-DS-nacher.cpp computeOrbits.cpp /home/stephen/code/nauty26r11/nauty.a -I /home/stephen/code/nauty26r11 -I /home/stephen/code/gurobi752/linux64/include/ -L /home/stephen/code/gurobi752/linux64/lib/ -lgurobi_c++ -lgurobi75 -lm


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
#include "computeOrbits.h"


using namespace std;

struct Classifications
{
   int MDSsize=0;
   int byOrbit=0;
   int byessentialrule=0;
   int bySubset=0;
   int intswap=0;
   int byredundantrule=0;
   int MDStest=0;
   boost::dynamic_bitset<> MDS;
   boost::dynamic_bitset<> essential;
   boost::dynamic_bitset<> intermittent;
   boost::dynamic_bitset<> redundant; 
   boost::dynamic_bitset<> done;
   boost::dynamic_bitset<> orbitComputed;
   std::vector<std::vector<int> > orbits;
   std::vector<std::vector<int> > propersubset;
   std::vector<std::vector<int> > propersuperset;
   std::vector<boost::dynamic_bitset<> > subset;
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

void lowDegree(Graph& g, int k, boost::dynamic_bitset<>& currentMDS, boost::dynamic_bitset<>& redundant)
{
    const int numV=g.numVertices();
    int max=g.maxDegree();

    for (int i = 0; i < numV; ++i)
    {
        //cout<<(k-1)*(max+1)+(g.degree(i)+1)<<" : "<<numV<<"\n";
        if (!currentMDS[i] && (k-1)*(max+1)+(g.degree(i)+1)<numV)
        {
            redundant[i]=true;
        }
    }
    return;
}

std::vector<int> KHighest(Graph& g, int k)
{
    const int numV=g.numVertices();
    std::vector<int> degrees(numV);
    for (int i = 0; i < numV; ++i)
    {
        degrees[i]=g.degree(i);
    }
    partial_sort(degrees.begin(),degrees.begin()+k,degrees.end(),maxSort);
    
    /*for (int i = 0; i < numV; ++i)
    {
        cout<<degrees[i]<<" ";
    }
    cout<<"\n";*/
    return degrees;
}

std::vector<std::vector<int> > createOrbits(std::vector<boost::dynamic_bitset<> >& bitsets)
{
    std::vector<std::vector<int> > orbitSets;
    cout<<"starting orbits\n";
    orbitSets=orbits(bitsets);
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

void createModel(GRBModel& model, std::vector<boost::dynamic_bitset<> >& bitsets, Classifications& classes, Graph& g)
{

    int numV=bitsets.size();
    for (int i = 0; i < numV; ++i)
    {
        classes.subset[i].resize(numV);
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
                classes.subset[i][index]=true;
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


void classifyOrbit(Classifications& classes, int v, int classification)
{
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
                }
                else if (classification==1)
                {
                    classes.intermittent[classes.orbits[v][i]]=true;
                    classes.done[classes.orbits[v][i]]=true;
                    classes.orbitComputed[classes.orbits[v][i]]=true;
                    classes.byOrbit++;
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
                    //cout<<"done\n";
                }
                else{
                    cout<<"error in classifyOrbit, improper classification given\n";
                }
            }
        }           
    }
    return ;
}


void essentialSubset(Classifications& classes, int v)
{
    size_t index;
    index=classes.subset[v].find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        if(!classes.done[index])
        {
            classes.redundant[index]=true;
            classes.bySubset++;
            classes.done[index]=true;
            classifyOrbit(classes,index,0);
            classes.variables[index].set(GRB_DoubleAttr_UB,0.0);
        }
        index=classes.subset[v].find_next(index);
    }
    return;
}

int getMaxDegree(boost::dynamic_bitset<>& subgraph, std::vector<boost::dynamic_bitset<> >& bitsets)
{
    int maxDegree=0;
    int degree;
    size_t index;
    index=subgraph.find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        degree=(subgraph&bitsets[index]).count();
        if (degree>maxDegree)
        {
            maxDegree=degree;
        }
        index=subgraph.find_next(index);
    }
    return maxDegree;
}


void essentialRule(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g, boost::dynamic_bitset<>& marked)
{
    int numV=bitsets.size();
    int psize,maxDegree, subgraphsize;
    bool isEssential;
    std::vector<int>  subgraph,disjointEssential;
    for (int i = 0; i < numV; ++i)
    {
        isEssential=false;
        subgraph.clear();
        psize=classes.propersubset[i].size();
        if (psize>1 && classes.propersuperset[i].size()==0)
        {
            for (int j = 0; j < psize; ++j)
            {
                if (classes.propersuperset[classes.propersubset[i][j]].size()==1)
                {
                    //cout<<classes.propersuperset[classes.propersubset[i][j]].size()<<"\n";
                    subgraph.push_back(classes.propersubset[i][j]);
                    //cout<<"pushing back\n";
                }
            }
            subgraphsize=subgraph.size();
            //printVec(subgraph);
            for (int k = 0; k <subgraphsize-1 ; ++k)
            {
                if (!isEssential)
                {
                    for (int p = k+1; p < subgraphsize; ++p)
                    {
                    //cout<<g.label(subgraph[k])<<" "<<g.label(subgraph[p])<<"\n";
                    //cout<<"intersection size: "<<(bitsets[subgraph[k]]&bitsets[subgraph[p]]).count()<<"\n";
                        if ((bitsets[subgraph[k]]&bitsets[subgraph[p]]).count()==1)
                        {
                            //cout<<g.label(i)<<" is essential with degree "<<bitsets[i].count()<<"\n";
                            disjointEssential.push_back(i);
                            marked[i]=true;
                            isEssential=true;
                            break;
                        }
                    }
                }
                else
                {
                    break;
                }
            }
        }
    }
    cout<<"Number disjointEssential: "<<disjointEssential.size()<<"\n";
 return;
}

bool twoPendantRule(std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g, int v)
{
    int numpendant;
    numpendant=0;
    for (int j = 0; j < bitsets.size(); ++j)
    {
        if (bitsets[v][j] && v!=j)
        {
            if(bitsets[j].count()==2)
            {
                numpendant++;
            }

        }
    }
    if(numpendant>1)
        return true;
    return false;
}

bool essentialRule2(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets, Graph& g, int v)
{
    bool isEssential;
    boost::dynamic_bitset<> subgraphvertices(g.numVertices()),nonsubgraphvertices(g.numVertices()),coveredvertices(g.numVertices()),uncovered(g.numVertices());
    std::vector<int> subgraph;
    for (int i = 0; i < classes.propersubset[v].size(); ++i)
    {
        subgraph.push_back(classes.propersubset[v][i]);
        subgraphvertices[classes.propersubset[v][i]]=true;
    }
    //cout<<v<<"'s subvertices ";
    //printVec(subgraph);
    nonsubgraphvertices=(g.neighbors(v)^subgraphvertices);
    //cout<<"for "<<v<<" these are not subsets: ";
    for (int i = 0; i < g.numVertices(); ++i)
    {
        if (nonsubgraphvertices[i])
        {
            //cout<<i<<" ";
            coveredvertices=(coveredvertices|(g.closedNeighbors(i)&subgraphvertices));
        }
    }
    //cout<<"\n";
    uncovered=(subgraphvertices^coveredvertices);
    std::vector<int> uncoverednames;
    for (int i = 0; i < g.numVertices(); ++i)
    {
        if (uncovered[i])
        {
            uncoverednames.push_back(i);
        }
    }
    //cout<<"for "<<v<<" these are uncovered: ";
    //printVec(uncoverednames);
    if (uncovered.count()>1)
    {
        isEssential=true;
        for (int i = 0; i < g.numVertices(); ++i)
        {
            if(subgraphvertices[i] && (uncovered.is_subset_of(g.closedNeighbors(i))))
            {
                isEssential=false;
            }
        }
        return isEssential;
    }
    else
    {
        return false;
    }

}


int dfsUtil(Graph& g, boost::dynamic_bitset<>& subgraphvertices, boost::dynamic_bitset<>& visited, std::vector<int>& component, int v)
{
    boost::dynamic_bitset<> neighbohood=g.neighbors(v);
    int privateFlag;
    if((g.closedNeighbors(v)&subgraphvertices).count()==g.closedNeighbors(v).count())
    {
        visited[v]=true;
        component.push_back(v);

        for (int i = 0; i < g.numVertices(); ++i)
        {
            if(neighbohood[i] && !visited[i])
            {
                privateFlag=dfsUtil(g,subgraphvertices,visited,component,i); 
                if (privateFlag==-1)
                {
                    return -1;
                }
                else
                {
                    return 1;
                }
            }
        }
        return 1;
    }
    else
    {
        return -1;
    }
}


std::vector<std::vector<int> > subsetComponents(Classifications& classes, Graph& g, int v)
{
    std::vector<std::vector<int> > subComponents;
    boost::dynamic_bitset<> subgraphvertices(g.numVertices());
    std::vector<int> candidate, component;
    boost::dynamic_bitset<> visited(g.numVertices());
    
    for (int i = 0; i < classes.propersubset[v].size(); ++i)
    {
        candidate.push_back(classes.propersubset[v][i]);
    }
    int candidatesize=candidate.size();

    //printVec(candidate);
    visited[v]=true;
    subgraphvertices[v]=true;
    for (int i = 0; i < candidatesize; ++i)
    {
        subgraphvertices[candidate[i]]=true;
    }
    int privateFlag;
    for (int i = 0; i < candidatesize; ++i)
    {
        //cout<<"i "<<i<<"\n";
        //cout<<"out "<<candidate[i]<<" ";
        if(!visited[candidate[i]])
        {
            //cout<<"in "<<candidate[i]<<"\n";
            privateFlag=dfsUtil(g,subgraphvertices,visited,component,candidate[i]);
            if (privateFlag==1)
            {
                subComponents.push_back(component);
            }
            component.clear();
        }
    }
    return subComponents;
}

void intermittentSwap(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets)
{
    //cout<<"doing the intermittent swap\n";
    int numV=bitsets.size();
    std::vector<boost::dynamic_bitset<> > exclusiveDominated(bitsets.size(),boost::dynamic_bitset<> (bitsets.size()));
    boost::dynamic_bitset<> exclusiveNeigh;
    size_t index,index2,index3;
    index=classes.MDS.find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        index2=bitsets[index].find_first();
        while(index2!=boost::dynamic_bitset<>::npos) 
        {
            exclusiveNeigh=bitsets[index2]&classes.MDS;
            if(exclusiveNeigh.count()==1)
            {
                exclusiveDominated[index][index2]=true;
            }
            index2=bitsets[index].find_next(index2);
        }
        index2=exclusiveDominated[index].find_first();
        while(index2!=boost::dynamic_bitset<>::npos)
        {
            index3=bitsets[index2].find_first();
            while(index3!=boost::dynamic_bitset<>::npos)
            {
                if(index!=index3 && exclusiveDominated[index].is_subset_of(bitsets[index3]))
                {
                    if (!classes.done[index])
                    {
                        classes.done[index]=true;
                        classes.intswap++;
                        classes.intermittent[index]=true;
                        //cout<<"classifying by orbits\n";
                        classifyOrbit(classes,index,1);
                        //cout<<"done classifying by orbits\n";
                    }
                    if (!classes.done[index3])
                    {
                        classes.done[index3]=true;
                        classes.intswap++;
                        classes.intermittent[index3]=true;
                        //cout<<"classifying by orbits\n";
                        classifyOrbit(classes,index3,1);
                        //cout<<"done classifying by orbits\n";
                        //classes.numInt++;
                    }
                }
                index3=bitsets[index2].find_next(index3);
            }
            index2=exclusiveDominated[index].find_next(index2);
        }

        index=classes.MDS.find_next(index);
    }
    
    return;
}

int solve_model(GRBModel& model)
{
    model.optimize();
    return model.get(GRB_IntAttr_Status);
}


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
                //num_classified=classifyOrbit(classes,index,1);
                //classes.byTrueOrbit+=num_classified;
            }
            else
            {
                classes.essential[index]=true;
                classes.variables[index].set(GRB_DoubleAttr_LB,1.0);
                classes.done[index]=true;
                //num_classified=classifyOrbit(classes,index,2);
                //classes.byTrueOrbit+=num_classified;
                //essentialSubset(classes,index);
            }
        }
        index=classes.MDS.find_next(index);
    }

    return;
}

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
                //classifyOrbit(classes,i,0);
                classes.variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    return;
}

void allNeighborsEssential(Classifications& classes, std::vector<boost::dynamic_bitset<> >& bitsets,Graph& g)
{
    int numV=bitsets.size();
    bool allAdjacentEssential;
    size_t index;
    for (int i = 0; i < numV; ++i)
    {
        if (!classes.done[i]&&!classes.MDS[i])
        {
            if ((g.neighbors(i)&classes.essential).count()==g.neighbors(i).count())
            {
                classes.redundant[i]=true;
                classes.byredundantrule++;
                classes.done[i]=true;
                //classifyOrbit(classes,i,0);
                classes.variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    return;
}


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
                //classes.byTrueOrbit+=classifyOrbit(classes,i,1);
            }
            else
            {
                classes.redundant[i]=true;
                classes.done[i]=true;
                //classes.byTrueOrbit+=classifyOrbit(classes,i,0);   
                classes.variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    return;
}






int main(int argc, char const *argv[])
{
    ifstream graphFile(argv[1]);
    ofstream outputFile(argv[2]);

    Graph g;
    cout<<"reading graph\n";
    g.readGraph(graphFile);
    graphFile.close();
    cout<<g.numVertices()<<" vertices and "<<g.numEdges()<<" edges\n";
    GRBEnv* env=NULL;
    env = new GRBEnv();

    GRBModel model = GRBModel(*env);

    const int numV=g.numVertices();
    //std::vector<GRBVar> variables(numV);
    std::vector<boost::dynamic_bitset<> > bitsets(numV);
    //std::vector<boost::dynamic_bitset<> > subset(numV);
    //std::vector<boost::dynamic_bitset<> > propersubset(numV);
    std::vector<GRBConstr> constraints(numV);
    struct Classifications classes;
    classes.MDS.resize(numV);
    classes.essential.resize(numV);
    classes.intermittent.resize(numV);
    classes.redundant.resize(numV);
    classes.done.resize(numV);
    classes.propersubset.resize(numV);
    classes.propersuperset.resize(numV);
    classes.subset.resize(numV);
    classes.variables.resize(numV);
    classes.orbitComputed.resize(numV);
    //boost::dynamic_bitset<> currentMDS(numV);
    //boost::dynamic_bitset<> essential(numV);
    //boost::dynamic_bitset<> intermittent(numV);
    //boost::dynamic_bitset<> redundant(numV);
    //boost::dynamic_bitset<> done(numV);
    //boost::dynamic_bitset<> orbitComputed(numV);
    std::vector<std::vector<int> > superset(numV);
    std::vector<int> numSuperset(numV,0);

    int numEss,numInt,numRed;
    numEss=numInt=numRed=0;

    createModel(model,bitsets,classes,g);
    //classes.orbits=createOrbits(bitsets);
    boost::dynamic_bitset<> marked(numV);
    //essentialRule(classes,bitsets,g,marked);
    //twoPendantRule(bitsets,g);
    int bycomponent=0;
    int byessrule2=0;
    bool ess;
    std::vector<std::vector<int> > subcomps;
    marked.reset();

    std::vector<float> MDS_compute_times;
    std::chrono::time_point<std::chrono::system_clock> start, end; 
    start=std::chrono::system_clock::now();
    
    for (int i = 0; i < numV; ++i)
    {
        if (!marked[i] && classes.propersubset[i].size()>0 && classes.propersuperset[i].size()==0)
        {
            if(twoPendantRule(bitsets,g,i))                                 //perform two pendant rule
            {
                //cout<<i<<" is essential by rule 2.0\n";
                classes.essential[i]=true;
                classes.variables[i].set(GRB_DoubleAttr_LB,1.0);
                classes.byessentialrule++;
                classes.done[i]=true;
                //essentialSubset(classes,i);
                //classifyOrbit(classes,i,2);
                byessrule2++;
            }   
        }
    }
    end=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout<<"Two pendant rule time: "<<elapsed_seconds.count() << "s\n";

    cout<<"by essential rule 2: "<<byessrule2<<"\n";

    start=std::chrono::system_clock::now();
    model.optimize();
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    MDS_compute_times.push_back(elapsed_seconds.count());
    classes.MDStest++;
    //cout<<MDStest<<" out of "<<numV<<" done\n";

    classes.MDSsize=model.get(GRB_DoubleAttr_ObjVal);
    cout<<"MDSsize: "<<classes.MDSsize<<"\n";

    for (int i = 0; i < numV; ++i)
    {
        if (classes.variables[i].get(GRB_DoubleAttr_X)==1)
        {
            classes.MDS[i]=true;
        }
    }
    //cout<<"intermittentSwap rule\n";
    //intermittentSwap(classes,bitsets);
    //cout<<classes.intswap<<"\n";
    cout<<"essential_or_intermittent\n";
    essential_or_intermittent(classes,bitsets,model,MDS_compute_times);
    cout<<"redundantRule\n";
    start=std::chrono::system_clock::now();
    allNeighborsEssential(classes,bitsets,g);                                   //perform all neighbors essential rule
    end=std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout<<"All neighbor essential time: "<<elapsed_seconds.count() << "s\n";
    cout<<"redundant_or_intermittent\n";
    redundant_or_intermittent(classes,bitsets,model,MDS_compute_times);

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
    cout<<"Number by orbit: "<<classes.byOrbit<<"\n";
    cout<<"Number by essential rule: "<<classes.byessentialrule<<"\n";
    cout<<"Number by subset: "<<classes.bySubset<<"\n";
    cout<<"Number by intermittent swap: "<<classes.intswap<<"\n";
    cout<<"Number by redundant rule: "<<classes.byredundantrule<<"\n";
    cout<<"Number of MDS tests: "<<classes.MDStest<<"\n";
    cout<<"\n |V|: "<<numV<<" |C| "<<classes.essential.count()+classes.redundant.count()+classes.intermittent.count()<<"\n";

#if 0
    

    

    //**********************************************************************************************//
    //**************************check redunancy by essential neighbor rule**************************//
    //**********************************************************************************************//
    boost::dynamic_bitset<> essentialNeigh(numV);
    bool allAdjacentToEssential;
    int numEssNeighTests=0;
    for (int i = 0; i < numV; ++i)
    {
        if (!done[i])
        {
            allAdjacentToEssential=true;
            index=bitsets[i].find_first();
            while(index!=boost::dynamic_bitset<>::npos) 
            {
                essentialNeigh.clear();
                essentialNeigh=bitsets[index]&essential;
                //cout<<essentialNeigh<<"\n";
                if (essentialNeigh.count()==0)
                {
                    allAdjacentToEssential=false;
                    //cout<<"breaking out\n";
                    break;

                }
                index=bitsets[i].find_next(index);
            }
            if(allAdjacentToEssential)
            {
                //cout<<"found redundant by essential neighbor\n";
                numEssNeighTests++;
                redundant[i]=true;
                done[i]=true;
                numRed++;
                variables[i].set(GRB_DoubleAttr_UB,0.0);
                for (int j = 0; j < orbitSets[i].size(); ++j)
                {

                    if (!done[orbitSets[i][j]])
                    {
                        redundant[orbitSets[i][j]]=true;
                        orbitComputed[orbitSets[i][j]]=true;
                        byOrbit++;
                        numRed++;
                        done[orbitSets[i][j]]=true;
                        variables[orbitSets[i][j]].set(GRB_DoubleAttr_UB,0.0);
                    }   
                }
            }            
        }
    }


    //orbitComputed.reset();
    //get redundant or interittent vertices by orbits out of solution
    for (int i = 0; i < numV; ++i)
    {
        if (!done[i] && !orbitComputed[i] && !currentMDS[i] && orbitSets[i].size()>1)
        {
            //cout<<"testing for redundant/intermittent orbit: "<<i<<"\n";
            variables[i].set(GRB_DoubleAttr_LB,1.0);
            model.optimize();
            MDStest++;
            //cout<<MDStest<<" out of "<<numV<<" done\n";
            variables[i].set(GRB_DoubleAttr_LB,0.0);
            if (MDSsize<model.get(GRB_DoubleAttr_ObjVal))
            {
                redundant[i]=true;
                done[i]=true;
                numRed++;
                variables[i].set(GRB_DoubleAttr_UB,0.0);
                for (int j = 0; j < orbitSets[i].size(); ++j)
                {

                    if (!done[orbitSets[i][j]])
                    {
                        redundant[orbitSets[i][j]]=true;
                        orbitComputed[orbitSets[i][j]]=true;
                        byOrbit++;
                        numRed++;
                        done[orbitSets[i][j]]=true;
                        variables[orbitSets[i][j]].set(GRB_DoubleAttr_UB,0.0);
                    }
                    
                }

            }
            else
            {
                intermittent[i]=true;
                done[i]=true;
                numInt++;
                for (int j = 0; j < orbitSets[i].size(); ++j)
                {
                    if (!currentMDS[orbitSets[i][j]] && !done[orbitSets[i][j]])
                    {
                        intermittent[orbitSets[i][j]]=true;
                        orbitComputed[orbitSets[i][j]]=true;
                        byOrbit++;
                        numInt++;
                        done[orbitSets[i][j]]=true;
                    }
                }

            }
            orbitComputed[i]=true;
        }
    }
    
     

    /*lowDegree(g,MDSsize,currentMDS,redundant);

    index=redundant.find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        cout<<index<<" ";
        index=redundant.find_next(index);
    }
    cout<<"\n";

    std::vector<int> sortedDeg=KHighest(g,MDSsize);
    */

    

    for (int i = 0; i < numV; ++i)
    {
        if (!done[i])
        {
            //cout<<"finding redundant\n";
            variables[i].set(GRB_DoubleAttr_LB,1.0);
            model.optimize();
            MDStest++;
            //cout<<MDStest<<" out of "<<numV<<" done\n";
            variables[i].set(GRB_DoubleAttr_LB,0.0);
            //cout<<model.get(GRB_DoubleAttr_ObjVal)<<"\n";
            if (MDSsize<model.get(GRB_DoubleAttr_ObjVal))
            {
                //cout<<i<<" is redundant\n";
                redundant[i]=true;
                numRed++;
                variables[i].set(GRB_DoubleAttr_UB,0.0);
            }
            else
            {
                //cout<<i<<" is intermittent\n";
                intermittent[i]=true;
                numInt++;
            }
        }
        done[i]=true;
    }

    cout<<"Number Essential: ";//<<numEss<<"\n";
    cout<<essential.count()<<"\n";
    cout<<"Number Redundant: ";//<<numRed<<"\n";
    cout<<redundant.count()<<"\n";
    cout<<"Number Intermittent: ";//<<numInt<<"\n";
    cout<<intermittent.count()<<"\n";
    cout<<"Number by orbit: "<<byOrbit<<"\n";
    cout<<"Number by subset: "<<bySubset<<"\n";
    cout<<"Number by intermittent swap: "<<intswap<<"\n";
    cout<<"Number by redundant test: "<<numEssNeighTests<<"\n";
    cout<<"Number of MDS tests: "<<MDStest<<"\n";

    size_t index1;
    index=essential.find_first();
    index1=intermittent.find_first();
    while(index!=boost::dynamic_bitset<>::npos) 
    {
        if(essential[index])
            cout<<g.label(index)<<" ";
        index=essential.find_next(index);    
    }
    cout<<"\n";

    return 0;
#endif
}