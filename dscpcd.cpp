#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cstdlib>
#include <time.h> 
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#define RndUniDevInt(Range) (rand() % Range) // Randomly generate an integer between 0 and b, [0,b), do not contain b.
#define RndUniDev() (rand() / double(RAND_MAX)) // Ramdomly generate a float between 0 and 1, 0~1.
using namespace std;

typedef vector<int> TIntV;
typedef vector<double> TDblV;
typedef vector<string> TStrV;
typedef vector<vector<int> > TIntVV;
typedef vector<vector<double> > TDblVV;
typedef vector<set<int> > TIntSetV;
typedef vector<set<double> > TDblSetV;
typedef vector<pair<double, int> > TDblIntPrV;
typedef set<int> TIntSet;
typedef set<int>::iterator TIntSetIter;
typedef set<double> TDblSet;
typedef map<int, int> TIntIntMap;

// record the cost time from tod1 to tod2
long long int todiff(struct timeval* tod1, struct timeval* tod2) {
  long long t1, t2;
  t1 = tod1->tv_sec * 1000000 + tod1->tv_usec;
  t2 = tod2->tv_sec * 1000000 + tod2->tv_usec;
  return t1 - t2;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// basic graph class
class GRAPH {
private:
  // dataset name (do not contain the suffix name)
  string DataName;
  /*---------basic topology-----------*/
  // network topology, each set of corresponding node contains its all neighbors
  TIntSetV Neighbor;
  // number of nodes
  int NumNodes;
  // number of edges
  int NumEdges;
  /*--------weighted graph------------*/
  // values on nodes (each node is associated with a single value (e.g., centrality))
  TDblV WgtNodes;
  // values on nodes (each node is associated with a vector,0-1, binary)
  TIntVV BinNodeFea;
  // length of feature of nodes when each node is associated with a vector
  int LenNodeFea;
  /*-------cluster gt information---------*/
  // number of real clusters
  int NumClus;
  // each set of corresponding node contains the cluster ids she belongs
  TIntSetV ClusInNode;
  // each set of corresponding cluster contains all nodes in this cluster
  TIntSetV ClusInClus;

public:
  // constructor (initialization, file name is needed)
  GRAPH(string dataname) { DataName = dataname; }
  // get data name
  string GetDataName() { return DataName; }
  
  // load a graph (directed) with topology information
  void LoadGraphTpl();
  // get graph topology
  TIntSet GetNeighbor(int NID) { return Neighbor[NID-1]; }
  TIntSetV GetNeighbor() { return Neighbor; }
  // get number of nodes
  int GetNumNodes() { return NumNodes; }
  // get number of edges
  int GetNumEdges() { return NumEdges; }
  // get degree of a node
  int GetDeg(int NID) { return Neighbor[NID-1].size(); }
  // get degree of all nodes
  TIntV GetDegV() { 
    TIntV DegAllNodes;
    for (int NID = 1; NID <= GetNumNodes(); NID++) { DegAllNodes.push_back(GetDeg(NID)); } 
    return DegAllNodes;
  }
  
  // get a node's weight (single value)
  double GetWgtNode(int NID) { return WgtNodes[NID-1]; }
  // load a graph node feature information
  void LoadNodeFea();
  // get length of features
  int GetLenNodeFea() { return LenNodeFea; }
  // get node features (in 0-1)
  TIntV GetBinNodeFea(int NID) { return BinNodeFea[NID-1]; }
  
  // load real cluster information
  void LoadClusGt();
  // get gt clusters listed in cluster
  TIntSetV GetClusInClus() { return ClusInClus; }
  // get a set of nodes in a specific cluster 
  TIntSet GetClusInClus(int CID) { return ClusInClus[CID-1]; }
  // get gt clusters listed in node
  TIntSetV GetClusInNode() { return ClusInNode; }
  // get a set of cluster ids to which a specific node belongs
  TIntSet GetClusInNode(int NID) { return ClusInNode[NID-1]; }
  // get number of clusters
  int GetNumClus() { return NumClus; }
  // calculate the number of common communities of all node pairs
  int CalNumCmnCmty(int NId1, int NId2) { return CalNumJointSets(ClusInNode[NId1-1], ClusInNode[NId2-1]); }
  
  // test the functionality of this Graph class
  void FuncTest();
  
  // calculate the number of elements in the intersection set of two sets
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
};

void GRAPH::LoadGraphTpl() {
  string SGraphFile = DataName + ".graph";
  const char* GraphFile = SGraphFile.c_str();
  if (access(GraphFile, R_OK|W_OK) != 0) {
    printf("No topology (.graph) file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finG;
  finG.open(GraphFile);
  string szLine;
  TIntSet tempset;
  while(getline(finG, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      tempset.insert(atoi((*iter).c_str()));
    }
    Neighbor.push_back(tempset);
    tempset.clear();  
  }
  finG.close();
  finG.clear();
  NumNodes = Neighbor.size();
  NumEdges = 0;
  for (int i = 0; i < Neighbor.size(); i++) { NumEdges += Neighbor[i].size(); }
  // each edge is twicely recorded (undirected graph)
  NumEdges /= 2;
}

// load a graph node feature information
void GRAPH::LoadNodeFea() {
  string SFeaFile = DataName + ".nf";
  const char* FeaFile = SFeaFile.c_str();
  if (access(FeaFile, R_OK|W_OK) != 0) {
    cout << "There is no node feature (.nf) file!!!" << endl;
    return;
  }
  // if there exists a .fea file
  int maxFeaId = 0;
  int tempd;
  ifstream finNF;
  finNF.open(FeaFile);
  string szLine; 
  while(getline(finNF,szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<TStrV >(tData));
    // find the maximum feature id
    for(TStrV::iterator iter = tData.begin(); iter!=tData.end(); ++iter) {
      tempd = atoi((*iter).c_str());
      if (maxFeaId <= tempd) { maxFeaId = tempd; }
    }  
  }
  finNF.close();
  finNF.clear();
  // length of each node's features
  LenNodeFea = maxFeaId;
  cout << "Length of node features: " << LenNodeFea << endl;
  BinNodeFea.resize(NumNodes);
  int i, j;
  for (i = 0; i < NumNodes; i++) { 
    for (j = 0; j < LenNodeFea; j++) { 
    BinNodeFea[i].push_back(0); 
  }
  }
  i = 0;
  finNF.open(FeaFile);
  while(getline(finNF,szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter<TStrV >(tData));
    // reopen and write to fea_nodes
    // if this node has a feature, then the value of corresponding position is 1, 0 otherwise.
    for(TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) {
      BinNodeFea[i][atoi((*iter).c_str()) - 1] = 1;
    }
    i++;  
  }
  finNF.close();
  finNF.clear();
}

// load real cluster information
void GRAPH::LoadClusGt() {
  string SClusFile = DataName + ".gt";
  const char* ClusFile = SClusFile.c_str();
  if (access(ClusFile, R_OK|W_OK) != 0) { printf("No cluster (.gt) file exists!!!\n"); return; }
  ifstream finGT;
  finGT.open(ClusFile);
  string szLine;
  TIntSet tempset;
  while(getline(finGT, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) {
      tempset.insert(atoi((*iter).c_str()));
    }
    ClusInClus.push_back(tempset);
    tempset.clear();
  }
  finGT.close();
  finGT.clear();
  NumClus = ClusInClus.size();
  ClusInNode.resize(NumNodes);
  for (int i = 1; i <= ClusInClus.size(); i++) {
    for (TIntSetIter it_set = ClusInClus[i - 1].begin(); it_set != ClusInClus[i - 1].end(); it_set++) {
      ClusInNode[*it_set - 1].insert(i);
    }
  }
}

// test the functionality of this Graph class
void GRAPH::FuncTest() {
  if (NumNodes > 50) {
    printf("Data scale is too big not suitable for full-print, please specilize it!!!!\n");
    return;
  }
  int i, j;
  TIntSetIter it_set;
  if (Neighbor.empty()) { 
    printf("No topology information!!!\n");
    return;  
  } else {
    printf("============== Number of Nodes: %d, Number of edges: %d =================\n", NumNodes, NumEdges);
    for (i = 0; i < Neighbor.size(); i++) {
      printf("node-%d: ", i);
      for (it_set = Neighbor[i].begin(); it_set != Neighbor[i].end(); it_set++) {
        printf("%d ", *it_set);
      }
      printf("\n");
  }
  }
  
  if (BinNodeFea.empty()) { 
    printf("No node feature information!!!\n"); 
  } else {
    printf("================= Length of node feature: %d =====================\n", LenNodeFea);
    for (i = 0; i < NumNodes; i++) {
      printf("node-%d: ", i);
      for (j = 0; j < LenNodeFea; j++) {
        printf("%d ", BinNodeFea[i][j]);
      }
      printf("\n");
    }
  }
  
  if (ClusInClus.empty()) { 
    printf("No ground-truth information!!!\n"); 
  } else {
    printf("============ Number of clusters: %d ==================\n", NumClus);
    for (i = 0; i < ClusInClus.size(); i++) {
      printf("cluster-%d: ", i + 1);
      for (it_set = ClusInClus[i].begin(); it_set != ClusInClus[i].end(); it_set++) {
        printf("%d ", *it_set);
      }
      printf("\n");
    }
  }
}


//////////////////////////////////
///////////////////////////////////////////////
// class of evaluating a graph's basic characteristics
class EVALGRAPH {
public:
  TIntV GetMaxDeg(GRAPH G); // return index from 1 and maxvalue 
  TIntV GetMinDeg(GRAPH G); // return index from 1 and minvalue 
  int GetAvgDeg(GRAPH G); // return average degree of a graph
  TIntIntMap GetDegreeDistribution(GRAPH G); // key is number of degree value, value is number of nodes
}; 

TIntV EVALGRAPH::GetMaxDeg(GRAPH G) {
  TIntV DegAllNodes = G.GetDegV();
  vector<int>::iterator MaxPos = max_element(DegAllNodes.begin(), DegAllNodes.end());
  TIntV MaxPosValue;
  MaxPosValue.push_back(distance(DegAllNodes.begin(), MaxPos) + 1); // position index, begin from 1
  MaxPosValue.push_back(*MaxPos); // value
  cout << "max degree, node-id: " << MaxPosValue[0] << ", value= " << MaxPosValue[1] << endl;
  return MaxPosValue;
}

TIntV EVALGRAPH::GetMinDeg(GRAPH G) {
  TIntV DegAllNodes = G.GetDegV();
  vector<int>::iterator MinPos = min_element(DegAllNodes.begin(), DegAllNodes.end());
  TIntV MinPosValue;
  MinPosValue.push_back(distance(DegAllNodes.begin(), MinPos) + 1); // position index, begin from 1
  MinPosValue.push_back(*MinPos); // value
  cout << "min degree, node-id: " << MinPosValue[0] << ", value= " << MinPosValue[1] << endl;
  return MinPosValue;
}

int EVALGRAPH::GetAvgDeg(GRAPH G) {
  int SumDeg = 0;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) { SumDeg += G.GetDeg(NID); }
  SumDeg /= G.GetNumNodes();
  cout << "average degree: " << SumDeg << endl;
  return SumDeg;
}

TIntIntMap EVALGRAPH::GetDegreeDistribution(GRAPH G) {
  TIntIntMap DegDistirb;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) {
    int NDeg = G.GetDeg(NID);
    //cout << "Node-id: " << NID << ", Degree: " << NDeg << endl; 
    TIntIntMap::iterator iter = DegDistirb.find(NDeg);
    if (iter != DegDistirb.end()) {
      DegDistirb[NDeg] = iter->second + 1; // if the key exists, corresponding value plus 1
    } else {
      DegDistirb.insert(pair<int, int>(NDeg, 1)); // add a new key if the key does not exist, set the value as 1
    }
  }
  /*cout << endl;
  for (TIntIntMap::iterator iter = DegDistirb.begin(); iter != DegDistirb.end(); iter++) { 
    cout << "Degree: " << iter->first << ", Number: " << iter->second << endl; 
  }*/
  return DegDistirb;
}


//////////////////////////////////
///////////////////////////////////////////////
// class of evaluating nodes' centrality
class CENTRALITY {
public:
  // conductance centrality of one node
  double CondOneNode(GRAPH G, int NID);
  // conductance centrality of all nodes in a graph
  TDblV CondAllNodes(GRAPH G);
};

double CENTRALITY::CondOneNode(GRAPH G, int NID) {
// conductance centrality (Gleich, et al. KDD12)
// last modified in 2021-01-14 by Yuyao Wang
  double ImptNode = 0.0; // small value indicates an important node.
  TIntSetV Neighbor = G.GetNeighbor();
  TIntSet NghNId(Neighbor[NID - 1]);
  NghNId.insert(NID); // did in SNAP
  int LenNghs = NghNId.size();
  if (LenNghs < 5) {
    ImptNode = 1.0;
    return ImptNode;
  }
  int Edges2 = 2 * G.GetNumEdges();
  int Vol = 0,  Cut = 0;
  for (TIntSetIter it_set = NghNId.begin(); it_set != NghNId.end(); it_set++) {
    for (TIntSetIter it_set1 = Neighbor[*it_set - 1].begin(); it_set1 != Neighbor[*it_set - 1].end(); it_set1++) {
      // whether her neighbor's neighbor is also her neighbor.
      if (NghNId.find(*it_set1) == NghNId.end()) { Cut += 1; }
    }
    // Vol store the summation of degree of all nodes inside this set 
    Vol += Neighbor[*it_set - 1].size();
  }
  // get conductance
  if (Vol != Edges2) {
    if (2 * Vol > Edges2) { ImptNode = Cut / double (Edges2 - Vol); }
    else if (Vol == 0) { ImptNode = 0.0; }
    else { ImptNode = Cut / double(Vol); }
  } else {
    if (Vol == Edges2) { ImptNode = 1.0; }
  }
  return ImptNode;
}

TDblV CENTRALITY::CondAllNodes(GRAPH G) {
// conductance centrality of all nodes in a graph
// small value indicates an important node.
  TDblV Conductance;
  for (int NID = 1; NID <= G.GetNumNodes(); NID++) {
    Conductance.push_back(CondOneNode(G, NID));
  }
  return Conductance;
}


//////////////////////////////////
///////////////////////////////////////////////
// class of file output
class IOCLASS {
public:
  // input a TwoDimVV (stored in TIntVV).
  void LoadIntVV(string DataName, string Suffix, TIntVV& TwoDimVV);
  // input a TwoDimVV (stored in TDblVV).
  void InputDblVV(string DataName, string Suffix, TDblVV& TwoDimVV);
  
public:
  // output a vector (stored in T[X]V).
  template <class T> void OutputVec(string FileName, string Suffix, vector<T> OneDimVec);
  // output a set (stored in T[X]Set).
  template <class T> void OutputSet(string FileName, string Suffix, set<T> OneDimSet);
  // output a TwoDimVV (stored in T[X]VV).
  template <class T> void OutputVV(string FileName, string Suffix, vector<vector<T> > TwoDimVV);
  // output a TwoDimSetV (stored in T[X]SetV).
  template <class T> void OutputSetV(string FileName, string Suffix, vector<set<T> > TwoDimSetV);
  // output a map (stored in T[X][X]Map)
  template <class T> void OutputMap(string FileName, string Suffix, map<T, T> FMap);
  // output a .gml file which can be fed into Cytoscape software for visualizing a graph.
  void OutputGML(string DataName, TIntSetV Neighbor, TIntSetV ClusInNode);
};

void IOCLASS::LoadIntVV(string DataName, string Suffix, TIntVV& TwoDimVV) {
  string FileName = DataName + "." + Suffix;
  const char* CharFileName = FileName.c_str();
  if (access(CharFileName, R_OK|W_OK) != 0) {
    printf("No IntVV file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finIntVV;
  finIntVV.open(CharFileName);
  string szLine;
  TIntV TempV;
  while(getline(finIntVV, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      TempV.push_back(atoi((*iter).c_str()));
    }
    TwoDimVV.push_back(TempV);
    TempV.clear();  
  }
  finIntVV.close();
  finIntVV.clear();
}

void IOCLASS::InputDblVV(string DataName, string Suffix, TDblVV& TwoDimVV) {
  string FileName = DataName + "." + Suffix;
  const char* CharFileName = FileName.c_str();
  if (access(CharFileName, R_OK|W_OK) != 0) {
    printf("No DblVV file exists, please check it!!!\n");
    exit(0);
  }
  ifstream finDblVV;
  finDblVV.open(CharFileName);
  string szLine;
  TDblV TempV;
  while(getline(finDblVV, szLine)) {
    TStrV tData;
    istringstream iss(szLine);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<TStrV >(tData));
    for (TStrV::iterator iter = tData.begin(); iter != tData.end(); ++iter) { 
      TempV.push_back(atof((*iter).c_str()));
    }
    TwoDimVV.push_back(TempV);
    TempV.clear();  
  }
  finDblVV.close();
  finDblVV.clear();
}

template <class T> 
void IOCLASS::OutputVec(string FileName, string Suffix, vector<T> OneDimVec) {
// output a vector.
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutV;
  foutV.open(File);
  for (int ind1 = 0; ind1 < OneDimVec.size(); ind1++) { foutV << OneDimVec[ind1] << endl; }
  foutV.close();
  foutV.clear();
}

template <class T> 
void IOCLASS::OutputSet(string FileName, string Suffix, set<T> OneDimSet) {
// output a vector.
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutSet;
  foutSet.open(File);
  typename std::set<T>::iterator TempIter;
  for (TempIter = OneDimSet.begin(); TempIter != OneDimSet.end(); TempIter++) {
    foutSet << *TempIter << " ";
  }
  foutSet.close();
  foutSet.clear();
}

template <class T> 
void IOCLASS::OutputVV(string FileName, string Suffix, vector<vector<T> > TwoDimVV) {
// output a TwoDimVV (stored in T[X]VV)
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutVV;
  foutVV.open(File);
  for (int ind1 = 0; ind1 < TwoDimVV.size(); ind1++) {
    for (int ind2 = 0; ind2 < TwoDimVV[ind1].size(); ind2++) { foutVV << TwoDimVV[ind1][ind2] << " "; }
    foutVV << endl;
  }
  foutVV.close();
  foutVV.clear();
}

template <class T> 
void IOCLASS::OutputSetV(string FileName, string Suffix, vector<set<T> > TwoDimSetV) {
// output a TwoDimSetV (stored in T[X]SetV)
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutSetV;
  foutSetV.open(File);
  typename std::set<T>::iterator TempIter;
  for (int ind1 = 0; ind1 < TwoDimSetV.size(); ind1++) {
    for (TempIter = TwoDimSetV[ind1].begin(); TempIter != TwoDimSetV[ind1].end(); TempIter++) {
      foutSetV << *TempIter << " ";
    }
  foutSetV << endl;
  }
  foutSetV.close();
  foutSetV.clear();
}

template <class T>
void IOCLASS::OutputMap(string FileName, string Suffix, map<T, T> FMap) {
  string FullName = FileName + "." + Suffix;
  const char* File = FullName.c_str();
  ofstream foutMap;
  foutMap.open(File);
  typename std::map<T, T>::iterator TempIter;
  for (TempIter = FMap.begin(); TempIter != FMap.end(); TempIter++) { 
    foutMap << TempIter->first << " " << TempIter->second << endl; 
  }
  foutMap.close();
  foutMap.clear();
}

void IOCLASS::OutputGML(string FileName, TIntSetV Neighbor, TIntSetV ClusInNode) {
// visualize the graph (output a .gml file which can be fed into software named Cytoscape
  string SGmlFile = FileName + ".gml";
  const char* GmlFile = SGmlFile.c_str();
  ofstream foutG;
  foutG.open(GmlFile);
  foutG << "graph" << endl;
  foutG << "[" << endl;
  foutG << "  directed 0" << endl;
  TIntSetIter it_set;
  int i;
  for (i = 0; i < Neighbor.size(); i++) {
    foutG << "  node" << endl;
    foutG << "  [" << endl;
    foutG << "    id " << i + 1 << endl;
    foutG << "    label " << "\"" << i + 1 << "\"" << endl;
    // non-overlapping clusters
    // if one node does not belong to any cluster, the flag equals 0;
    it_set = ClusInNode[i].begin();
    foutG << "    cluster " << "\"" << *it_set << "\"" << endl;
    foutG << "  ]" << endl;
  }
  
  for (i = 0; i < Neighbor.size(); i++) {
    for (it_set = Neighbor[i].begin(); it_set != Neighbor[i].end(); it_set++) {
      if (i < *it_set) {
        foutG << "  edge" << endl;
        foutG << "  [" << endl;
        foutG << "    source " << i + 1 << endl;
        foutG << "    target " << *it_set << endl;
        foutG << "  ]" << endl;
      }
    }
  }
  foutG << "]" << endl;
  foutG.close();
  foutG.clear();
}


//////////////////////////////////
///////////////////////////////////////////////
//Class of performance test on nodes clustering
/*Tips:
1. The cluster structure of calculating AvgF1, NMI, and ARI is the second Type (each line denotes a set of nodes in a common community)
2. The cluster structure of calculating Omega Index is the first Type (each line denotes a set of clusters a specific node belongs)
3. The number of nodes is needed when calculate NMI, Omega Index, and ARI.
*/
class ClusterTest {
public: 
  // no ground-truth information 
  // Modularity of non-overlapping cluster strucuture. (Newman, et al. PRE, 2004)
  // Our1 and Our2 are same cluster structure with different format.
  // each set in Our1 denotes a cluster, containing a set of nodes who belong it.
  // each set in Our2 is a set of clusters to which a specific node belongs.
  double CalNonOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges);
  // Modularity of overlapping cluster structure (Shen, et al. Physica A, 2009)
  double CalOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges);
  // Tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
  double CalTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor);
  // Adjusted tightness (overlapping/non-overlapping)(Bu, et al. Information Fusion, 2017)
  // penalizing very small and very large clusters and produces well-balanced solutions
  double CalAdjTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor);
  
public:
  // has ground-truth (.gt file) information
  // 1. Average of f1-score (AvgF1)
  double CalF1(TIntSet SetA, TIntSet SetB);
  double CalMaxF1(TIntSet SetA, TIntSetV TargetClus);
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalAvgF1(TIntSetV Our, TIntSetV GT);
  // 2. Normalized mutual information (NMI)
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalNMI(TIntSetV Our, TIntSetV GT, int num_nodes);
  // 3. Omega Index
  // each line in ClusInNode_Our and ClusInNode_GT is a set of clusters to which a specific node belongs
  double CalOmegaIndex(TIntSetV ClusInNode_Our, TIntSetV ClusInNode_GT, int num_nodes);
  // 4. Adjusted Rand Index (ARI)
  // each set in Our and GT denotes a cluster, containing a set of nodes who belong it
  double CalARI(TIntSetV Our, TIntSetV GT, int num_nodes);
  
  // calculate the number of elements in the intersection set of two sets
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
  TIntSet inline CalJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    return JointSet;
  }
};

double ClusterTest::CalNonOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
// last modified in 2021-01-13 by Yuyao Wang.
// calculate modularity of non-overlapping cluster structure. (Newman, PRE, 2004.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// if input is a overlapping cluster structure, this function will only take the first cluster for each node as her FuzzyMembership.
// num_edges is the number of all edges in the graph.
// the value of modularity lies between -1 and 1.
  if (neighbor.size() != Our2.size()) { 
    cout << "inappropriate inputs!!!" << endl;
    exit(1); 
  }
  TDblV theta;
  theta.resize(Our1.size());
  TIntSetIter it_set;
  for (int CID = 0; CID < theta.size(); CID++) {
    theta[CID] = 0.0;
    // total degree (in-degree + out-degree) of kth community, for the pre-computation of \sum_j d_j in the equation
    // here we sum all nodes including the target node who is removed sometimes in some papers.
    for (it_set = Our1[CID].begin(); it_set != Our1[CID].end(); it_set++) {
      theta[CID] += neighbor[*it_set - 1].size();
    }
    theta[CID] /= 2.0 * num_edges;
  }
  // modularity of single node
  TDblV singleModul;
  singleModul.resize(neighbor.size());
#pragma omp parallel for
  for (int NID = 1; NID <= neighbor.size(); NID++) {
    double tempValue = 0.0;
    singleModul[NID-1] = 0.0;
    // determine whether NID does not belong to any cluster.
    if (Our2[NID-1].empty()) { 
      singleModul[NID-1] = 0.0; 
    } else if (neighbor[NID-1].empty() == 0) {
      it_set = Our2[NID-1].begin(); // adopt the first value for non-overlapping communities (if there are multiple values)
      tempValue = CalNumJointSets(neighbor[NID-1], Our1[*it_set - 1]); // first term of equation
      singleModul[NID-1] = (tempValue - theta[*it_set - 1] * neighbor[NID-1].size()) / 2.0 / num_edges;
    }    
  }
  double SumModul = 0.0;
  // total value of modularity
  for (int NID = 1; NID <= neighbor.size(); NID++) { SumModul += singleModul[NID-1]; }
  return SumModul;
}

double ClusterTest::CalOverlapModul(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor, int num_edges) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating modularity of overlapping and hierarchical cluster structure.(Huawei Shen, et al. Eqn.(2) in Physica A, 2009.)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// num_edges is the number of all edges in the graph.
// here we tranform the defintion of community-oriented to node-oriented.
// In particular, we divide 1/O_iO_j.
  if (neighbor.size() != Our2.size()) { 
    cout << "inappropriate inputs!!!" << endl;
    exit(1); 
  }
  // modularity of single node
  TDblV singleModul;
  singleModul.resize(neighbor.size());
//#pragma omp parallel for
  for (int NID = 1; NID <= neighbor.size(); NID++) {
    double tempValue1 = 0.0, tempValue2 = 0.0;
    singleModul[NID-1] = 0.0;
    // determine whether NID does not belong to any cluster.
    if (Our2[NID-1].empty()) { 
      singleModul[NID-1] = 0.0; 
    } else if (neighbor[NID-1].empty() == 0) {
      for (TIntSetIter it_set = neighbor[NID-1].begin(); it_set != neighbor[NID-1].end(); it_set++) {
        if ((Our2[*it_set - 1].empty() == 0) && (CalNumJointSets(Our2[NID-1], Our2[*it_set - 1]) > 0)) {
          tempValue1 += (1.0 / Our2[*it_set - 1].size());
          tempValue2 += (neighbor[NID-1].size() * neighbor[*it_set - 1].size() / 2.0 / num_edges / Our2[*it_set - 1].size());
        }
      }
      singleModul[NID-1] = (tempValue1 - tempValue2) / Our2[NID-1].size();
    }    
  }
  double SumModul = 0.0;
  // total value of modularity
  for (int NID = 1; NID <= neighbor.size(); NID++) { SumModul += singleModul[NID-1]; }
  SumModul /= (2.0 * num_edges);
  return SumModul;
}

double ClusterTest::CalTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating tightness of cluster structure (Zhan Bu, et al. Information Fusion, 2017)
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
  double FirTerm = 0.0, SecTerm = 0.0;
  int num_clusters = Our1.size();
  int num_nodes = neighbor.size();
  double SumTgt = 0.0;
  TDblV SinglTgt;
  SinglTgt.resize(num_clusters);
  TIntV numInEdges;
  numInEdges.resize(num_clusters);
  TIntV numOutEdges;
  numOutEdges.resize(num_clusters);
  // a cluster k
#pragma omp parallel for
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      //id_i = *it_set;
      int TempNumJointSets = CalNumJointSets(neighbor[*it_set - 1], Our1[k]);
      numInEdges[k] += TempNumJointSets;
      numOutEdges[k] += (neighbor[*it_set - 1].size() - TempNumJointSets);
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
  }
  for (int k = 0; k < num_clusters; k++) {
    if (Our1[k].empty() == 0 && Our1[k].size() != num_nodes) {
      FirTerm = 2.0 * numInEdges[k] / Our1[k].size() / (double)Our1[k].size();
      SecTerm = numOutEdges[k] / Our1[k].size() / (double)(num_nodes - Our1[k].size());
      SinglTgt[k] = FirTerm - SecTerm;
    }
    SumTgt += SinglTgt[k];
  }
  return SumTgt;
}

double ClusterTest::CalAdjTgt(TIntSetV Our1, TIntSetV Our2, TIntSetV neighbor) {
// last modified in 2021-01-14 by Yuyao Wang.
// calculating adjusted tightness of cluster structure (Zhan Bu, et al. Information Fusion, 2017)
// penalizing very small and very large clusters and produces well-balanced solutions
// each set in Our1 denotes a cluster, containing a set of nodes who belongs it.
// each set in Our2 denotes a node, containing a set of clusters to which a specific node belongs.
// neighbor is the topology of graph.
  double FirTerm = 0.0, SecTerm = 0.0;
  int num_clusters = Our1.size();
  int num_nodes = neighbor.size();
  double SumTgt = 0.0;
  TDblV SinglTgt;
  SinglTgt.resize(num_clusters);
  TIntV numInEdges;
  numInEdges.resize(num_clusters);
  TIntV numOutEdges;
  numOutEdges.resize(num_clusters);
  // a cluster k
#pragma omp parallel for
  for (int k = 0; k < num_clusters; k++) {
    SinglTgt[k] = 0.0;
    numInEdges[k] = 0;
    numOutEdges[k] = 0;
    // a node id_i in the cluster-k
    for (TIntSetIter it_set = Our1[k].begin(); it_set != Our1[k].end(); it_set++) {
      //id_i = *it_set;
      int TempNumJointSets = CalNumJointSets(neighbor[*it_set - 1], Our1[k]);
      numInEdges[k] += TempNumJointSets;
      numOutEdges[k] += (neighbor[*it_set - 1].size() - TempNumJointSets);
    }
    // each linked node pair in cluster-k is counted twice.
    // but for the whole graph, each node pair who is not in common cluster is also counted twice,
    // so it is unecessary for the next command.
    // numInEdges[k] /= 2;
  }
  for (int k = 0; k < num_clusters; k++) {
    if (Our1[k].empty() == 0) {
      FirTerm = 2.0 * (num_nodes - Our1[k].size()) * numInEdges[k] /  (double)Our1[k].size();
      SinglTgt[k] = FirTerm - numOutEdges[k];
    }
    SumTgt += SinglTgt[k];
  }
  return SumTgt;
}

// calculate F1
double ClusterTest::CalF1(TIntSet SetA, TIntSet SetB) {
  int K_A = SetA.size();
  int K_B = SetB.size();
  int NumJoint = CalNumJointSets(SetA, SetB);
  double precision = 0.0, recall = 0.0;
  double f1 = 0.0;
  if (K_A != 0 && K_B != 0) {
    precision = (double)NumJoint / K_A;
    recall = (double)NumJoint / K_B;
  }
  if ((precision + recall) != 0 ) {
    f1 = 2.0 * precision * recall / (precision + recall);
    return f1;
  } else { return f1; }
}

// calculate maximum of F1
double ClusterTest::CalMaxF1(TIntSet SetA, TIntSetV TargetClus) {
  int K_A = SetA.size();
  int Size_T = TargetClus.size();
  double maxf1 = 0.0;
  for (int q = 0; q < Size_T; q++) {
    // double ClusterTest::CalF1(int* A, int* B)
    double tempf1 = CalF1(SetA, TargetClus[q]);
    if (tempf1 >= maxf1) { maxf1 = tempf1; }
  }
  return maxf1;
}

// calculating the average of f1-score
double ClusterTest::CalAvgF1(TIntSetV Our, TIntSetV GT) {
  int K_1 = Our.size();
  int K_2 = GT.size();
  TDblV MaxF1_1;
  MaxF1_1.resize(K_1);
  TDblV MaxF1_2;
  MaxF1_2.resize(K_2);
  int p = 0, q = 0;
  double sumMAXF1_1 = 0.0, sumMAXF1_2 = 0.0;
  for (p = 0; p < K_1; p++) {
    MaxF1_1[p] = CalMaxF1(Our[p], GT);
    sumMAXF1_1 += MaxF1_1[p];
  }
  for (q = 0; q < K_2; q++) {
    MaxF1_2[q] = CalMaxF1(GT[q], Our);
    sumMAXF1_2 += MaxF1_2[q];
  }
  double AvgF1 = sumMAXF1_1 / 2.0 / K_1 + sumMAXF1_2 / 2.0 / K_2;
  return AvgF1;
}

// calculating the normalized mutual information (NMI)
// [S. Fortunato. Community detection in graphs. Physics Reports, 486(3-5):75 – 174, 2010.]
double ClusterTest::CalNMI(TIntSetV Our, TIntSetV GT, int num_nodes) {
  int LenA = Our.size();
  int LenB = GT.size();
  int i = 1, j = 1;
  TDblVV number_joint;
  number_joint.resize(LenA);
  for (i = 0; i < LenA; i++) {
    number_joint[i].resize(LenB);
    for (j = 0; j < LenB; j++) {
      number_joint[i][j] = (double)CalNumJointSets(Our[i], GT[j]);
    }
  }
  double NMI = 0.0;
  TDblV NMI_A;
  NMI_A.resize(LenA);
  TDblV NMI_B;
  NMI_B.resize(LenB);
  for (i = 0; i < LenA; i++) { NMI_A[i] = 0.0; }
  for (j = 0; j < LenB; j++) { NMI_B[j] = 0.0; }
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) { NMI_A[i] += number_joint[i][j]; }
  }
  for (j = 0; j < LenB; j++) {
    for (i = 0; i < LenA; i++) { NMI_B[j] += number_joint[i][j]; }
  }
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) {
      if (number_joint[i][j] != 0.0) {
        NMI += number_joint[i][j] * log(number_joint[i][j] * num_nodes / NMI_A[i] / NMI_B[j]) / log(2);
      }
    }
  }
  NMI = NMI * (-1 * 2.0);
  double sum_NMI_A = 0.0;
  double sum_NMI_B = 0.0;
  for(i = 0; i < LenA; i++) {
    if (NMI_A[i] != 0.0) {
      sum_NMI_A += NMI_A[i] * log(NMI_A[i] / num_nodes) / log(2);
    }
  }
  for(j = 0; j < LenB; j++) {
    if(NMI_B[j] != 0) { sum_NMI_B += NMI_B[j] * log(NMI_B[j] / num_nodes) / log(2); }    
  }
  NMI = NMI / (sum_NMI_A + sum_NMI_B);
  return NMI;
}

// calculating the value of Omega index
// [S. Gregory. Fuzzy overlapping communities in networks. J. of Stat. Mech.: Theory and Experiment, 2011.]
// each set in ClusInNode_Our and ClusInNode_GT denotes a set of clusters to which a specific node belongs.
double ClusterTest::CalOmegaIndex(TIntSetV ClusInNode_Our, TIntSetV ClusInNode_GT, int num_nodes) {
  int SumTemp = 0;
  double results = 0.0;
  for (int i = 0; i < num_nodes; i++) {
    for (int j = i + 1; j < num_nodes; j++) {
      if (CalNumJointSets(ClusInNode_Our[i], ClusInNode_Our[j]) == CalNumJointSets(ClusInNode_GT[i], ClusInNode_GT[j])) { SumTemp++; }
    }
  }
  results = (double)SumTemp / num_nodes / num_nodes;
  return results;
}

// calculating Adjusted Rand Index (ARI)
// each set in Our1 and GT denotes a cluster, containing a set of nodes who belong it.
double ClusterTest::CalARI(TIntSetV Our, TIntSetV GT, int num_nodes) {
  int LenA = Our.size();
  int LenB = GT.size();
  int i, j;
  double numerator = 0.0, denominator = 0.0;
  double TempVal = 0.0;
  double numerator_1 = 0.0, numerator_2 = 0.0;
  for (i = 0; i < LenA; i++) {
    for (j = 0; j < LenB; j++) {
      TempVal = CalNumJointSets(Our[i], GT[j]);
      numerator_1 = numerator_1 + TempVal * (TempVal - 1.0) / 2.0;
    }
  }
  double TempVal_1 = 0.0;
  for (i = 0; i < LenA; i++) {
    TempVal_1 += Our[i].size() * (Our[i].size() - 1.0) / 2.0;
  }
  double TempVal_2 = 0.0;
  for (j = 0; j < LenB; j++) {
    TempVal_2 += GT[j].size() * (GT[j].size() - 1.0) / 2.0;
  }
  numerator_2 = (TempVal_1 + TempVal_2) * 2.0 / num_nodes / (num_nodes - 1);
  denominator = (TempVal_1 + TempVal_2) / 2.0 - numerator_2;
  double results = (numerator_1 - numerator_2) / denominator;
  return results;
}


///////////////////////////////////
///////////////////////////////////////////////
class CD4AT {
private:
  GRAPH G;
  int NumNodes, NumEdges; // number of nodes and edges
  int NumComs; // number of communities
  TIntV GtComSize; // size of communities, gt file
  TDblVV F; //* membership matrix (size: NumNodes * NumClus)
  TDblV SumF; // community size (Sum_i F_ic for community c)
  TDblVV PI; //* payoff matrix (size: NumClus * NumClus)
  TDblVV Theta; //* interaction matrix (size: NumClus * NumClus)
  TDblV Lambda; //* multiplier vector (size: 1 * NumClus)
  float alpha; // a positive constant
  float beta, gamma, zeta; // learning rate
  int MaxTotalIter; // maximum number of iterations
  
public:
  CD4AT(GRAPH graph, float Alpha, float stepsize, int maxiter): G(graph.GetDataName()), alpha(Alpha), gamma(stepsize), beta(stepsize), zeta(stepsize), MaxTotalIter(maxiter) {}
  void Initialization(); // load data
  void NeighborComInit(); // initialize F with best neighborhood communities
  void KeyParameterInit();
  
  double ComputeLikelihood();
  void UpdateF();
  void UpdateTheta();
  void UpdateLambda();
  void UpdateSinglePI(int NID1, int NID2);
  void UpdatePI();
  void TransF2Com(TIntSetV& TempClusInClus, TIntSetV& TempClusInNode);
  void PeformCD4AT(); // perform CD4AT
  
  void inline UpdateCom(int NID, int CID, double Val) {
    // update F and SumF
    SumF[CID-1] -= F[NID-1][CID-1];
    F[NID-1][CID-1] = Val;
    SumF[CID-1] += Val;
  }
  double inline DotProduct(TDblV Vec1, TDblV Vec2) {
    double ValDP = 0.0;
    if (Vec1.size() != Vec2.size()) {
      printf("Two Vectors have different size!!!\n");
      exit(0);
    } else {
      for (int VID = 0; VID < Vec1.size(); VID++) {
        ValDP += Vec1[VID] * Vec2[VID];
      }
    }
    return ValDP;
  }
  void inline NormVec(TDblV& Vec) {
    if(Vec.empty()) { 
      printf("Illegal for normalize a empty vector!!!\n"); 
      exit(1);
    }
    double SumVec = 0.0;
    for (int VID = 1; VID <= Vec.size(); VID++) { SumVec += Vec[VID-1]; }
    for (int VID = 1; VID <= Vec.size(); VID++) { 
      if (SumVec != 0) { Vec[VID-1] /= SumVec; }
    }
  }
};

void CD4AT::Initialization() {
  G.LoadGraphTpl(); // load .graph file (pure directed graph)
  NumNodes = G.GetNumNodes();
  NumEdges = G.GetNumEdges();
  printf("Number of Nodes: %d, Number of edges: %d\n", NumNodes, NumEdges);
  G.LoadClusGt(); // load .gt file (cluster ground-truth)
  NumComs = G.GetNumClus();
  printf("Number of clusters: %d\n", NumComs);
}

void CD4AT::NeighborComInit() {
  // initialize F with best neighborhood communities (Gleich et.al. KDD'12)
  int i, j;
  F.resize(NumNodes);
  for (j = 0; j < NumComs; j++) { SumF.push_back(0.0); }
  for (i = 0; i < F.size(); i++) {
    for (j = 0; j < NumComs; j++) { F[i].push_back(0.0); }
  }
  TDblIntPrV NIdPhiV;
  TIntSet InvalidNIDS;
  // compute conductance of neighborhood community
  CENTRALITY myCentrality;
  for (int NId = 1; NId <= F.size(); NId++) {
    double Phi = myCentrality.CondOneNode(G, NId);
    NIdPhiV.push_back(make_pair(Phi, NId));
  }
  sort(NIdPhiV.begin(), NIdPhiV.end());
  printf("Conductance computation completed!!!\n");
  // choose nodes with local minimum in conductance
  int CurCID = 1;
  for (int ui = 0; ui < NIdPhiV.size(); ui++) {
    int UID = NIdPhiV[ui].second;
    if (InvalidNIDS.find(UID) != InvalidNIDS.end()) { continue; }
    // add the node and its neighbors to the current community
    UpdateCom(UID, CurCID, 1.0);
    TIntSet Ngh(G.GetNeighbor(UID)); // G.GetNeighbor(UID).begin() is a bug. (abandoned)
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
      UpdateCom(*it_set, CurCID, 1.0);
      // exclude its neighbors from the next considerations
      InvalidNIDS.insert(*it_set);
    }
    CurCID++;
    if (CurCID > NumComs) { break; }
  }
  if (NumComs >= CurCID) {
    printf("%d communities needed to fill randomly\n", NumComs - CurCID + 1);
  }
  //assign a member to zero-member community (if any)
  for (int CID = 1; CID <= SumF.size(); CID++) {
    if (SumF[CID-1] == 0.0) {
      int ComSz = 10;
      for (int u = 0; u < ComSz; u++) {
        int UID = RndUniDevInt(NumNodes) + 1;
        UpdateCom(UID, CID, RndUniDev());
      }
    }
  }
  // normalize F
  for (int NID = 1; NID <= NumNodes; NID++) { NormVec(F[NID-1]); }
}

void CD4AT::KeyParameterInit() {
  // compute SumF
  SumF.resize(NumComs);
  for (int CID = 1; CID <= NumComs; CID++) {
    SumF[CID-1] = 0.0;
    for (int NID = 1; NID <= NumNodes; NID++) {
      SumF[CID-1] += F[NID-1][CID-1]; 
    }
  }
  srand((unsigned)time(0));
  int CID1, CID2; // community ID
  // initialize PI, Theta
  PI.resize(NumComs);
  Theta.resize(NumComs);
  for (CID1 = 1; CID1 <= NumComs; CID1++) { 
    for (CID2 = 1; CID2 <= NumComs; CID2++) { 
      PI[CID1-1].push_back(RndUniDev());
      Theta[CID1-1].push_back(RndUniDev()); 
    }
  }
  // initialize Lambda
  for (CID1 = 1; CID1 <= NumComs; CID1++) { 
    Lambda.push_back(RndUniDev()); 
    //cout << "node-" << CID << ": " << Lambda[CID-1] << endl;
  }
  TIntSetV ClusInClus = G.GetClusInClus();
  for (CID1 = 1; CID1 <= NumComs; CID1++) { GtComSize.push_back(ClusInClus[CID1-1].size()); }
}

double CD4AT::ComputeLikelihood() {
  double likelihood = 0.0;
  double TempVal1 = 0.0, TempVal2 = 0.0, TempVal3 = 0.0;
  // first term
  for (int NID = 1; NID <= NumNodes; NID++) {
    TDblV TempV;
    for (int CID = 1; CID <= NumComs; CID++) { TempV.push_back(DotProduct(F[NID-1], PI[CID-1])); }
    TIntSet Ngh(G.GetNeighbor(NID));
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) { likelihood += DotProduct(TempV, F[*it_set-1]); }
  }
  // second term
  for (int CID = 1; CID <= NumComs; CID++) {
    TempVal1 = SumF[CID-1] - GtComSize[CID-1];
    likelihood = likelihood + Lambda[CID-1] / NumNodes * TempVal1;
    likelihood = likelihood + alpha * 0.5 * pow((TempVal1 / NumNodes),2);
  }
  return likelihood;
}

void CD4AT::UpdateSinglePI(int NID1, int NID2) {
  for(int CID1 = 1; CID1 <= NumComs; CID1++) {
    for(int CID2 = 1; CID2 <= NumComs; CID2++) {
      if (CID1==CID2) {
        PI[CID1-1][CID2-1] += exp(F[NID1-1][CID1-1] * Theta[CID1-1][CID2-1] * F[NID2-1][CID1-1] - SumF[CID1-1]);
      } else {
        PI[CID1-1][CID2-1] += exp(F[NID1-1][CID1-1] * Theta[CID1-1][CID2-1] * F[NID2-1][CID1-1]);
      }
    }
  }
}

void CD4AT::UpdatePI() {
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet Ngh(G.GetNeighbor(NID));
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
      UpdateSinglePI(NID, *it_set);
    }
  }
  for(int CID1 = 1; CID1 <= NumComs; CID1++) {
    for(int CID2 = 1; CID2 <= NumComs; CID2++) {
      PI[CID1-1][CID2-1] /= (2.0 * NumEdges);
    }
  }
}

void CD4AT::UpdateF() {
  double PartialF = 0.0;
  double TempVal = 0.0, TempVal1 = 0.0, TempVal2 = 0.0;
  double UpdateVal = 0.0;
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet Ngh(G.GetNeighbor(NID));
    for (int CID1 = 1; CID1 <= NumComs; CID1++) {
      for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
        TempVal1 = F[NID-1][CID1-1] * Theta[CID1-1][CID1-1] * F[*it_set-1][CID1-1] -SumF[CID1-1];
        double FirTerm = F[*it_set-1][CID1-1] * pow(exp(TempVal1),2) * (Theta[CID1-1][CID1-1] * F[*it_set-1][CID1-1] - 1);
        PartialF += FirTerm; // add first
        for (int CID2 = 1; CID2 <= NumComs; CID2++) {
          TempVal = 0.0;
          if (CID2 != CID1) {
            TempVal2 = F[NID-1][CID1-1] * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1];
            double SecTerm = F[*it_set-1][CID2-1] * pow(exp(TempVal2),2) * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1];
            TempVal -= SecTerm;
          }
        } // end CID2 (l /neq k)
        PartialF += TempVal; // add second
      } // end i's neighbor
      PartialF = PartialF * 0.5 + (Lambda[CID1-1] - alpha * GtComSize[CID1-1]) / NumNodes + alpha * SumF[CID1-1] / NumNodes / NumNodes;
      UpdateVal = F[NID-1][CID1-1] + PartialF * beta;
      if (UpdateVal > 0) { F[NID-1][CID1-1] = UpdateVal; } // update
    } // end CID1
    NormVec(F[NID-1]); // normalization
  } // end NID
}

void CD4AT::UpdateTheta() {
  double PartialTheta = 0.0;
  double TempVal1 = 0.0, TempVal2 = 0.0;
  for (int CID1 = 1; CID1 <= NumComs; CID1++) { 
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      if (CID1==CID2) {
        for (int NID = 1; NID <= NumNodes; NID++) {
          TIntSet Ngh(G.GetNeighbor(NID));
          for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
            TempVal1 = F[NID-1][CID1-1] * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1] - SumF[CID1-1];
            PartialTheta = exp(TempVal1) * pow(F[NID-1][CID1-1],2) * pow(F[*it_set-1][CID2-1],2); // partial
            TempVal2 = Theta[CID1-1][CID2-1] + PartialTheta * beta;
            if (TempVal2 > 0) { Theta[CID1-1][CID2-1] = TempVal2; } // update
          }
        }
      } else {
        for (int NID = 1; NID <= NumNodes; NID++) {
          TIntSet Ngh(G.GetNeighbor(NID));
          for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
            TempVal1 = F[NID-1][CID1-1] * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1];
            PartialTheta = -exp(TempVal1) * pow(F[NID-1][CID1-1],2) * pow(F[*it_set-1][CID2-1],2); // partial
            TempVal2 = Theta[CID1-1][CID2-1] + PartialTheta * beta;
            if (TempVal2 > 0) { Theta[CID1-1][CID2-1] = TempVal2; } // update
          }
        }
      } // endif
    }
  }
}

void CD4AT::UpdateLambda() {
  double PartialLambda = 0.0;
  for (int CID = 1; CID <= NumComs; CID++) {
    PartialLambda = (SumF[CID-1] - GtComSize[CID-1]) / NumNodes;
    Lambda[CID-1] += (PartialLambda * zeta);
  }
}

void CD4AT::TransF2Com(TIntSetV& TempClusInClus, TIntSetV& TempClusInNode) {
  TempClusInClus.clear();
  TempClusInNode.clear();
  TempClusInClus.resize(NumComs);
  TempClusInNode.resize(NumNodes);
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      if (F[NID-1][CID-1] > 0.001) {
        TempClusInClus[CID-1].insert(F[NID-1][CID-1]);
        TempClusInNode[NID-1].insert(F[NID-1][CID-1]);
      }
    }
  }
}

void CD4AT::PeformCD4AT() {
  EVALGRAPH myEvalGraph;  
  Initialization();  
  NeighborComInit();
  KeyParameterInit();
  TDblV LikeliSeq;
  double OldLikelihood = ComputeLikelihood();
  UpdateF();
  UpdateTheta();
  UpdateLambda();
  UpdatePI();
  double NewLikelihood = ComputeLikelihood();
  int iter = 0;
  TIntSetV TempClusInClus, TempClusInNode;
  while ((abs(NewLikelihood - OldLikelihood) > 0.0001) && (iter < MaxTotalIter)) {
    LikeliSeq.push_back(NewLikelihood);
    OldLikelihood = NewLikelihood;
    TransF2Com(TempClusInClus, TempClusInNode);
    UpdatePI();
    UpdateF();
    UpdateTheta();
    UpdateLambda();
    NewLikelihood = ComputeLikelihood();
  }
}


//////////////////////////////////
///////////////////////////////////////////////
class DSCPCD {
private:
  GRAPH G;
  int NumNodes, NumEdges;
  int NumComs; // number of communities
  //TIntSetV Neighbor1; // topology of first tier graph
  
  int TypeSecOrdProximity; // type of calculation of implicit relationship
  double ThresGenSecOrdGraph; // [0,1], threshold for generating implicit relationship graph.
  
  TIntSetV Neighbor2; // topology of implicit relationship graph
  int NumEdges2; // number of edges of implicit relationship graph
  TIntVV NumJointNgh; // store the number of joint neighbors
  TIntVV NumUnionNgh; // store the number of union neighbors
  TDblV XI; // store the intensity of implicit relationship (numnodes * 1)
  
  TDblVV F, F1, F2; //* membership matrix (size: NumNodes * NumClus)
  TDblV SumF, SumF1, SumF2; // community size (Sum_i F_ic for community c)
  vector<vector<TDblVV > > WPI; // payoff matrices for each node pair
  TDblVV HPI; //* payoff matrix for all node pairs (size: NumClus * NumClus)
  TDblVV Theta; //* interaction matrix (size: NumClus * NumClus)
  TDblV Lambda1, Lambda2; //* multiplier vector (size: 1 * NumClus)

  float alpha; // a positive constant
  float beta, gamma, zeta; // learning rate
  int MaxTotalIter; // maximum number of iterations
  
  // for update
  TDblVV FHO, F1HO, F2HO; 
  TDblV SumFHO, SumF1HO, SumF2HO;
  vector<vector<TDblVV > > WPIHO;
  TDblVV HPIHO;
  TDblVV ThetaHO;
  TDblV Lambda1HO, Lambda2HO;
  
public:
  DSCPCD(GRAPH graph): G(graph.GetDataName()), TypeSecOrdProximity(1), ThresGenSecOrdGraph(0.0) {}
  void Initialization(); // load data
  void GetNumJointNgh(IOCLASS myIO);
  double NghBsdProx(int NID1, int NID2, int type);
  void GenerateSecondOrderGraph(TIntSetV& Neighbor2, int TypeSecOrdProximity, double ThresGenSecOrdGraph);
  void PlotMotivation();
  double CondOneNode(TIntSetV Neighbor, int NumEdges, int NID);
  void NeighborComInit();
  void BackupCurrentState();
  void ComputeWPI();
  void ComputeXI();
  void ComputeHPI();
  void ComputeNetState();
  double ComputeHappinessOneNode(int NID);
  double ComputeHappinessAllNodes();
  double ComputeParF(int NID, int CID, int TNet);
  void ComputeF();
  double ComputeParTheta(int CID);
  double ComputeParTheta(int CID1, int CID2);
  void ComputeTheta();
  double ComputeParLambda(int CID, int TNet); 
  void ComputeLambda();
  void RecoverLastState();
  void PeformDSCPCD();
  
  int inline CalNumJointSets(TIntSet FirSet, TIntSet SecSet) {
    TIntSet JointSet;
    set_intersection(FirSet.begin(), FirSet.end(), SecSet.begin(), SecSet.end(), insert_iterator<TIntSet >(JointSet, JointSet.begin()));
    int NumJoint = JointSet.size();
    return NumJoint;
  }
  void inline UpdateCom(TDblVV& TempF, TDblV& TempSumF, int NID, int CID, double Val) {
    // update F and SumF
    TempSumF[CID-1] -= TempF[NID-1][CID-1];
    TempF[NID-1][CID-1] = Val;
    TempSumF[CID-1] += Val;
  }
  void inline NormVec(TDblV& Vec) {
    if(Vec.empty()) { 
      printf("Illegal of normalizing a empty vector!!!\n"); 
      exit(1);
    }
    double SumVec = 0.0;
    for (int VID = 1; VID <= Vec.size(); VID++) { SumVec += Vec[VID-1]; }
    for (int VID = 1; VID <= Vec.size(); VID++) { 
      if (SumVec != 0) { Vec[VID-1] /= SumVec; }
    }
  }
  double inline DotProduct(TDblV Vec1, TDblV Vec2) {
    double ValDP = 0.0;
    if (Vec1.size() != Vec2.size()) {
      printf("Two Vectors have different size!!!\n");
      exit(0);
    } else {
      for (int VID = 0; VID < Vec1.size(); VID++) {
        ValDP += Vec1[VID] * Vec2[VID];
      }
    }
    return ValDP;
  }
};

void DSCPCD::Initialization() {
  G.LoadGraphTpl(); // load .graph file (pure directed graph)
  printf("Number of Nodes: %d, Number of edges: %d\n", G.GetNumNodes(), G.GetNumEdges());
  G.LoadClusGt(); // load .gt file (cluster ground-truth)
  printf("Number of clusters: %d\n", G.GetNumClus());
  
  F.resize(NumNodes);
  cout << "hello" << endl;
  F1.resize(NumNodes);
  F2.resize(NumNodes);
  for (int CID = 1; CID <= NumComs; CID++) { 
    SumF.push_back(0.0); 
    SumF1.push_back(0.0);
    SumF2.push_back(0.0);
  }
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) { 
       F[NID-1].push_back(0.0); 
       F1[NID-1].push_back(0.0);
       F2[NID-1].push_back(0.0);
    }
  }
  NeighborComInit();
  Theta.resize(NumComs);
  for (int CID1 = 1; CID1 <= NumComs; CID1++) { 
    Lambda1.push_back(RndUniDev());
    Lambda2.push_back(RndUniDev());
    for (int CID2 = 1; CID2 <= NumComs; CID2++) { 
      Theta[CID1-1].push_back(RndUniDev()); 
    }
  }
  WPI.resize(NumNodes);
  for (int NID = 1; NID <= NumNodes; NID++) {
    WPI[NID-1].resize(NumNodes);
  }
  HPI.resize(NumComs);
  for(int CID = 1; CID <= NumComs; CID++) {
    HPI[CID-1].resize(NumComs);
  }
  
  WPIHO.resize(NumNodes);
  for (int NID = 1; NID <= NumNodes; NID++) {
    WPIHO[NID-1].resize(NumNodes);
  }
  HPIHO.resize(NumComs);
  for(int CID = 1; CID <= NumComs; CID++) {
    HPIHO[CID-1].resize(NumComs);
  }
  FHO.resize(NumNodes);
  F1HO.resize(NumNodes);
  F2HO.resize(NumNodes);
  for (int NID = 1; NID <= NumNodes; NID++) {
    FHO[NID-1].resize(NumComs);
    F1HO[NID-1].resize(NumComs);
    F2HO[NID-1].resize(NumComs);
  }
  for (int CID = 1; CID <= NumComs; CID++) { 
    SumFHO.push_back(0.0); 
    SumF1HO.push_back(0.0);
    SumF2HO.push_back(0.0);
  }
  ThetaHO.resize(NumComs);
  for (int CID = 1; CID <= NumComs; CID++) {
    ThetaHO[CID-1].resize(NumComs);
  }
  Lambda1HO.resize(NumComs);
  Lambda2HO.resize(NumComs);
}

void DSCPCD::GetNumJointNgh(IOCLASS myIO) {
  // allocate space for NumJointNgh.
  NumJointNgh.resize(G.GetNumNodes());
  NumUnionNgh.resize(G.GetNumNodes());
  for(int NID = 1; NID <= NumJointNgh.size(); NID++) { 
    NumJointNgh[NID-1].resize(NumJointNgh.size()-NID); 
    NumUnionNgh[NID-1].resize(NumJointNgh.size()-NID); 

  }
  // calculate NumJointNgh
  int CompletedNodes = 0;
//#pragma omp parallel for
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      TIntSet TempSet1 = G.GetNeighbor(NID1);
      TIntSet TempSet2 = G.GetNeighbor(NID2);
      NumJointNgh[NID1-1][NID2-NID1-1] = CalNumJointSets(TempSet1, TempSet2);
      NumUnionNgh[NID1-1][NID2-NID1-1] = TempSet1.size() + TempSet2.size() - NumJointNgh[NID1-1][NID2-NID1-1];
    }
    //printf("Number of Completed Nodes for calculating NumJointNgh: %d\n", CompletedNodes++);
    //fflush(stdout);
  }
  cout << endl;
  cout << "Calculation of Number of Joint Neighbors is fininshed!!!" << endl; 
  myIO.OutputVV(G.GetDataName(), "NumJointNgh", NumJointNgh); // pre-compute and store it for following fast computation
}

double DSCPCD::NghBsdProx(int NID1, int NID2, int type) {
  TIntSet FirNode = G.GetNeighbor(NID1);
  TIntSet SecNode = G.GetNeighbor(NID2);
  int LenFirNode = FirNode.size();
  int LenSecNode = SecNode.size();
  if (LenFirNode == 0 || LenSecNode == 0) { return 0.0; }
  int NumJoint = NumJointNgh[NID1-1][NID2-NID1-1]; // pre-computed
  int NumUnion = LenFirNode + LenSecNode - NumJoint;
  NumJoint += 2;
  NumUnion += 2;
  LenFirNode++;
  LenSecNode++;
  // Jaccard Coefficient
  if (type == 1) { return (double)NumJoint / NumUnion; }
  // Salton Index
  else if (type == 2) { return (double)NumJoint / (sqrt(LenFirNode * LenSecNode)); }
  // Sorensen Index
  else if (type == 3) { return 2.0 * NumJoint / (LenFirNode + LenSecNode); }
  // Hub Promoted Index
  else if (type == 4) {
    if (LenFirNode < LenSecNode) { return 2.0 * NumJoint / LenFirNode; }
    else { return 2.0 * NumJoint / LenSecNode; }
  }
  // Hub Depressed Index
  else if (type == 5) {
    if (LenFirNode < LenSecNode) { return 2.0 * NumJoint / LenSecNode; }
    else { return 2.0 * NumJoint / LenFirNode; }
  }
  // Leicht-Holme-Newman Index
  else if (type == 6) { return 2.0 * NumJoint / LenFirNode / LenSecNode; }
  return 0.0;
}

void DSCPCD::GenerateSecondOrderGraph(TIntSetV& Neighbor2, int TypeSecOrdProximity, double ThresGenSecOrdGraph) {
  Neighbor2.resize(G.GetNumNodes()); // core dumped if the space is not pre-allocated.
  NumEdges2 = 0;  
  // allocate space for SecOrdProx
  TDblVV SecOrdProx;
  SecOrdProx.resize(G.GetNumNodes());
  for(int NID = 1; NID <= SecOrdProx.size(); NID++) { SecOrdProx[NID-1].resize(SecOrdProx.size()-NID); } 
#pragma omp parallel for
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      SecOrdProx[NID1-1][NID2-NID1-1] = NghBsdProx(NID1, NID2, TypeSecOrdProximity);
      //SecOrdProx[NID1-1][NID2-NID1-1] = NumJointNgh[NID1-1][NID2-NID1-1]; // set the number of joint number as second-order proximity.
    }
  }
  
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) {
    for (int NID2 = NID1+1; NID2 <= G.GetNumNodes(); NID2++) {
      if(SecOrdProx[NID1-1][NID2-NID1-1] > ThresGenSecOrdGraph) {
        // construct a link when two nodes share at least one neighbor.
        Neighbor2[NID1-1].insert(NID2); // this opearation can not be parallized
        Neighbor2[NID2-1].insert(NID1);
      }
    }
  }
  for (int NID1 = 1; NID1 < G.GetNumNodes(); NID1++) { NumEdges2 += Neighbor2[NID1-1].size(); }
  NumEdges2 = NumEdges2 / 2;
  cout << "Number of Edges in Second-Tier Graph: " << NumEdges2 << endl;
}

void DSCPCD::PlotMotivation() {
  // for motivation plot
  Initialization();
  
  IOCLASS myIO; // for output
  
  ClusterTest mytest; // for test
  TDblV StoreMetric; // store value of metrics
  
  double Modul1 = mytest.CalNonOverlapModul(G.GetClusInClus(), G.GetClusInNode(), G.GetNeighbor(), G.GetNumEdges());
  double Tgt1 = mytest.CalTgt(G.GetClusInClus(), G.GetClusInNode(), G.GetNeighbor());
  //StoreMetric[0].push_back(Modul1);
  //StoreMetric.push_back(AdjTgt1);
  
  myIO.LoadIntVV(G.GetDataName(), "NumJointNgh", NumJointNgh); // load NumJointNgh pre-computed by GetNumJointNgh(myIO);
  
  for (TypeSecOrdProximity = 4; TypeSecOrdProximity <= 4; TypeSecOrdProximity++) {
    if (TypeSecOrdProximity == 1) { printf("Second-order proximity: Jaccard Coefficient.\n"); }
      else if (TypeSecOrdProximity == 2) { printf("Second-order proximity: Salton Index.\n"); }
      else if (TypeSecOrdProximity == 3) { printf("Second-order proximity: Sorensen Index.\n");  }
      else if (TypeSecOrdProximity == 4) { printf("Second-order proximity: Hub Promoted Index.\n"); }
      else if (TypeSecOrdProximity == 5) { printf("Second-order proximity: Hub Depressed Index.\n"); }
      else if (TypeSecOrdProximity == 6) { printf("Second-order proximity: Leicht-Holme-Newman Index.\n"); }
    for (int IndThresGenSecOrdGraph = 0; IndThresGenSecOrdGraph <= 5; IndThresGenSecOrdGraph++) {
      ThresGenSecOrdGraph = IndThresGenSecOrdGraph / 2.0;
      TIntSetV TempNeighbor2;
      GenerateSecondOrderGraph(TempNeighbor2, TypeSecOrdProximity, IndThresGenSecOrdGraph);
      Neighbor2 = TempNeighbor2;
      double Modul2 = mytest.CalNonOverlapModul(G.GetClusInClus(), G.GetClusInNode(), TempNeighbor2, NumEdges2);
      double Tgt2 = mytest.CalTgt(G.GetClusInClus(), G.GetClusInNode(), TempNeighbor2);
    }
    cout << endl;
  }
  //myIO.OutputVec(G.GetDataName(), "MetricResults", StoreMetric);
  /*myIO.OutputGML(G.GetDataName(), G.GetNeighbor(), G.GetClusInNode());
  string DataName = G.GetDataName() + "2";
  myIO.OutputGML(DataName, Neighbor2, G.GetClusInNode());
  */
  
  EVALGRAPH myEvalGraph;
  myEvalGraph.GetDegreeDistribution(G);
  myEvalGraph.GetMaxDeg(G);
  myEvalGraph.GetMinDeg(G);
  myEvalGraph.GetAvgDeg(G);
}

double DSCPCD::CondOneNode(TIntSetV Neighbor, int NumEdges, int NID) {
  double ImptNode = 0.0; // small value indicates an important node.
  //TIntSetV Neighbor = G.GetNeighbor();
  TIntSet NghNId(Neighbor[NID - 1]);
  NghNId.insert(NID); // did in SNAP
  int LenNghs = NghNId.size();
  if (LenNghs < 5) {
    ImptNode = 1.0;
    return ImptNode;
  }
  int Edges2 = 2 * NumEdges;
  int Vol = 0,  Cut = 0;
  for (TIntSetIter it_set = NghNId.begin(); it_set != NghNId.end(); it_set++) {
    for (TIntSetIter it_set1 = Neighbor[*it_set - 1].begin(); it_set1 != Neighbor[*it_set - 1].end(); it_set1++) {
      // whether her neighbor's neighbor is also her neighbor.
      if (NghNId.find(*it_set1) == NghNId.end()) { Cut += 1; }
    }
    // Vol store the summation of degree of all nodes inside this set 
    Vol += Neighbor[*it_set - 1].size();
  }
  // get conductance
  if (Vol != Edges2) {
    if (2 * Vol > Edges2) { ImptNode = Cut / double (Edges2 - Vol); }
    else if (Vol == 0) { ImptNode = 0.0; }
    else { ImptNode = Cut / double(Vol); }
  } else {
    if (Vol == Edges2) { ImptNode = 1.0; }
  }
  return ImptNode;
}

void DSCPCD::NeighborComInit() {
  TDblIntPrV NIdPhiV1, NIdPhiV2;
  TIntSet InvalidNIDS1, InvalidNIDS2;
  // compute conductance of neighborhood community
  for (int NId = 1; NId <= F.size(); NId++) {
    double Phi1 = CondOneNode(G.GetNeighbor(), G.GetNumEdges(), NId);
    NIdPhiV1.push_back(make_pair(Phi1, NId));
    double Phi2 = CondOneNode(Neighbor2, NumEdges2, NId);
    NIdPhiV2.push_back(make_pair(Phi2, NId));
  }
  sort(NIdPhiV1.begin(), NIdPhiV1.end());
  sort(NIdPhiV2.begin(), NIdPhiV2.end());
  printf("Conductance computation completed!!!\n");

  // choose nodes with local minimum in conductance
  int CurCID1 = 1;
  for (int ui = 0; ui < NIdPhiV1.size(); ui++) {
    int UID = NIdPhiV1[ui].second;
    if (InvalidNIDS1.find(UID) != InvalidNIDS1.end()) { continue; }
    // add the node and its neighbors to the current community
    UpdateCom(F1, SumF1, UID, CurCID1, 1.0);
    TIntSet Ngh(G.GetNeighbor(UID)); // G.GetNeighbor(UID).begin() is a bug. (abandoned)
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
      UpdateCom(F1, SumF1, *it_set, CurCID1, 1.0);
      // exclude its neighbors from the next considerations
      InvalidNIDS1.insert(*it_set);
    }
    CurCID1++;
    if (CurCID1 > NumComs) { break; }
  }
  if (NumComs >= CurCID1) {
    printf("%d communities needed to fill randomly\n", NumComs - CurCID1 + 1);
  }
  //assign a member to zero-member community (if any)
  for (int CID = 1; CID <= SumF1.size(); CID++) {
    if (SumF1[CID-1] == 0.0) {
      int ComSz = 10;
      for (int u = 0; u < ComSz; u++) {
        int UID = RndUniDevInt(NumNodes) + 1;
        UpdateCom(F1, SumF1, UID, CID, RndUniDev());
      }
    }
  }
  // normalize F1
  for (int NID = 1; NID <= NumNodes; NID++) { NormVec(F1[NID-1]); }
    
   int CurCID2 = 1;
  for (int ui = 0; ui < NIdPhiV2.size(); ui++) {
    int UID = NIdPhiV2[ui].second;
    if (InvalidNIDS2.find(UID) != InvalidNIDS2.end()) { continue; }
    // add the node and its neighbors to the current community
    UpdateCom(F2, SumF2, UID, CurCID2, 1.0);
    TIntSet Ngh(Neighbor2[UID-1]); // G.GetNeighbor(UID).begin() is a bug. (abandoned)
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
      UpdateCom(F2, SumF2, *it_set, CurCID2, 1.0);
      // exclude its neighbors from the next considerations
      InvalidNIDS2.insert(*it_set);
    }
    CurCID2++;
    if (CurCID2 > NumComs) { break; }
  }
  if (NumComs >= CurCID2) {
    printf("%d communities needed to fill randomly\n", NumComs - CurCID2 + 1);
  }
  //assign a member to zero-member community (if any)
  for (int CID = 1; CID <= SumF2.size(); CID++) {
    if (SumF2[CID-1] == 0.0) {
      int ComSz = 10;
      for (int u = 0; u < ComSz; u++) {
        int UID = RndUniDevInt(NumNodes) + 1;
        UpdateCom(F2, SumF2, UID, CID, RndUniDev());
      }
    }
  }
  // normalize F2
  for (int NID = 1; NID <= NumNodes; NID++) { NormVec(F2[NID-1]); }
  
  // generate F
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      if (F1[NID-1][CID-1] >= F2[NID-1][CID-1]) {
        F[NID-1][CID-1] = F1[NID-1][CID-1];
      }
      else {
        F[NID-1][CID-1] = F2[NID-1][CID-1];
      }
    }
  }
  // normalize F
  for (int NID = 1; NID <= NumNodes; NID++) { NormVec(F[NID-1]); }
  // compute SumF
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      UpdateCom(F, SumF, NID, CID, F[NID-1][CID-1]);
    }
  } 
}

void DSCPCD::BackupCurrentState() {
  for (int NID1 = 1; NID1 <= NumNodes; NID1++) {
    for (int NID2 = NID1; NID2 <= NumNodes; NID2++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        for (int CID2 = 1; CID2 <= NumComs; CID2++) {
          WPIHO[NID1-1][NID2-1] = WPI[NID1-1][NID2-1];
        }
      }
    }
  }
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      HPIHO[CID1-1][CID2-1] = HPI[CID1-1][CID2-1];
    }
  }
  for (int CID = 1; CID <= NumComs; CID++) {
    SumF1HO[CID-1] = SumF1[CID-1];
    SumF2HO[CID-1] = SumF2[CID-1];
  }
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      F1HO[NID-1][CID-1] = F1[NID-1][CID-1];
      F2HO[NID-1][CID-1] = F2[NID-1][CID-1];
      FHO[NID-1][CID-1] = F[NID-1][CID-1];
    }
  }
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      ThetaHO[CID1-1][CID2-1] = Theta[CID1-1][CID2-1];
    }
  }
  for (int CID = 1; CID <= NumComs; CID++) {
    Lambda1HO[CID-1] = Lambda1[CID-1];
    Lambda2HO[CID-1] = Lambda2[CID-1];
  }
}

void DSCPCD::ComputeWPI() {
  TDblVV TempPayoffMatrix;
  TempPayoffMatrix.resize(NumComs);
  for (int CID = 1; CID <= NumComs; CID++) {
    TempPayoffMatrix[CID-1].resize(NumComs);
  }
  for (int NID1 = 1; NID1 <= NumNodes; NID1++) {
    for (int NID2 = NID1; NID2 <= NumNodes; NID2++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        for (int CID2 = 1; CID2 <= NumComs; CID2++) {
          if (CID1==CID2) {
            TempPayoffMatrix[CID1-1][CID2-1] = exp(-SumF[CID1-1]) * exp(F[NID1-1][CID1-1]*Theta[CID1-1][CID2-1]*F[NID2-1][CID2-1]);
          } else {
            TempPayoffMatrix[CID1-1][CID2-1] = -exp(F[NID1-1][CID1-1]*Theta[CID1-1][CID2-1]*F[NID2-1][CID2-1]);
          }
        }
      }
      WPI[NID1-1][NID2-1] = TempPayoffMatrix;
    }
  }
}

void DSCPCD::ComputeXI() {
  XI.resize(NumNodes);
  for (int NID = 1; NID <= NumNodes; NID++) {
    XI[NID-1] = 0.0;
  }
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet TempNghSet = G.GetNeighbor(NID);
    for (TIntSetIter it_set = TempNghSet.begin(); it_set != TempNghSet.end(); it_set++) {
      double xi_ij = 0.0;
      if (NumUnionNgh[NID-1][*it_set-NID-1] != 0) {
        xi_ij = NumJointNgh[NID-1][*it_set-NID-1] / NumUnionNgh[NID-1][*it_set-NID-1];
      } else {
        xi_ij = 0.0;
      }
      XI[NID-1] += xi_ij;
    }
  }
}

void DSCPCD::ComputeHPI() {
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet Ngh(G.GetNeighbor(NID));
    for (TIntSetIter it_set = Ngh.begin(); it_set != Ngh.end(); it_set++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        for (int CID2 = 1; CID2 <= NumComs; CID2++) {
          double TempPI = WPI[NID-1][*it_set-1][CID1-1][CID2-1];
          double fraction = (NumJointNgh[NID-1][*it_set-NID-1] / NumUnionNgh[NID-1][*it_set-NID-1]) / XI[NID-1];
          HPI[CID1-1][CID2-1] += (fraction * TempPI);
        }
      }
    }
  }
}

void DSCPCD::ComputeNetState() {
  TDblV TempVec;
  TempVec.resize(NumComs);
  for (int CID2 = 1; CID2 <= NumComs; CID2++) {
    double TempVal = 0.0;
    for (int CID3 = 1; CID3 <= NumComs; CID3++) {
      TempVal += (SumF1[CID3-1] * HPI[CID3-1][CID2-1]);
    }
    TempVec[CID2-1] = TempVal;
  } 
  // second term
  double secTerm = DotProduct(TempVec, SumF1);
 
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    SumF1[CID1-1] = SumF1[CID1-1] + SumF1[CID1-1] * (DotProduct(HPI[CID1-1], SumF1) - secTerm);
  }
  
  for (int CID2 = 1; CID2 <= NumComs; CID2++) {
    double TempVal = 0.0;
    for (int CID3 = 1; CID3 <= NumComs; CID3++) {
      TempVal += (SumF2[CID3-1] * HPI[CID3-1][CID2-1]);
    }
    TempVec[CID2-1] = TempVal;
  } 
  // second term
  secTerm = DotProduct(TempVec, SumF2);
 
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    SumF2[CID1-1] = SumF2[CID1-1] + SumF2[CID1-1] * (DotProduct(HPI[CID1-1], SumF2) - secTerm);
  }
}

double DSCPCD::ComputeHappinessOneNode(int NID) {
  double H, H_ex, H_im, H_stru_cons;
  TDblV TempVec;
  TempVec.resize(NumComs);
  // explicit network
  H_ex = 0.0;
  TIntSet TempNghSet1 = G.GetNeighbor(NID);
  for (TIntSetIter it_set = TempNghSet1.begin(); it_set != TempNghSet1.end(); it_set++) {
    TDblVV TempPayoffMatrix1 = WPI[NID-1][*it_set-1];
    for (int CID = 1; CID <= NumComs; CID++) {
      double tempVal = 0.0;
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        tempVal += (F1[NID-1][CID1-1] * TempPayoffMatrix1[CID1-1][CID-1]);
      }
      TempVec[CID-1] = tempVal;
    }
    H_ex += DotProduct(TempVec, F1[*it_set-1]);
  }
  // implicit network
  H_im = 0.0;
  TIntSet TempNghSet2 = Neighbor2[NID-1];
  for (TIntSetIter it_set = TempNghSet2.begin(); it_set != TempNghSet2.end(); it_set++) {
    TDblVV TempPayoffMatrix2 = WPI[NID-1][*it_set-1];
    for (int CID = 1; CID <= NumComs; CID++) {
      double tempVal = 0.0;
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        tempVal += (F2[NID-1][CID1-1] * TempPayoffMatrix2[CID1-1][CID-1]);
      }
      TempVec[CID-1] = tempVal;
    }
    H_im += DotProduct(TempVec, F2[*it_set-1]);
  }
  // 2nd term
  H_stru_cons = 0.0;
  double TempVal1 = 0.0;
  for (int CID = 1; CID <= NumComs; CID++) {
    TempVal1 = F1[NID-1][CID-1] - F2[NID-1][CID-1];
    H_stru_cons += (TempVal1 * TempVal1);
  }
  H_stru_cons = sqrt(H_stru_cons);
  
  H = H_ex + H_im - H_stru_cons;
  return H;
}

double DSCPCD::ComputeHappinessAllNodes() {
  double TotalH = 0.0;
  for (int NID = 1; NID <= NumNodes; NID++) {
    TotalH += ComputeHappinessOneNode(NID);
  }
  return TotalH;
}

double DSCPCD::ComputeParF(int NID, int CID, int TNet) {
  double ParHF = 0.0;
  double TempVal1 = F1[NID-1][CID-1] + F1[NID-1][CID-1];
  double TempVal2 = F1[NID-1][CID-1] - F1[NID-1][CID-1];
  if (TNet==1) {
    if (TempVal1==0) {
      if (TempVal2==0) {
        ParHF = 0.0;
      } else {
        ParHF = 0.5 * abs(TempVal2) / TempVal2;
      }
    } else {
      if (TempVal2==0) {
        ParHF = 0.5 * abs(TempVal1) / TempVal1;
      } else {
        ParHF = 0.5 * (abs(TempVal1) / TempVal1 + abs(TempVal2) / TempVal2);
      }
    }
  }
  if (TNet==2) {
    if (TempVal1==0) {
      if (TempVal2==0) {
        ParHF = 0.0;
      } else {
        ParHF = -0.5 * abs(TempVal2) / TempVal2;
      }
    } else {
      if (TempVal2==0) {
        ParHF = 0.5 * abs(TempVal1) / TempVal1;
      } else {
        ParHF = 0.5 * (abs(TempVal1) / TempVal1 - abs(TempVal2) / TempVal2);
      }
    }
  }
  
  double Delta_1 = 0.0;
  if (TNet==1) { 
    TIntSet TempNghSet1 = G.GetNeighbor(NID);
    for (TIntSetIter it_set = TempNghSet1.begin(); it_set != TempNghSet1.end(); it_set++) {
      double TempVal3 = exp(F[NID-1][CID-1]*Theta[CID-1][CID-1]*F[*it_set-1][CID-1]-SumF[CID-1]);
      Delta_1 += F1[*it_set-1][CID-1] * TempVal3 * TempVal3 * (Theta[CID-1][CID-1] * F[*it_set-1][CID-1] - 1) * ParHF;
    }
  }
  if (TNet==2) { 
    TIntSet TempNghSet2 = Neighbor2[NID-1];
    for (TIntSetIter it_set = TempNghSet2.begin(); it_set != TempNghSet2.end(); it_set++) {
      double TempVal4 = exp(F[NID-1][CID-1]*Theta[CID-1][CID-1]*F[*it_set-1][CID-1]-SumF[CID-1]);
      Delta_1 += F2[*it_set-1][CID-1] * TempVal4 * TempVal4 * (Theta[CID-1][CID-1] * F[*it_set-1][CID-1] - 1) * ParHF;
    }
  }
  
  double Delta_2 = 0.0;
  if (TNet==1) {
    TIntSet TempNghSet1 = G.GetNeighbor(NID);
    for (TIntSetIter it_set = TempNghSet1.begin(); it_set != TempNghSet1.end(); it_set++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        if (CID1==CID) { 
          continue;
        }
        double TempVal5 = exp(F[NID-1][CID-1]*Theta[CID-1][CID1-1]*F[*it_set-1][CID1-1]);
        Delta_2 += (-F1[*it_set-1][CID1-1]) * TempVal5 * TempVal5 * Theta[CID-1][CID1-1] * F[*it_set-1][CID1-1] * ParHF;
      }
    }
  }
  if (TNet==2) {
    TIntSet TempNghSet2 = Neighbor2[NID-1];
    for (TIntSetIter it_set = TempNghSet2.begin(); it_set != TempNghSet2.end(); it_set++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        if (CID1==CID) { 
          continue;
        }
        double TempVal6 = exp(F[NID-1][CID-1]*Theta[CID-1][CID1-1]*F[*it_set-1][CID1-1]);
        Delta_2 += (-F2[*it_set-1][CID1-1]) * TempVal6 * TempVal6 * Theta[CID-1][CID1-1] * F[*it_set-1][CID1-1] * ParHF;
      }
    }
  }
  
  double Term_3 = 0.0;
  if (TNet==1) {
    Term_3 = 2.0 * (F1[NID-1][CID-1] - F2[NID-1][CID-1]);
  }
  if (TNet==2) {
    Term_3 = 2.0 * (F2[NID-1][CID-1] - F1[NID-1][CID-1]);
  }
  
  double Term_4 = 0.0;
  if (TNet==1) {
    Term_4 = (1 / NumNodes) * (Lambda1[CID-1] + alpha * (SumF1[CID-1] / NumNodes - SumF1[CID-1]));
  }
  if (TNet==2) {
    Term_4 = (1 / NumNodes) * (Lambda2[CID-1] + alpha * (SumF2[CID-1] / NumNodes - SumF2[CID-1]));
  }
  
  double ParF = Delta_1 + Delta_2 + Term_3 + Term_4;
  return ParF;
}

void DSCPCD::ComputeF() {
  // for explicit network
  int TNet = 1;
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      double ParF = ComputeParF(NID, CID, TNet);
      F1[NID-1][CID-1] = F1[NID-1][CID-1] + beta * ParF;
      if (F1[NID-1][CID-1] < 0) {
        F1[NID-1][CID-1] = 0;
      }
    }
  }
  // for implicit network
  TNet = 2;
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      double ParF = ComputeParF(NID, CID, TNet);
      F2[NID-1][CID-1] = F2[NID-1][CID-1] + beta * ParF;
      if (F2[NID-1][CID-1] < 0) {
        F2[NID-1][CID-1] = 0;
      }
    }
  }
    // update F
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      if (F1[NID-1][CID-1] >= F2[NID-1][CID-1]) {
        F[NID-1][CID-1] = F1[NID-1][CID-1];
      }
      else {
        F[NID-1][CID-1] = F2[NID-1][CID-1];
      }
    }
  }
}

double DSCPCD::ComputeParTheta(int CID) {
  double ParTheta = 0.0;
  // (TNet==1)
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet TempNghSet1 = G.GetNeighbor(NID);
    for (TIntSetIter it_set = TempNghSet1.begin(); it_set != TempNghSet1.end(); it_set++) {
      ParTheta += (exp(F[NID-1][CID-1] * Theta[CID-1][CID-1] * F[*it_set-1][CID-1] - SumF[CID-1]) * F[NID-1][CID-1] * F[*it_set-1][CID-1] * F1[NID-1][CID-1] * F1[*it_set-1][CID-1]);
    }
  }
  // (TNet==2)
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet TempNghSet2 = Neighbor2[NID-1];
    for (TIntSetIter it_set = TempNghSet2.begin(); it_set != TempNghSet2.end(); it_set++) {
      ParTheta += (exp(F[NID-1][CID-1] * Theta[CID-1][CID-1] * F[*it_set-1][CID-1] - SumF[CID-1]) * F[NID-1][CID-1] * F[*it_set-1][CID-1] * F2[NID-1][CID-1] * F2[*it_set-1][CID-1]);
    }
  }
  return ParTheta;
}

double DSCPCD::ComputeParTheta(int CID1, int CID2) {
  double ParTheta = 0.0;
  // (TNet==1)
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet TempNghSet1 = G.GetNeighbor(NID);
    for (TIntSetIter it_set = TempNghSet1.begin(); it_set != TempNghSet1.end(); it_set++) {
      ParTheta += (-exp(F[NID-1][CID1-1] * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1]) * F[NID-1][CID1-1] * F[*it_set-1][CID2-1] * F1[NID-1][CID1-1] * F1[*it_set-1][CID2-1]);
    }
  }
  // (TNet==2)
  for (int NID = 1; NID <= NumNodes; NID++) {
    TIntSet TempNghSet2 = Neighbor2[NID-1];
    for (TIntSetIter it_set = TempNghSet2.begin(); it_set != TempNghSet2.end(); it_set++) {
      ParTheta += (exp(F[NID-1][CID1-1] * Theta[CID1-1][CID2-1] * F[*it_set-1][CID2-1]) * F[NID-1][CID1-1] * F[*it_set-1][CID2-1] * F2[NID-1][CID1-1] * F2[*it_set-1][CID2-1]);
    }
  }
  return ParTheta;
}

void DSCPCD::ComputeTheta() {
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      if (CID1==CID2) {
        Theta[CID1-1][CID2-1] = Theta[CID1-1][CID2-1] + gamma * ComputeParTheta(CID1);
      } else {
        Theta[CID1-1][CID2-1] = Theta[CID1-1][CID2-1] + gamma * ComputeParTheta(CID1, CID2);
      }
      if (Theta[CID1-1][CID2-1] < 0) {
        Theta[CID1-1][CID2-1] = 0;
      }
    }
  }
}

double DSCPCD::ComputeParLambda(int CID, int TNet) {
  double ParLambda = 0.0;
  if (TNet==1) {
   ParLambda = SumF1[CID-1] / NumNodes - SumF1[CID-1];
  }
  if (TNet==2) {
   ParLambda = SumF2[CID-1] / NumNodes - SumF2[CID-1];
  }
  return ParLambda;
}

void DSCPCD::ComputeLambda() {
  // for explicit network
  int TNet = 1;
  for (int CID = 1; CID <= NumComs; CID++) {
    Lambda1[CID-1] = Lambda1[CID-1] + ComputeParLambda(CID, TNet) * zeta;
  }
  // for implicit network
  TNet = 2;
  for (int CID = 1; CID <= NumComs; CID++) {
    Lambda2[CID-1] = Lambda2[CID-1] + ComputeParLambda(CID, TNet) * zeta;
  }
}

void DSCPCD::RecoverLastState() {
  // recover WPI
  for (int NID1 = 1; NID1 <= NumNodes; NID1++) {
    for (int NID2 = NID1; NID2 <= NumNodes; NID2++) {
      for (int CID1 = 1; CID1 <= NumComs; CID1++) {
        for (int CID2 = 1; CID2 <= NumComs; CID2++) {
          WPI[NID1-1][NID2-1] = WPIHO[NID1-1][NID2-1];
        }
      }
    }
  }
  // recover HPI
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      HPI[CID1-1][CID2-1] = HPIHO[CID1-1][CID2-1];
    }
  }
  // recover net state
  for (int CID = 1; CID <= NumComs; CID++) {
    SumF1[CID-1] = SumF1HO[CID-1];
    SumF2[CID-1] = SumF2HO[CID-1];
  }
  // recover F
  for (int NID = 1; NID <= NumNodes; NID++) {
    for (int CID = 1; CID <= NumComs; CID++) {
      F1[NID-1][CID-1] = F1HO[NID-1][CID-1];
      F2[NID-1][CID-1] = F2HO[NID-1][CID-1];
      F[NID-1][CID-1] = FHO[NID-1][CID-1];
    }
  }
  // recover theta
  for (int CID1 = 1; CID1 <= NumComs; CID1++) {
    for (int CID2 = 1; CID2 <= NumComs; CID2++) {
      Theta[CID1-1][CID2-1] = ThetaHO[CID1-1][CID2-1];
    }
  }
  // recover lambda
  for (int CID = 1; CID <= NumComs; CID++) {
    Lambda1[CID-1] = Lambda1HO[CID-1];
    Lambda2[CID-1] = Lambda2HO[CID-1];
  }
}

void DSCPCD::PeformDSCPCD() {
  Initialization();
  IOCLASS myIO;
  GetNumJointNgh(myIO);
  GenerateSecondOrderGraph(Neighbor2, TypeSecOrdProximity, ThresGenSecOrdGraph);
  ComputeWPI();
  ComputeXI();
  ComputeHPI();
  ComputeNetState();
  double H, TempH;
  int NumIter = 0;
  while (NumIter <= 200) {
    H = ComputeHappinessAllNodes();
    BackupCurrentState();
    ComputeWPI();
    ComputeHPI();
    ComputeNetState();
    ComputeF();
    ComputeTheta();
    ComputeLambda();
    TempH = ComputeHappinessAllNodes();
    if (TempH >= H) {
      continue;
    } else {
      RecoverLastState();
    }
    NumIter++;
  }
}

int main(int argc, char **argv)
{
  string head = argv[1]; // name of dataset
  //float Alpha = atof(argv[2]); 
  //float stepsize = atof(argv[3]);
  //int maxiter = atoi(argv[4]);
  //float delta = atof(argv[5]);
  
  GRAPH myGraph(head); // initialize an object
  struct timeval tod1, tod2; // record time
  gettimeofday(&tod1, NULL);
  //CD4AT CD4AT(myGraph, Alpha, stepsize, maxiter);
  DSCPCD DSCPCD(myGraph);
  DSCPCD.PeformDSCPCD();
  gettimeofday(&tod2, NULL);
  cout << "========= " << "Time cost: " << todiff(&tod2, &tod1) / 1000000.0 << "s" << " =========" << endl;
  cout << endl;

  return 0;
}
