# Dual Structural Consistency Preserving Community Detection on Social Networks (DSCPCD)

# Abstract
Community detection on social networks is a fundamental and crucial research topic in the field of social computing. Here we propose \emph{SCPCDG}---a structural consistency preserving community detection game model for community detection, which is designed regarding two criteria: 1) users within a social network should interact with each other in a manner combing uncertainty and certainty; 2) for a given social network, original explicit network (two linked users are friends) and potential implicit network (two linked users have common friends) should have a consistent community structure. In particular, \emph{SCPCDG} formulates each user in a social network as an individual in an evolutionary game associated with community-aware payoff settings, where the community sate evolves under the guide of replicator dynamics. The adoption of evolutionary game makes it possible to introduce the uncertainty to \emph{SCPCDG}. To further seek each user's membership, we develop a \emph{happiness} index to measure all user's satisfaction towards two community structures in explicit and implicit networks respectively, meanwhile the consistency of community structure in two networks is also characterized. Specifically, each user is assumed to maximize the \emph{happiness} under the constrains of evolutionary community state. This \emph{happiness} maximization operation is also based on the users' rationality, i.e., certainty. We evaluate \emph{SCPCDG} on several real-world and synthetic datasets, and the results show that it can yield substantial performance gains in terms of detection accuracy over several baselines. Also, we carry out an ablation test cutting the structural consistency criterion from \emph{happiness}, and a case study on a \emph{polbook} network to further verify \emph{SCPCDG}'s performance.
  
# Commands:
g++ dscpcd.cpp -o dscpcd -fopenmp -std=c++11
./dscpcd DataName

# 1. Format of input graph topological files.  
SONTA supports two format of input topology files. All nodes will be renamed and the correspondence between nodes' new id and original name is stored in a map<string, int> structure.

1) edge list (.edgelist).  
NodeName1 NodeName2  
NodeName3 NodeName4  
...  ...  
...  ...  
...  ...  

2) adjacency list (.adjlist). (first name is target node)  
TargetNodeName1 NeighborNodeName1 NeighborNodeName2 NeighborNodeName3 ...... NeighborNodeNameI  
TargetNodeName2 NeighborNodeName4  
TargetNodeName3 NeighborNodeName5 NeighborNodeName6 NeighborNodeName7 ...... NeighborNodeNameJ  
......  
......  
......  

# 2. Format of input graph cluster ground truth files.  
SONTA supports two format of input cluster ground truth files. Similarly, all clusters will be renamed and the correspondence between clusters' new id and their original name is stored in a map<string, int> structure.

1) each line denotes a node's affilication (.ngt). (first name is target node, others are cluster name)  
NodeName1 ClusterName1 ClusterName2 ...... ClusterNameI  
NodeName2 ClusterName3  
NodeName3 ClusterName4 ClusterName5 ...... ClusterNameJ  
......  
......  
......  

2) edge line denotes a certain cluster to which a set of nodes belong (.cgt). (first name is target cluster, others are node name)  
ClusterName1 NodeName1 NodeName2 ...... NodeNameI  
ClusterName2 NodeName3  
ClusterName3 NodeName4 NodeName5 ...... NodeNameJ  
......  
......  
......  
