#include <iostream>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <cstdlib>
#include <ctime> 
#include "tree.h"

using namespace std;

class Container{
	 public:
         int nInstances;
         int nVariables;
	 variable **variables;
         double **data;
         int *weights;
         int* elitismList;
         int nTrees;
         int minbucket;
         int minsplit;
         int maxNode; // maximum number of nodes a tree with maxdepth level can have
         int maxCat;  // maximum number of categories of nominal variables
         int nIterations;
         double probMutateMajor;
         double probMutateMinor;
         double probSplit;
         double probPrune;
         double probCrossover;
         double *performanceHistory; // average performance of the trees in the elitism list, is used for the termination of the algorithm
         Tree **trees;
         int elitismRange; // number of trees in the elitism list
         int *treesAge; // number of iterations a tree has not been replace
         int method; // method=1 for classification, and method=6 for regression / other criteria are not sufficently tested yet
         double alpha;  // weights the complexity part of the cost function
         double agePenalty;  // used for the age penalty. value is changed during the run
         double acceptProb;  // used for the age penalty. with a value of 0.25 an average tree is replaced in 4 iterations
         int sumWeights; // weighted variance for regression, sum of weights for classification
         double populationMSE; // population variance; used for regression trees only
         public:
         Container(int* nInst, int* nVar, int *varType, double* ndata, int* weights, int* prediction, int *splitN, int *splitV, double *splitP, int* csplit, int *maxNode, int *minbucket, int* minsplit,
              int* nIter, int* nTrees, int* pMutateMajor, int* pMutateMinor, int* pCrossover, int *pSplit, int* pPrune, int* evCriteria, double* evParameter );
         ~Container();
         void initVariables(int* varType);
         void evolution();
         bool evaluateTree(int treeNumber, bool pruneIfInvalid, int nodeNumber);
         double initMutateNode(int treeNumber, bool minorChange);
         double mutateNode(int treeNumber, int node, bool minorChange);
         int calculateNoOfNodesInSubtree(int treeNumber, int nodeNumber);
         double splitNode(int treeNumber);
         double pruneNode(int treeNumber);
         int pruneAllNodes(int treeNumber);
         bool removeNode(int treeNumber, int nodeNumber);
         double crossover(int treeNumber);
         int getRandomTree(bool elitismTree);
         int getGenitor(void);
         int randomSplitNode(int treeNumber);
         int randomTerminalNode(int treeNumber);
         void randomSplitVariable(int treeNumber, int nodeNumber);
         bool randomSplitPoint(int treeNumber, int nodeNumber);
         bool changeSplitPoint(int treeNumber, int nodeNumber);
         bool changeRandomCategories(int treeNumber, int nodeNumber);
         void overwriteTree(int targetPos);
         void overwriteTree(int sourcePos, int targetPos);
         int evaluateNewSolution(int treeNumber, double* oldPerformance);
         bool updatePerformanceList(int tree);
         int initNVPCrossoverTree1(int treeNumber, int node, int randomNode1, int* tempN, int* tempV, double* tempP, int** csplit);
         int initNVPCrossoverTree2(int treeNumber, int randomNode2, int randomNode1, int* tempN, int* tempV, double* tempP, int** csplit);
};

