#include "node.h"
using namespace std;

class Tree{
	 public:
         int *nInstances;
         int *nVariables;
	 variable **variables;
         double **data;
	 int *maxNode;
         int * maxCat;
	 int *splitV;
	 double *splitP;
         int *weights;
   	 int *splitN;
         int** csplit;
	 int nNodes;
	 int *classification;	
	 int *posDependendVar;
	 Node **nodes;
	 double performance;  

	 public:
	 Tree(int* nInst, int* nVariables, double** data, int* weights, int *splitN, int *splitV, double *splitP, int** csplit, int* maxCat, int* nNodes, variable** variables, int* maxNode);
	 Tree(int* nInst, int* nVariables, double** data, int* weights, int* maxCat, variable** variables, int* maxNode, int* minbucket, int* minsplit);
         ~Tree();
 	 int predictClass(int minbucket, int minsplit, bool crossover, int nodeNumber);
	 void initNode(int nodeNo);
 	 double calculateTotalMC(int nodeNumber);
         bool calculateTotalCosts(int method, double alpha, int sumWeights, double populationMSE);
	 double calculateTotalBIC(int nodeNumber);
         double calculateTotalMDL(int nodeNumber);
         double calculateTotalSE(int nodeNumber);
  	 void printTree(int evCriteria);
	 void printNode(int nodeNo, int evCriteria);
  	 void printTerminalsOf(int nodeNo);
         bool deleteChildNodes(int nodeNo);
         bool reverseClassification(int startNode, int nodeNumber);
         void randomizeCategories(int nodeNumber);
};
