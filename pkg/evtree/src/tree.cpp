#include "tree.h"

Tree::Tree(int* nInstances, int* nVariables, double** data, int* weights, int *splitN, int *splitV, double *splitP, int** csplit, int* maxCat, int* nNodes, variable** variables, int* maxNode){
        // constructor used by crossover
        this->nInstances= nInstances;
        this->nVariables= nVariables;
        this->nNodes= *nNodes;
        this->maxNode= maxNode;
        this->maxCat= maxCat;
        this->splitV= new int[*this->maxNode];
        this->splitP= new double[*this->maxNode];
        this->splitN= new int[*this->maxNode];
        this->variables= variables;
        this->nodes= new Node*[*this->maxNode];
        this->classification = new int[*this->nInstances];
        this->data= data;
        this->performance= 999999;
        this->csplit= new int*[*this->maxCat];
        this->weights= weights;
	for (int i = 0; i < *this->maxCat; i++)
            this->csplit[i] = new int[(*this->maxNode)];
        for (int v = 0; v < *this->maxNode ; v++){
            for(int i=0; i< *this->maxCat; i++){
                this->csplit[i][v]= csplit[i][v];
            }
            this->splitV[v] = splitV[v];
            this->splitP[v] = splitP[v];
            this->splitN[v] = splitN[v];
        }
        for (int nodeNumber = *this->maxNode-1; nodeNumber >= 0; nodeNumber--){
               this->nodes[nodeNumber]= NULL;
               this->initNode(nodeNumber);
        }

}


Tree::Tree(int* nInstances, int* nVariables, double** data, int* weights, int* maxCat, variable** variables, int* maxNode, int* minbucket, int* minsplit){
         // initializes a tree with a random root node
        this->nInstances= nInstances;
        this->nVariables= nVariables;
        this->nNodes= 1;
        this->maxNode= maxNode;
        this->maxCat= maxCat;
        this->splitV= new int[*this->maxNode];
        this->splitP= new double[*this->maxNode];
        this->splitN= new int[*this->maxNode];
        this->variables= variables;
        this->nodes= new Node*[*this->maxNode];
        this->classification = new int[*this->nInstances];
        this->data= data;
        this->performance= 999999;
        this->csplit= new int*[*this->maxCat];
        this->weights= weights;
        
	for (int i = 0; i < *this->maxCat; i++)
            this->csplit[i] = new int[(*this->maxNode)];
        for (int v = 0; v < *this->maxNode ; v++){
            for(int i=0; i< *this->maxCat; i++){
                this->csplit[i][v]= 2;
            }
            this->splitN[v]= -999999;
            this->splitV[v]= -999999;
            this->splitP[v]= -999999;
            this->nodes[v]= NULL;
        }
        this->splitN[0] = 0;
        this->splitV[0]= (rand()%(*this->nVariables-1));
        this->nodes[0]= NULL;
        this->initNode(0);

        if(variables[this->splitV[0]]->isCat==false){
              if((this->variables[this->splitV[0]]->nCats-1) > 1 )
                   this->splitP[0]= variables[this->splitV[0]]->sortedValues[(rand()%(this->variables[this->splitV[0]]->nCats-1))+1];
              else
                   this->splitP[0]= variables[this->splitV[0]]->sortedValues[0];
        }else{
              this->randomizeCategories(0);
        }
        int dint;
        int check=0;
        while(this->predictClass(*minbucket, *minsplit, false, 0)==false){
            dint= (rand()%(*this->nVariables-1));
            this->splitV[0]= dint;
            if(variables[this->splitV[0]]->isCat==false){
                  if((this->variables[this->splitV[0]]->nCats-1) > 1 )
                       this->splitP[0]= variables[this->splitV[0]]->sortedValues[(rand()%(this->variables[this->splitV[0]]->nCats-1))+1];
                  else
                       this->splitP[0]= variables[this->splitV[0]]->sortedValues[0];
            }else{
                 this->randomizeCategories(0);
             }
            check++;
            if(check==1000){
                cout << "tree could not be initialized!" << endl;
                exit(0);
            }
        }
}


Tree::~Tree(){
	for (int i = 0; i < *this->maxNode; i++) {
           delete nodes[i];
	}
	delete [] nodes;
        nodes= NULL;
	delete [] classification;
        classification=NULL;
	delete [] splitP;
        splitP=NULL;
	delete [] splitV;
        splitV=NULL;
	delete [] splitN;
        splitN= NULL;
        for (int i = 0; i < *this->maxCat; i++)
            delete [] csplit[i];
        delete [] csplit;
        csplit= NULL;
        variables= NULL;
        data=NULL;
	maxNode= NULL;
        maxCat= NULL;
        nInstances= NULL;
        nVariables= NULL;
        weights= NULL;
}


void Tree::initNode(int nodeNumber){
    // initializes a node
    if(this->splitN[nodeNumber] >= 0 ){
        int leftChild= -1;
        int rightChild=-1;
        // is leaf node?
        if( nodeNumber*2+2 < *this->maxNode){
            if( (this->splitN[nodeNumber*2+1]  ) ==  nodeNumber*2+1 ){
                    leftChild=nodeNumber*2+1;
            }
            if( (this->splitN[nodeNumber*2+2]) ==  nodeNumber*2+2){
                    rightChild=nodeNumber*2+2;
            }
        }

        if( leftChild <= 0 && rightChild <= 0){
                this->nodes[nodeNumber]= new Node(&this->splitN[nodeNumber], & this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit, NULL, NULL,
this->data, this->nInstances, this->nVariables,  this->variables);
        }else if( leftChild <= 0 ){
                this->nodes[nodeNumber]= new Node(&this->splitN[nodeNumber], &this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit, NULL, this->nodes[rightChild] ,
this->data, this->nInstances, this->nVariables,   this->variables);
        }else if( rightChild <= 0){
                this->nodes[nodeNumber]= new Node(&this->splitN[nodeNumber], &this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit,this->nodes[leftChild], NULL,
this->data, this->nInstances, this->nVariables,  this->variables);
        }else{
                this->nodes[nodeNumber]= new Node(&this->splitN[nodeNumber], & this->splitV[nodeNumber], &this->splitP[nodeNumber],this->csplit, this->nodes[leftChild], this->nodes[rightChild],
this->data, this->nInstances, this->nVariables, this->variables );
        }
    }else{
        this->nodes[nodeNumber] = NULL;
    }

}


int Tree::predictClass(int minbucket, int minsplit, bool pruneIfInvalid, int nodeNumber){
    // predict the class membership
    // if pruneIfInvalid == TRUE non-valid nodes a pruned
    // otherwise -1 is returned for non-valid nodes
    if(nodeNumber == 0){
        for(int i =0; i<*this->nInstances; i++){
            this->classification[i]=0;
        }
    }else{
        this->reverseClassification(nodeNumber, nodeNumber);
    }

    int returnValue=  this->nodes[nodeNumber]->partition( this->classification, this->weights, this->variables, &this->nNodes, minbucket, minsplit );
    if(returnValue == -1){
        return -1;
    }else if(returnValue < -1 || returnValue == 0 || pruneIfInvalid == false){
        return returnValue;
    }else{
       this->deleteChildNodes(returnValue); // call recursion delete node and everything below it
       return predictClass(minbucket, minsplit, true, 0 );
    }
}


bool Tree::reverseClassification(int startNode, int nodeNumber){
    // all observations which are in the subtree below "startNode" are classified to be in node "startNode"
    // this saves computation time when the split-rule in startNode is changed. after the call of "reverseClassification" only the instances in
    // node "startnode" are newly evaluated. The rest of the tree stays the same
    for(int i=0; i<*this->nInstances; i++){
        if(this->classification[i] == (nodeNumber*2+1) || this->classification[i] == (nodeNumber*2+2) ){
            this->classification[i]= startNode;
        }
    }
    if(nodeNumber*2+1 < *this->maxNode){
        if(splitN[nodeNumber] == nodeNumber){
             reverseClassification(startNode, nodeNumber*2+1);
        }
    }
    if(nodeNumber*2+2 < *this->maxNode){
        if(splitN[nodeNumber] == nodeNumber){
            reverseClassification(startNode, nodeNumber*2+2);
        }
    }
    return true;
}


bool Tree::deleteChildNodes(int nodeNumber){
    // used by predictClass() and mutadeMode() to delete a child node
    if(this->splitN[nodeNumber] > 0  && nodeNumber > 0){
        if(this->nodes[nodeNumber]->leftChild != NULL ){
            this->deleteChildNodes(nodeNumber*2+1);
        }
        if(this->nodes[nodeNumber]->rightChild != NULL ){
            this->deleteChildNodes(nodeNumber*2+2);
        }
        if(nodeNumber%2 == 0){
            this->nodes[(int) ((nodeNumber-1) /2)  ]->rightChild= NULL;
        }else{
            this->nodes[(int) ((nodeNumber-1) /2)  ]->leftChild= NULL;
        }
        this->splitN[nodeNumber]= -999999;
        this->splitV[nodeNumber]= -999999;
        this->splitP[nodeNumber]= -999999;
        this->nNodes--;
        delete this->nodes[nodeNumber];
        this->nodes[nodeNumber]= NULL;
        return true;
    }else{
        cout << "warning: mode could not be deleted " << endl;
        return false;
    }
}


void Tree::randomizeCategories(int nodeNumber){
    // assigns random categories of a categorical variable
    bool left=false;
    bool right=false;
    for(int i=0; i< this->variables[ *this->nodes[nodeNumber]->splitV ]->nCats ; i++){
        if( i==this->variables[ *this->nodes[nodeNumber]->splitV ]->nCats-1 && left == false){
            this->csplit[i][nodeNumber]=1;
        }else if( i==this->variables[ *this->nodes[nodeNumber]->splitV ]->nCats-1 && right == false){
            this->csplit[i][nodeNumber]=3;
        }else if( rand()%2==1 ){
            this->csplit[i][nodeNumber]=1;
            left=true;
        }else{
            this->csplit[i][nodeNumber]=3;
            right=true;
        }
    }

}


bool Tree::calculateTotalCosts(int method, double alpha, int sumWeights, double populationMSE){
    // only 1 and 6 are implemented with weights and sufficently tested
 //   if(method == 0){
 //   this->performance= 2*(double(*this->nInstances)-this->calculateTotalMC(0)*(double(*this->nInstances)/((double) sumWeights)) +  alpha*(this->nNodes+1.0)*log((double)(*this->nInstances));
     //  this->performance= 2*(double(*this->nInstances)-this->calculateTotalMC(0)*(double(*this->nInstances)/dataSSE)) +  alpha*(*this->nodes[0]->nClassesDependendVar-1)*(this->nNodes+1)*log((double)(*this->nInstances));
    if(method == 1){
       this->performance= 2.0*(((double) sumWeights)-this->calculateTotalMC(0)) + alpha*(this->nNodes+1.0)*log(((double)sumWeights));
 /*   }else if(method == 2){
       this->performance= 2*this->calculateTotalBIC(0)+  ((*this->nodes[0]->nClassesDependendVar+1)*(this->nNodes+1)-1)*log((double)(*this->nInstances));
    }else if(method == 3) {
       this->performance= 2*this->calculateTotalBIC(0)+  ((*this->nodes[0]->nClassesDependendVar-1)*(this->nNodes+1))*log((double)(*this->nInstances));
    }else if(method == 4) {
       this->performance= this->calculateTotalMDL(0)+
            this->nNodes*2+1+
            this->nNodes*log2((double)(*this->nVariables));
    }else if(method == 5){
       this->performance= 1.0-this->calculateTotalMC(0)+ alpha*(this->nNodes+1);*/
    }else if(method == 6){
               double SMSE= max(this->calculateTotalSE(0)/(populationMSE), 0.001);
               this->performance= (
                  ((double) sumWeights)*log(SMSE)+alpha*4.0*log(((double) sumWeights))*((double)this->nNodes+2.0)
               +  ((double) sumWeights)*7.0  // constant such that formula is alway positive
               );
    }
    return true;
}


double Tree::calculateTotalSE(int nodeNumber){
    double performance=0;
    if(this->nodes[nodeNumber]->leftChild != NULL)
        performance += this->calculateTotalSE(nodeNumber*2+1);
    if(this->nodes[nodeNumber]->rightChild != NULL)
        performance += this->calculateTotalSE(nodeNumber*2+2);
    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->leftChild == NULL){
        performance += this->nodes[nodeNumber]->calculateChildNodeSE(true, this->weights);
    }
    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->rightChild == NULL){
        performance +=  this->nodes[nodeNumber]->calculateChildNodeSE(false, this->weights);
    }
    return performance;
}


double Tree::calculateTotalMC(int nodeNumber){
    double performance=0;
    if(this->nodes[nodeNumber]->leftChild != NULL)
        performance += this->calculateTotalMC(nodeNumber*2+1);
    if(this->nodes[nodeNumber]->rightChild != NULL)
        performance += this->calculateTotalMC(nodeNumber*2+2);

    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->leftChild == NULL){
        performance += this->nodes[nodeNumber]->calculateChildNodePerf(true, 1, this->weights);
    }
    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->rightChild == NULL){
        performance +=  this->nodes[nodeNumber]->calculateChildNodePerf(false, 1, this->weights);
    }

    return performance;
}


double Tree::calculateTotalBIC(int nodeNumber){
    double performance=0;
    if(this->nodes[nodeNumber]->leftChild != NULL)
        performance += this->calculateTotalBIC(nodeNumber*2+1);
    if(this->nodes[nodeNumber]->rightChild != NULL)
        performance += this->calculateTotalBIC(nodeNumber*2+2);

    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->leftChild == NULL){
        performance += this->nodes[nodeNumber]->calculateChildNodePerf(true, 2, this->weights);
    }
    if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->rightChild == NULL){
        performance +=  this->nodes[nodeNumber]->calculateChildNodePerf(false, 2, this->weights);
    }

    return performance;
}


double Tree::calculateTotalMDL(int nodeNumber){
      double performance=0;
      if(this->nodes[nodeNumber]->leftChild != NULL)
          performance += this->calculateTotalMDL(nodeNumber*2+1);
      if(this->nodes[nodeNumber]->rightChild != NULL)
          performance += this->calculateTotalMDL(nodeNumber*2+2);

      if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->leftChild == NULL){
          performance += this->nodes[nodeNumber]->calculateChildNodePerf(true, 4, this->weights);
      }
      if( this->splitN[nodeNumber] == nodeNumber && this->nodes[nodeNumber]->rightChild == NULL){
          performance +=  this->nodes[nodeNumber]->calculateChildNodePerf(false, 4, this->weights);
      }
      return performance;
}



void Tree::printTree(int method){ /// for debugging
    cout << "------------------------------------------------------------" << endl;
    if(method == 6){
        cout << "MSE: " << this->calculateTotalSE(0)/(*this->nInstances) << ", Number of Nodes " << this->nNodes << ", Quality:" << this->performance << endl;
        printNode(0, method); // init recursive call of node 0
    }else{
        cout << "Correct Classified: " << this->calculateTotalMC(0)/(*this->nInstances) << ", Number of Nodes " << this->nNodes << ", Quality:" << this->performance*100 << endl;
        printNode(0, method); // init recursive call of node 0
    }
}


void Tree::printNode(int i, int method){ /// for debugging
    int exp=0;
    while( pow(2,exp) <= i+1 ){ ///
        cout << " ";
        exp++;
    }
    cout << i+1 << ") "; ///
    cout << splitV[i];
    if( variables[abs(this ->splitV[i])]->isCat == false){
        if( this->splitV[i] >= 0 ){
               cout << " < " << (double)this->splitP[i];
        }else{
             cout << " >= " << (double)this->splitP[i];
        }
    }else{
        cout << " = ";
        bool printComma= false;
        for(int k=0; k<variables[this->splitV[i]]->nCats; k++){
            if(this->csplit[k][i]==1){
                if (printComma == true)
                       cout << ", " ;
              cout  << variables[ this->splitV[i]]->sortedValues[k];
                printComma= true;
            }
        }
    }
    if(method < 6)
        cout << ", NI: " << this->nodes[i]->sumLocalWeights << ", Class: " << this->nodes[i]->predictionInternalNode << ", CC: " << this->nodes[i]->calculateNodeMC(this->weights) << "%" << endl;
    else
        cout << ", NI: " << this->nodes[i]->sumLocalWeights << ", Mean: " << this->nodes[i]->predictionInternalNode << ", MSE: " << this->nodes[i]->calculateNodeSE(this->weights) << "" << endl;


    if(this->splitN[(i*2)+1]==i*2+1 && i*2+1 < *this->maxNode){
            printNode(i*2+1, method);
    }
    if(this->splitN[(i*2)+2]==i*2+2 && i*2+2 < *this->maxNode){
            printNode(i*2+2, method);
    }
}


