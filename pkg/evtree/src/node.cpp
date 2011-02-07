#include "node.h"

Node::Node(int* splitN, int* splitV, double* splitP, int** csplit, Node* leftChild, Node* rightChild, double** data,
int* nInst, int* nVar, variable** variables){
    this->splitN = splitN;
    this->splitV= splitV;
    this->splitP= splitP;
    this->nInst = nInst;
    this->nVar = nVar;
    this->data = data;
    this->leftChild= leftChild;
    this->rightChild= rightChild;
    this->variables= variables;
    this->localClassification = new int[*this->nInst];
    for(int i=0; i< *this->nInst; i++)
        this->localClassification[i]= 0;
    this->sumLocalWeights=0;
    this->sumLeftLocalWeights=0;
    this->sumRightLocalWeights=0;
    this->predictionInternalNode= 0;
    this->predictionLeftTerminal= 0;
    this->predictionRightTerminal= 0;
    this->nClassesDependendVar= &variables[*this->nVar-1]->nCats;
    this->csplit = csplit;
}


int Node::partition( int* classification, int* weights, variable** variables, int* nNodes, int minbucket, int minsplit ){
    // assigns instances to belong to the right or the left child node
    for(int i=0; i< *this->nInst; i++)
            this->localClassification[i]= classification[i];
    this->sumLeftLocalWeights= 0;
    this->sumRightLocalWeights= 0;
    if( this->variables[*this->splitV]->isCat == true){ // categorical split-variable
        bool flag = false;
        for(int i=0; i < *this->nInst; i++){
            if(classification[i] == *this->splitN){
                flag = false;
                for(int k=0; k < variables[ *this->splitV]->nCats && flag==false ; k++){
                    if( variables[ *this->splitV ]->sortedValues[k]==this->data[i][*this->splitV] ){
                       if( this->csplit[k][*this->splitN] == 1 ){
                           classification[i]= (*this->splitN)*2+1;
                           this->localClassification[i]=classification[i];
                           this->sumLeftLocalWeights++;
                       }else{
                           classification[i]= (*this->splitN)*2+2;
                           this->localClassification[i]=classification[i];
                           this->sumRightLocalWeights++;
                       }
                       flag=true;
                    }
               }
           }
       }
    }else if( variables[*this->splitV]->isCat == false){  // numeric split-variable
        for(int i=0; i < *this->nInst; i++){
            if(classification[i] == *this->splitN){
                    if( (double)this->data[i][ *this->splitV ] < (double)*this->splitP  ){
                            classification[i]= (*this->splitN)*2+1;
                            this->sumLeftLocalWeights += weights[i];
                            this->localClassification[i]=classification[i];
                    }else{
                            classification[i]= (*this->splitN)*2+2;
                            this->sumRightLocalWeights += weights[i];
                            this->localClassification[i]=classification[i];
                    }
            }
        }

    }
    this->sumLocalWeights = this->sumLeftLocalWeights + this->sumRightLocalWeights ;
    if( this->sumLocalWeights < minsplit  && *this->splitN > 0){  // checks if there are enough instances in the nodes; otherwise the node is pruned
        return (int) *this->splitN;
    }

    // recursive call of partition until there are no further internal nodes
    int temp1=-1, temp2=-1;
    if( this->leftChild != NULL){
       temp1 = this->leftChild->partition( classification, weights, variables, nNodes, minbucket, minsplit );
    }
    if( this->rightChild != NULL){
       temp2= this->rightChild->partition( classification, weights, variables, nNodes, minbucket, minsplit );
    }

    if( temp1 == -2 || temp2 == -2){
           return -2;
    }else if( temp1 == 0 || temp2 == 0){
           return 0;
    }else if(temp1 != -1){
        return temp1;
    }else if(temp2 != -1){
        return temp2;
    }

    if(this->sumLeftLocalWeights < minbucket){   // if there are to few instances in the left terminal-node the internal node is pruned
        return *this->splitN;
    }else if(this->sumRightLocalWeights < minbucket){ // if there are to few instances in the right terminal-node the internal node is pruned
        return *this->splitN;
    }
    return -1;
} // end partition


double Node::calculateNodeMC(int* weights){
    // calculate the fraction of correctly classified instances
    double sumWeights= 0;
    double *sumsClassification = new double[*this->nClassesDependendVar];
    for(int i=0; i<*this->nClassesDependendVar; i++){
        sumsClassification [i]= 0.0;
    }
    for(int i=0; i<*this->nInst; i++){
        if( localClassification[i] == (*this->splitN)*2+1 || localClassification[i] == (*this->splitN)*2+2) {
            sumsClassification[ (int)data[i][ *this->nVar-1] -1] += weights[i];
            sumWeights += weights[i];
        }
    }
    double correctClassified= sumsClassification [0];
    this->predictionInternalNode= 0;
    for(int j=1; j<*this->nClassesDependendVar; j++){
        if(  sumsClassification [j] > correctClassified  ){
            correctClassified= sumsClassification [j];
            this->predictionInternalNode= j;
        }
    }

    delete [] sumsClassification;
    sumsClassification= NULL;
    return ((double)correctClassified) / ((double)sumWeights);
}


double Node::calculateNodeSE(int* weights){
        double nodeMean=0;
        double squaredSum=0;
        int sumWeights= 0;

        for(int i=0; i<*this->nInst; i++){
            if( this->localClassification[i] == (*this->splitN)*2+1 || this->localClassification[i] == (*this->splitN)*2+2 ){
                 nodeMean += data[i][*this->nVar-1]*weights[i];
                 squaredSum += data[i][*this->nVar-1]*data[i][*this->nVar-1]*weights[i];
                 sumWeights += weights[i];
            }
        }

        nodeMean /= ((double)sumWeights);
        this->predictionInternalNode= nodeMean;
        return 1.0/((double)sumWeights)*squaredSum-(nodeMean*nodeMean);
}


double Node::calculateChildNodePerf(bool leftNode, int method, int* weights){
        // calculate the performance depending of the cost criterium
        double performance_cc= 0;
        double performance= 0;
        int sumWeights=0;
        double *sumsClassification = new double[*this->nClassesDependendVar];

        for(int i=0; i<*this->nClassesDependendVar; i++){
                sumsClassification [i]= 0.0;
        }
        if(leftNode==true){
            for(int i=0; i<*this->nInst; i++){
                if( localClassification[i] == (*this->splitN)*2+1){
                     sumsClassification[ (int)data[i][*this->nVar-1]-1]  += weights[i];
                     sumWeights += weights[i];
                }
            }
        }else{
            for(int i=0; i<*this->nInst; i++){
                if(localClassification[i] == (*this->splitN)*2+2 ){
                     sumsClassification[ (int)data[i][*this->nVar-1]-1]  += weights[i];
                     sumWeights += weights[i];
                }
            }
        }
        int localMajorityClassVariable=0;
        performance_cc= sumsClassification[0];
        if(method == 4){  // MDL criteria
            if( sumsClassification[0] != 0 ){
                performance = sumsClassification[0]*log2(((double)sumWeights)/sumsClassification[0]);
            }

            for(int i=1; i<*this->nClassesDependendVar; i++){
                if(  sumsClassification [i] > performance_cc  ){  // CC part
                    performance_cc= sumsClassification [i];
                    localMajorityClassVariable=i;
                }
                if( sumsClassification[i] != 0 ){   // BIC part
                    performance += sumsClassification[i]*log2(((double)sumWeights)/sumsClassification[i]) ;
                }
            }
            performance  += (double)(*this->nClassesDependendVar-1)/2.0*log2(((double)sumWeights)/2.0)
                         + log2(pow(pi, (double)(*this->nClassesDependendVar)/2.0) / (double)factorial( int(ceil((double)(*this->nClassesDependendVar)/2.0)) )) ;


            if(variables[*this->splitV]->isCat){
                performance+= log2((double)pow(2,variables[*this->splitV]->nCats)-2);///MDLTEST//
            }else{
                performance+= log2((double)(variables[*this->splitV]->nCats-1));
            }
        }else{
            if( method > 1 && sumsClassification[0] != 0 ){   // BIC crieteria only
                performance = sumsClassification[0]*log(((double)sumWeights)/sumsClassification[0]);
            }

            for(int i=1; i<*this->nClassesDependendVar; i++){
                if(  sumsClassification [i] > performance_cc  ){  // CC part
                    performance_cc= sumsClassification [i];
                    localMajorityClassVariable=i;
                }
                if(method > 1 &&  sumsClassification[i] != 0 ){   // BIC part
                    performance += sumsClassification[i]*log(((double)sumWeights)/sumsClassification[i]);
                }
            }
        }
        delete [] sumsClassification;
        sumsClassification= NULL;

        if(leftNode==true){
            this->predictionLeftTerminal= localMajorityClassVariable;
            this->performanceLeftTerminal= performance_cc/((double)sumWeights);
        }else{

            this->predictionRightTerminal= localMajorityClassVariable;
            this->performanceRightTerminal= performance_cc/((double)sumWeights);
        }

        if(method<=1)
            return performance_cc;
        else
            return performance;
}


double Node::calculateChildNodeSE(bool leftNode, int* weights){
        // calculate performance for Criteria MSE
        double performance=0;
        double nodeMean=0;
        double squaredSum=0;
        int sumWeights= 0;
        if(leftNode==true){
            for(int i=0; i<*this->nInst; i++){
                if( localClassification[i] == (*this->splitN)*2+1){
                     nodeMean += data[i][*this->nVar-1]*weights[i];
                     squaredSum += data[i][*this->nVar-1]*data[i][*this->nVar-1]*weights[i];
                     sumWeights += weights[i];
                }
            }
        }else{
            for(int i=0; i<*this->nInst; i++){
                if(localClassification[i] == (*this->splitN)*2+2 ){
                     nodeMean += data[i][*this->nVar-1]*weights[i];
                     squaredSum += data[i][*this->nVar-1]*data[i][*this->nVar-1]*weights[i];
                     sumWeights += weights[i];
                }
            }
        }

        nodeMean /= ((double)sumWeights);

        performance= ((double)sumWeights)*(1.0/(double)(sumWeights)*squaredSum-(nodeMean*nodeMean));

        if(leftNode==true){
            this->performanceLeftTerminal= performance/((double)sumWeights);
            this->predictionLeftTerminal= nodeMean;
        }else{
            this->performanceRightTerminal= performance/((double)sumWeights);
            this->predictionRightTerminal= nodeMean;
        }
        return performance;
}


int Node::factorial( int n ){
    // recursively calculates the factorial of a number
    if ( n <= 1 )
        return 1;
    else
        return  n * factorial( n-1 );
}


Node::~Node(){
       delete [] localClassification;
       localClassification= NULL;
       leftChild= NULL;
       rightChild= NULL;
       splitN= NULL;
       splitV= NULL;
       splitP= NULL;
       csplit= NULL;
       leftChild= NULL;
       rightChild= NULL;
       nInst= NULL;
       nVar= NULL;
       data= NULL;
}

