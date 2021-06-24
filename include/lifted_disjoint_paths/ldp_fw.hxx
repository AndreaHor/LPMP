#ifndef LDP_FW_HXX
#define LDP_FW_HXX


#include "FW-MAP.h"
#include <cassert>
#include <array>
#include <iostream>
#include <cmath>
#include "ldp_instance.hxx"

namespace LPMP{

//struct InOutFlowSubproblem {

//  size_t nodeID;
//  bool isOutflow;
//  //std::vector<size_t> traverseOrder;
//};


struct OptimalSolution{
    size_t baseEdgeIndex;
    std::vector<size_t> liftedEdgesIndices;
    double optValue;
};


class LdpProblem{
public:
    LdpProblem(const lifted_disjoint_paths::LdpInstance* pLdpInstance){
        pInstance=pLdpInstance;
    }

    size_t getIndexInYLifted(size_t centralNodeID)const{
        assert(centralNodeID<numberOfNodes);
        size_t index=centralNodeID*(maxTimeGap+1)+1;
        return index;
    }

    size_t getIndexInYBase(size_t centralNodeID) const{
        size_t index=centralNodeID*(maxTimeGap+1);
        return index;
    }

    //Can be called just to get the index of the first neighbor and then iterate without calling this function always
    size_t getIndexInWILifted(size_t centralNodeID,size_t neighborIndex)const{
        size_t index=0;
        if(isOutFlow){
            index=pInstance->getMyGraphLifted().getIndexForward(centralNodeID,neighborIndex);
        }
        else{
            index=pInstance->getMyGraphLifted().getIndexBackward(centralNodeID,neighborIndex);
        }
        return index+numberOfNodes+numberOfBaseEdges;
    }

    //Can be called just to get the index of the first neighbor and then iterate without calling this function always
    size_t getIndexInWIBase(size_t centralNodeID,size_t neighborIndex)const{
        size_t index=0;
        if(isOutFlow){
            index=pInstance->getMyGraph().getIndexForward(centralNodeID,neighborIndex);
        }
        else{
            index=pInstance->getMyGraph().getIndexBackward(centralNodeID,neighborIndex);
        }
        return index+numberOfNodes;
    }

    size_t getIndexInWINode(size_t centralNodeID){
        return centralNodeID;
    }


    size_t getYLength()const{
        return yLength;
    }

    LdpProblem(const lifted_disjoint_paths::LdpInstance* _pInstance,bool _isOutFlow)
    {
        pInstance=_pInstance;
        numberOfNodes=pInstance->getNumberOfVertices()-2;
        numberOfBaseEdges=pInstance->getMyGraph().getNumberOfEdges();
        numberOfLiftedEdges=pInstance->getMyGraph().getNumberOfEdges();
        isOutFlow=_isOutFlow;
        maxTimeGap=pInstance->parameters.getMaxTimeLifted();
        yLength=numberOfNodes*(maxTimeGap+1);
        xLength=numberOfNodes+pInstance->getMyGraph().getNumberOfEdges()+pInstance->getMyGraphLifted().getNumberOfEdges();

        nodeIDToIndex=std::vector<std::map<size_t,size_t>>(numberOfNodes);
        if(isOutFlow){
            for (size_t i = 0; i < numberOfNodes; ++i) {
                size_t indexCounter=0;
                const LdpDirectedGraph::edge* edgeIt=pInstance->getMyGraphLifted().forwardNeighborsBegin(i);
                for (;edgeIt!=pInstance->getMyGraphLifted().forwardNeighborsEnd(i);edgeIt++) {
                    size_t nodeID=edgeIt->first;
                    nodeIDToIndex[i].insert(nodeID,indexCounter);
                    indexCounter++;
                }
            }
        }
        else{
            for (size_t i = 0; i < numberOfNodes; ++i) {
                size_t indexCounter=0;
                const LdpDirectedGraph::edge* edgeIt=pInstance->getMyGraphLifted().backwardNeighborsBegin(i);
                for (;edgeIt!=pInstance->getMyGraphLifted().backwardNeighborsEnd(i);edgeIt++) {
                    size_t nodeID=edgeIt->first;
                    nodeIDToIndex[i].insert(nodeID,indexCounter);
                    indexCounter++;
                }
            }

        }

        traverseOrders=std::vector<std::vector<size_t>>(numberOfNodes);
        if(isOutFlow){
            for (size_t i = 0; i < numberOfNodes; ++i) {
                const std::set<size_t>& reachableVertices=pInstance->reachableFromVertex(i);
                for (auto it=reachableVertices.begin();it!=reachableVertices.end();it++) {
                    traverseOrders[i].push_back(*it);
                }
            }
        }
        else{
            assert(numberOfNodes>0);
            for (size_t i = numberOfNodes-1; i >=0; --i) {
                const std::set<size_t>& reachableVertices=pInstance->reachableFromVertex(i);
                for (auto it=reachableVertices.rbegin();it!=reachableVertices.rend();it++) {
                    assert(*it<numberOfNodes);
                    traverseOrders[*it].push_back(i);
                }
            }
        }
    }

    size_t getNumberOfNodes() const{
        return numberOfNodes;
    }

    bool getIsOutFlow()const{
        return isOutFlow;
    }

    size_t nodeIDtoNeighborIndex(const size_t& centralNodeID, const size_t& neighborID){
        assert(centralNodeID<numberOfNodes);
        auto f=nodeIDToIndex[centralNodeID].find(neighborID);
        if(f!=nodeIDToIndex[centralNodeID].end()){
            return f->second;
        }
        else{
            throw std::invalid_argument("Required index of non-existing neighbor of node "+std::to_string(centralNodeID));
        }
    }

    double topDownMethod(size_t centralNodeID,double* wi,size_t* y){
    //TODO: top down method here

        const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
        const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();

       assert(centralNodeID<numberOfNodes);
            for(size_t v:traverseOrders[centralNodeID]){
                pInstance->sncTDStructure[v]=0;
                pInstance->sncNeighborStructure[v]=getVertexToReach();
                pInstance->isBSF[v]=0;

            }
            size_t numberOfLiftedNeighbors=0;
            if(isOutFlow){
                numberOfLiftedNeighbors=liftedGraph.getNumberOfEdgesFromVertex(centralNodeID);
                size_t counter=0;
                for (auto it=liftedGraph.forwardNeighborsBegin(centralNodeID);it!=liftedGraph.forwardNeighborsEnd(centralNodeID);it++) {
                    pInstance->sncTDStructure[it->first]=it->second+wi[getIndexInWILifted(centralNodeID,counter)];
                    counter++;
                }
            }
            else{
                numberOfLiftedNeighbors=liftedGraph.getNumberOfEdgesToVertex(centralNodeID);
                size_t counter=0;
                for (auto it=liftedGraph.backwardNeighborsBegin(centralNodeID);it!=liftedGraph.backwardNeighborsEnd(centralNodeID);it++) {
                    pInstance->sncTDStructure[it->first]=it->second+wi[getIndexInWILifted(centralNodeID,counter)];
                    counter++;
                }
            }



        double bsfValue=std::numeric_limits<double>::max();
        size_t i=0;


        for (; i < traverseOrders[centralNodeID].size(); ++i) {

            size_t currentNode=traverseOrders[centralNodeID][i];

            double bestDescValue=0;
            size_t bestDescVertexID=getVertexToReach();

            //Search for best descendant
            if(isOutFlow){
                const LdpDirectedGraph::edge* vertexIt=baseGraph.forwardNeighborsBegin(currentNode);
                const LdpDirectedGraph::edge* end=baseGraph.forwardNeighborsEnd(currentNode);


                for (;vertexIt!=end;vertexIt++) {

                    size_t desc=vertexIt->first;

                    if(desc>=numberOfNodes) continue;

                    if(desc<=*traverseOrders[centralNodeID].rbegin()){

                        double value=pInstance->sncTDStructure[desc];
                        if(bestDescValue>value){
                            bestDescValue=value;
                            bestDescVertexID=desc;

                        }
                        if(pInstance->isBSF[desc]){
                            break;
                        }
                    }
                    else{
                        break;
                    }
                }
            }
            else{
                const LdpDirectedGraph::edge* vertexIt=baseGraph.backwardNeighborsEnd(currentNode);

                const LdpDirectedGraph::edge* begin=baseGraph.backwardNeighborsBegin(currentNode);
                bool doSearch=vertexIt!=begin;
                if(doSearch) vertexIt--;
                while(doSearch){
                    size_t desc=vertexIt->first;

                    if(desc>=numberOfNodes){
                        if(vertexIt==begin){
                            doSearch=false;
                        }
                        else{
                            vertexIt--;
                        }
                    }
                    else{

                        if(desc>=*traverseOrders[currentNode].begin()){

                            double value=pInstance->sncTDStructure[desc];
                            if(bestDescValue>value){
                                bestDescValue=value;
                                bestDescVertexID=desc;

                            }
                            if(pInstance->isBSF[desc]){
                                doSearch=false;
                            }

                            if(vertexIt==begin){
                                doSearch=false;
                            }
                            else{
                                vertexIt--;
                            }

                        }
                        else{
                            doSearch=false;
                        }
                    }


                }

            }


            double value=pInstance->sncTDStructure[currentNode]+bestDescValue;
            if(value<bsfValue){
              //  bsfVector.push_back(currentNode);
                bsfValue=value;
                pInstance->isBSF[currentNode]=true;
            }

            pInstance->sncTDStructure[currentNode]+=bestDescValue;
            pInstance->sncNeighborStructure[currentNode]=bestDescVertexID;

        }


        //all nodes closed, compute solution values
       const LdpDirectedGraph::edge* baseNeighborsIt=nullptr;
       const LdpDirectedGraph::edge* baseNeighborsEnd=nullptr;

       size_t nodeNotActive=0;
       if(isOutFlow){
          baseNeighborsIt=baseGraph.forwardNeighborsBegin(centralNodeID);
          baseNeighborsEnd=baseGraph.forwardNeighborsEnd(centralNodeID);
          nodeNotActive=baseGraph.getNumberOfEdgesFromVertex(centralNodeID);
       }
       else{
           baseNeighborsIt=baseGraph.backwardNeighborsBegin(centralNodeID);
           baseNeighborsEnd=baseGraph.backwardNeighborsEnd(centralNodeID);
           nodeNotActive=baseGraph.getNumberOfEdgesToVertex(centralNodeID);
       }


        double bestSolutionValue=0;
        size_t bestSolutionIndex=nodeNotActive;
        size_t bestSolutionID=std::numeric_limits<size_t>::max();

        std::vector<double> solutionCosts(nodeNotActive+1);
        solutionCosts[nodeNotActive]=0;
        size_t counter=0;
        for (; baseNeighborsIt!=baseNeighborsEnd; baseNeighborsIt++) {

            double baseCost=baseNeighborsIt->second;
            size_t index=getIndexInWIBase(centralNodeID,counter);
            baseCost+=wi[index];

            size_t neighborID=baseNeighborsIt->first;


                double valueToAdd=0;
                if(neighborID<numberOfNodes){
                    valueToAdd=pInstance->sncTDStructure[neighborID];
                }
                double nodeCost=pInstance->getVertexScore(neighborID);
                double value=baseCost+nodeCost+valueToAdd;

                solutionCosts[counter]=value;
                if(value<bestSolutionValue){
                    bestSolutionValue=value;
                    bestSolutionIndex=i;
                    bestSolutionID=neighborID;
                }

            counter++;
        }

        std::vector<size_t> optLiftedIndices;
        if(bestSolutionIndex!=nodeNotActive){
//            assert(bestSolutionID<numberOfNodes);
            while(bestSolutionID<numberOfNodes){
                if((isOutFlow&&pInstance->existLiftedEdge(centralNodeID,bestSolutionID))||(!isOutFlow&&pInstance->existLiftedEdge(bestSolutionID,centralNodeID))){
                    size_t liftedIndex=nodeIDToIndex[centralNodeID][bestSolutionID];
                    optLiftedIndices.push_back(liftedIndex);
                }
                bestSolutionID=pInstance->sncNeighborStructure[bestSolutionID];
            }
        }
        OptimalSolution optSolution;
        optSolution.optValue=bestSolutionValue;
        optSolution.baseEdgeIndex=bestSolutionIndex;
        optSolution.liftedEdgesIndices=optLiftedIndices;

        size_t yIndex=getIndexInYBase(centralNodeID);
        y[yIndex]=bestSolutionIndex;

        size_t j = 0;
        for (; j < optLiftedIndices.size(); ++j) {
            size_t liftedIndex=optLiftedIndices[j];
            y[yIndex+j+1]=liftedIndex;

        }
        for(;j<maxTimeGap;j++){
            y[yIndex+j+1]=numberOfNodes;
        }

        return bestSolutionValue;
//        myStr.optBaseIndex=bestSolutionIndex;
//        myStr.optValue=bestSolutionValue;


    }



    void copyYToX(double* x,size_t* y) const{
        const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
        const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
        for (size_t i = 0; i < numberOfNodes; ++i) {
            size_t indexInY=getIndexInYBase(i);
            size_t baseEdgeIndex=y[indexInY];
            size_t numberOfBase=baseGraph.getNumberOfEdgesFromVertex(i);
            if(!isOutFlow) numberOfBase=baseGraph.getNumberOfEdgesToVertex(i);
            size_t firstBaseIndex=getIndexInWIBase(i,0);
            size_t firstILiftedndexInX=getIndexInWILifted(i,0);
            size_t numberOfLifted=liftedGraph.getNumberOfEdgesFromVertex(i);
            if(!isOutFlow) numberOfLifted=liftedGraph.getNumberOfEdgesToVertex(i);
            if(baseEdgeIndex==numberOfBase){//inactive node
                x[i]=0;
                for (size_t j = 0; j < numberOfBase; ++j) {
                    x[firstBaseIndex+j]=0;
                }
                for (size_t j = 0; j < numberOfLifted; ++j) {
                    x[firstILiftedndexInX+j]=0;
                }
            }
            else{
                x[i]=1;
                for (size_t j = 0; j < numberOfBaseEdges; ++j) {
                    if(j!=baseEdgeIndex){
                        x[firstBaseIndex+j]=0;
                    }
                    else{
                        x[firstBaseIndex+j]=1;
                    }
                }
                size_t indexInY=getIndexInYLifted(i);
                size_t firstIndexInX=getIndexInWILifted(i,0);
                size_t nextOptLiftedIndex=y[indexInY];
                for (size_t j = 0; j < numberOfLifted; ++j) {
                    if(j==nextOptLiftedIndex){
                        x[firstIndexInX+j]=1;
                        indexInY++;
                        if(indexInY<maxTimeGap) nextOptLiftedIndex=y[indexInY];
                        else nextOptLiftedIndex=numberOfNodes;
                    }
                    else{
                        x[firstIndexInX+j]=0;
                    }

                }

            }
        }




    }





private:
    size_t getVertexToReach()const{
        if(isOutFlow){
            return pInstance->getTerminalNode();
        }
        else{
            return pInstance->getSourceNode();
        }
    }

    const lifted_disjoint_paths::LdpInstance* pInstance;
    size_t numberOfNodes;
    size_t numberOfLiftedEdges;
    size_t numberOfBaseEdges;
    bool isOutFlow;
    std::vector<std::map<size_t,size_t>> nodeIDToIndex;
    size_t maxTimeGap;
    size_t yLength;
    size_t xLength;
    std::vector<std::vector<size_t>> traverseOrders; //Obtained from reachability structure

};

double minInOutFlow(double* wi, FWMAP::YPtr _y, FWMAP::TermData term_data) // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
{
  LdpProblem* ldpProblem = (LdpProblem*) term_data;
  size_t* y = (size_t*) _y;
  memset(y, 0, ldpProblem->getYLength()*sizeof(char));
  double optValue=0;
  for (size_t i = 0; i < ldpProblem->getNumberOfNodes(); ++i) {
      optValue+=ldpProblem->topDownMethod(i,wi,y);

  }


  return optValue;
}

static void copyYtoX(double* xi, FWMAP::YPtr _y, FWMAP::TermData term_data)
{
    LdpProblem* ldpProblem = (LdpProblem*) term_data;
    size_t* y = (size_t*) _y;

    ldpProblem->copyYToX(xi,y);
}

static double dot_product_fn_test(double* wi, FWMAP::YPtr _y, FWMAP::TermData term_data)
{
  sub_problem_test* sp = (sub_problem_test*) term_data;
  size_t* y = (size_t*) _y;
  assert(*y == 0 || *y == 1);
  return wi[*y];
}

int main()
{
  FWMAP s (2, 2, max_fn_test, copy_fn_test, dot_product_fn_test);//int d, int n, MaxFn max_fn, CopyFn copy_fn, DotProductFn dot_product_fn;

  sub_problem_test sp1;
  sp1.cost = {9.0,10.0};
  sub_problem_test sp2;
  sp2.cost = {20.0,0.0};

  std::array<int,2> mapping = {0,1};
  s.SetTerm(0, &sp1, 2, &mapping[0], sizeof(size_t));
  s.SetTerm(1, &sp2, 2, &mapping[0], sizeof(size_t));

  s.init();
  double cost;
  for(size_t i=0; i<100; ++i) {
    cost = s.do_descent_step();
  }
  //double cost = s.Solve();

  double* lambda1 = s.GetLambda(0);
  double* lambda2 = s.GetLambda(1);

  std::cout << "\n\nsolution = ";
  for(int i=0; i<2; ++i) { std::cout << lambda1[i] << ","; }
  std::cout << "\n";
  for(int i=0; i<2; ++i) { std::cout << lambda2[i] << ","; }
  std::cout << "\n";

  assert(std::abs(cost - 10.0) < 10e-6);
  return 0;
}

} //End of namespace LPMP

#endif // LDP_FW_HXX
