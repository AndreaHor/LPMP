#include"lifted_disjoint_paths/ldp_fw.hxx"

namespace LPMP{

LdpProblem:: LdpProblem(const lifted_disjoint_paths::LdpInstance* _pInstance,bool _isOutFlow)
{
    pInstance=_pInstance;
    numberOfInnerNodes=pInstance->getNumberOfVertices()-2;
    numberOfAllNodes=numberOfInnerNodes+2;
    numberOfLocalBaseEdges=pInstance->getMyGraph().getNumberOfEdges();
    assert(numberOfLocalBaseEdges>numberOfInnerNodes);
    numberOfLocalBaseEdges-=numberOfInnerNodes; //Edges from s (to t) are not present in outflow (inflow)
    numberOfLiftedEdges=pInstance->getMyGraphLifted().getNumberOfEdges();
    isOutFlow=_isOutFlow;
    maxTimeGap=std::max(pInstance->parameters.getMaxTimeLifted(),pInstance->parameters.getMaxTimeLifted());
    yLength=numberOfInnerNodes*(maxTimeGap+1);
    xLength=numberOfInnerNodes+numberOfLocalBaseEdges+numberOfLiftedEdges;

    nodeIDToIndex=std::vector<std::map<size_t,size_t>>(numberOfInnerNodes);
    if(isOutFlow){
        for (size_t i = 0; i < numberOfInnerNodes; ++i) {
            size_t indexCounter=0;
            const LdpDirectedGraph::edge* edgeIt=pInstance->getMyGraphLifted().forwardNeighborsBegin(i);
            for (;edgeIt!=pInstance->getMyGraphLifted().forwardNeighborsEnd(i);edgeIt++) {
                size_t nodeID=edgeIt->first;
                nodeIDToIndex[i][nodeID]=indexCounter;
                assert(pInstance->getMyGraphLifted().getForwardEdgeVertex(i,indexCounter)==nodeID);
                indexCounter++;
            }
        }
    }
    else{
        for (size_t i = 0; i < numberOfInnerNodes; ++i) {
            size_t indexCounter=0;
            const LdpDirectedGraph::edge* edgeIt=pInstance->getMyGraphLifted().backwardNeighborsBegin(i);
            for (;edgeIt!=pInstance->getMyGraphLifted().backwardNeighborsEnd(i);edgeIt++) {
                size_t nodeID=edgeIt->first;
                nodeIDToIndex[i][nodeID]=indexCounter;
                assert(pInstance->getMyGraphLifted().getBackwardEdgeVertex(i,indexCounter)==nodeID);
                indexCounter++;
            }
        }

    }

    traverseOrders=std::vector<std::vector<size_t>>(numberOfInnerNodes);
    if(!isOutFlow){
        for (size_t i = 0; i < numberOfInnerNodes; ++i) {
            const std::set<size_t>& reachableVertices=pInstance->reachableFromVertex(i);
            for (auto it=reachableVertices.begin();it!=reachableVertices.end();it++) {
                if(*it<numberOfInnerNodes){
                    traverseOrders[*it].push_back(i);
                }
            }
        }
    }
    else{
        assert(numberOfInnerNodes>0);
        for (size_t i = 0; i <numberOfInnerNodes; ++i) {
            const std::set<size_t>& reachableVertices=pInstance->reachableFromVertex(i);
            for (auto it=reachableVertices.rbegin();it!=reachableVertices.rend();it++) {
                //std::cout<<size_t(*it)<<std::endl;
                //assert(<numberOfNodes);
                if(*it<numberOfInnerNodes){
                    traverseOrders[size_t(i)].push_back(*it);
                }
            }
        }
    }
}

size_t LdpProblem::getIndexInYLifted(size_t centralNodeID)const{
    assert(centralNodeID<numberOfInnerNodes);
    size_t index=centralNodeID*(maxTimeGap+1)+1;
    return index;
}

size_t LdpProblem::getIndexInYBase(size_t centralNodeID) const{
     assert(centralNodeID<numberOfInnerNodes);
    size_t index=centralNodeID*(maxTimeGap+1);
    return index;
}

//Can be called just to get the index of the first neighbor and then iterate without calling this function always
size_t LdpProblem::getIndexInWILifted(size_t centralNodeID)const{
    size_t index=0;
    if(isOutFlow){
        //index=pInstance->getMyGraphLifted().getIndexForward(centralNodeID,0);
        index=pInstance->getMyGraph().getOffsetForward(centralNodeID);
    }
    else{
        //index=pInstance->getMyGraphLifted().getIndexBackward(centralNodeID,0);
        index=pInstance->getMyGraph().getOffsetBackward(centralNodeID);
    }
    return index+numberOfInnerNodes+numberOfLocalBaseEdges;
}

//Can be called just to get the index of the first neighbor and then iterate without calling this function always
size_t LdpProblem::getIndexInWIBase(size_t centralNodeID)const{
    size_t index=0;
    if(isOutFlow){
        index=pInstance->getMyGraph().getIndexForward(centralNodeID,0);
    }
    else{
        index=pInstance->getMyGraph().getIndexBackward(centralNodeID,0);
    }
    return index+numberOfInnerNodes;
}

size_t LdpProblem::getIndexInWINode(size_t centralNodeID){
    return centralNodeID;
}


size_t LdpProblem::nodeIDtoNeighborIndex(const size_t& centralNodeID, const size_t& neighborID){
    assert(centralNodeID<numberOfInnerNodes);
    auto f=nodeIDToIndex[centralNodeID].find(neighborID);
    if(f!=nodeIDToIndex[centralNodeID].end()){
        return f->second;
    }
    else{
        throw std::invalid_argument("Required index of non-existing neighbor of node "+std::to_string(centralNodeID));
    }
}


double LdpProblem::dotProduct(double* wi,size_t* y) const{
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
    double product=0;

    for (size_t i = 0; i < numberOfInnerNodes; ++i) {
        size_t indexInYBase=getIndexInYBase(i);
        size_t baseEdgeIndex=y[indexInYBase];
        size_t numberOfBase=baseGraph.getNumberOfEdgesFromVertex(i);
        if(!isOutFlow) numberOfBase=baseGraph.getNumberOfEdgesToVertex(i);
        size_t firstBaseIndexInX=getIndexInWIBase(i);
        size_t firstLiftedIndexInX=getIndexInWILifted(i);
        size_t numberOfLifted=liftedGraph.getNumberOfEdgesFromVertex(i);
        if(!isOutFlow) numberOfLifted=liftedGraph.getNumberOfEdgesToVertex(i);
        if(baseEdgeIndex<numberOfBase){//active node
            product+=pInstance->getVertexScore(i);
            assert(firstBaseIndexInX+baseEdgeIndex<xLength);
            product+=wi[firstBaseIndexInX+baseEdgeIndex];

            size_t indexInYLifted=getIndexInYLifted(i);
            size_t optLiftedCounter=0;
           // size_t firstIndexInX=getIndexInWILifted(i);
            size_t nextOptLiftedIndex=y[indexInYLifted+optLiftedCounter];
            while(nextOptLiftedIndex<numberOfInnerNodes){
                assert(firstLiftedIndexInX+nextOptLiftedIndex<xLength);
                product+=wi[firstLiftedIndexInX+nextOptLiftedIndex];
                optLiftedCounter++;

                if(optLiftedCounter<maxTimeGap){
                    nextOptLiftedIndex=y[indexInYLifted+optLiftedCounter];
                }
                else break;

            }
        }
    }
    return product;

}



void LdpProblem::copyYToXi(double* x,size_t* y) const{
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
    for (size_t i = 0; i < numberOfInnerNodes; ++i) {
        size_t indexInYBase=getIndexInYBase(i);
        size_t activeBaseEdgeIndex=y[indexInYBase];
        size_t numberOfBase=baseGraph.getNumberOfEdgesFromVertex(i);
        if(!isOutFlow) numberOfBase=baseGraph.getNumberOfEdgesToVertex(i);
        size_t firstBaseIndexInX=getIndexInWIBase(i);
        size_t firsILiftedIndexInX=getIndexInWILifted(i);
        size_t numberOfLifted=liftedGraph.getNumberOfEdgesFromVertex(i);
        if(!isOutFlow) numberOfLifted=liftedGraph.getNumberOfEdgesToVertex(i);
        if(activeBaseEdgeIndex==numberOfBase){//inactive node
            x[i]=0;
            for (size_t j = 0; j < numberOfBase; ++j) {
                x[firstBaseIndexInX+j]=0;
            }
            for (size_t j = 0; j < numberOfLifted; ++j) {
                x[firsILiftedIndexInX+j]=0;
            }
        }
        else{
            x[i]=1;
            for (size_t j = 0; j < numberOfBase; ++j) {
                if(j!=activeBaseEdgeIndex){
                    x[firstBaseIndexInX+j]=0;
                }
                else{
                    x[firstBaseIndexInX+j]=1;
                }
            }
            size_t indexInY=getIndexInYLifted(i);
            size_t activeLiftedCounter=0;
            //size_t firstIndexInX=getIndexInWILifted(i);
            size_t nextOptLiftedIndex=y[indexInY];
            for (size_t j = 0; j < numberOfLifted; ++j) {
                if(j==nextOptLiftedIndex){
                    x[firsILiftedIndexInX+j]=1;
                    activeLiftedCounter++;
                    if(activeLiftedCounter<maxTimeGap) nextOptLiftedIndex=y[indexInY+activeLiftedCounter];
                    else nextOptLiftedIndex=numberOfInnerNodes;
                }
                else{
                    x[firsILiftedIndexInX+j]=0;
                }

            }

        }
    }

}



double LdpProblem::topDownMethod(size_t centralNodeID,double* wi,size_t* y){
    //TODO: Fix extraction of edge cost from instance. They must be divided by 2 unless they are edges from s or to t!

    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
    if(centralNodeID==4496){
        std::cout<<"weird case"<<std::endl;
    }

    assert(centralNodeID<numberOfInnerNodes);
    for(size_t v:traverseOrders[centralNodeID]){  //Initialize the structures
        pInstance->sncTDStructure[v]=0;
        pInstance->sncNeighborStructure[v]=getVertexToReach();
        pInstance->isBSF[v]=0;

    }
    size_t numberOfLiftedNeighbors=0;
    size_t firstIndexInWiLifted=getIndexInWILifted(centralNodeID);

    if(isOutFlow){  //Prefill sncTDStrucute with lifted edges costs
        numberOfLiftedNeighbors=liftedGraph.getNumberOfEdgesFromVertex(centralNodeID);
        size_t counter=0;
        for (auto it=liftedGraph.forwardNeighborsBegin(centralNodeID);it!=liftedGraph.forwardNeighborsEnd(centralNodeID);it++) {
            double origCost=it->second*0.5;
            pInstance->sncTDStructure[it->first]=origCost+wi[firstIndexInWiLifted+counter];
            counter++;
        }
    }
    else{
        numberOfLiftedNeighbors=liftedGraph.getNumberOfEdgesToVertex(centralNodeID);
        size_t counter=0;
        for (auto it=liftedGraph.backwardNeighborsBegin(centralNodeID);it!=liftedGraph.backwardNeighborsEnd(centralNodeID);it++) {
            double origCost=it->second*0.5;
            pInstance->sncTDStructure[it->first]=origCost+wi[firstIndexInWiLifted+counter];
            counter++;
        }
    }



    double bsfValue=std::numeric_limits<double>::max();
    size_t i=0;
    for (; i < traverseOrders[centralNodeID].size(); ++i) {//The last node in the traverse order should be the central node itself. Nodes s and t are excluded.

        size_t currentNode=traverseOrders[centralNodeID][i];

        size_t boundaryRelevantVertex=*traverseOrders[centralNodeID].begin();

        double bestDescValue=0;
        size_t bestDescVertexID=getVertexToReach();

        //Search for best descendant
        if(centralNodeID==4496&&(currentNode==5273||currentNode==5304)){
            std::cout<<"weird case"<<std::endl;
        }
        if(isOutFlow){
            const LdpDirectedGraph::edge* vertexIt=baseGraph.forwardNeighborsBegin(currentNode);
            const LdpDirectedGraph::edge* end=baseGraph.forwardNeighborsEnd(currentNode);


            for (;vertexIt!=end;vertexIt++) {

                size_t desc=vertexIt->first;

                if(desc>=numberOfInnerNodes) continue;

                if(desc<=boundaryRelevantVertex){ //First node in traverse order is the one with the highest ID within time gap

                    assert(pInstance->isReachable(centralNodeID,desc));
                    double value=pInstance->sncTDStructure[desc];
                    if(bestDescValue>value){
                        bestDescValue=value;
                        bestDescVertexID=desc;

                    }
                    if(pInstance->isBSF[desc]){ //All neighbors have higher ID, all nodes with higher ID have worse value
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

                if(desc>=numberOfInnerNodes){
                    if(vertexIt==begin){
                        doSearch=false;
                    }
                    else{
                        vertexIt--;
                    }
                }
                else{

                    if(desc>=boundaryRelevantVertex){//First in the traverse order is the most distant interesting node


                        assert(pInstance->isReachable(desc,centralNodeID));
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

    double nodeCost=pInstance->getVertexScore(centralNodeID)+wi[getIndexInWINode(centralNodeID)];
    std::vector<double> solutionCosts(nodeNotActive+1);
    solutionCosts[nodeNotActive]=0;
    size_t counter=0;
    size_t firstBaseIndexInW=getIndexInWIBase(centralNodeID);
    for (; baseNeighborsIt!=baseNeighborsEnd; baseNeighborsIt++) {

        double baseCost=baseNeighborsIt->second;
        size_t neighborID=baseNeighborsIt->first;
        if(neighborID<numberOfInnerNodes) baseCost*=0.5;
        size_t index=firstBaseIndexInW+counter;
        baseCost+=wi[index];


        double valueToAdd=0;
        if(neighborID<numberOfInnerNodes){
            valueToAdd=pInstance->sncTDStructure[neighborID];
        }

        double value=baseCost+nodeCost+valueToAdd;

        solutionCosts[counter]=value;
        if(value<bestSolutionValue){
            bestSolutionValue=value;
            bestSolutionIndex=counter;
            bestSolutionID=neighborID;
        }

        counter++;
    }

    std::vector<size_t> optLiftedIndices;
    if(bestSolutionIndex!=nodeNotActive){
        //            assert(bestSolutionID<numberOfNodes);
        while(bestSolutionID<numberOfInnerNodes){
            if((isOutFlow&&pInstance->existLiftedEdge(centralNodeID,bestSolutionID))||(!isOutFlow&&pInstance->existLiftedEdge(bestSolutionID,centralNodeID))){
                size_t liftedIndex=nodeIDToIndex[centralNodeID][bestSolutionID];
                optLiftedIndices.push_back(liftedIndex);
            }
            bestSolutionID=pInstance->sncNeighborStructure[bestSolutionID];
        }
    }
    //    OptimalSolution optSolution;
    //    optSolution.optValue=bestSolutionValue;
    //    optSolution.baseEdgeIndex=bestSolutionIndex;
    //    optSolution.liftedEdgesIndices=optLiftedIndices;

    size_t yIndex=getIndexInYBase(centralNodeID);
    y[yIndex]=bestSolutionIndex;

    assert(optLiftedIndices.size()<maxTimeGap);
    size_t j = 0;
    for (; j < optLiftedIndices.size(); ++j) {
        size_t liftedIndex=optLiftedIndices[j];
        y[yIndex+j+1]=liftedIndex;

    }
    for(;j<maxTimeGap;j++){
        y[yIndex+j+1]=numberOfInnerNodes;
    }

    return bestSolutionValue;
    //        myStr.optBaseIndex=bestSolutionIndex;
    //        myStr.optValue=bestSolutionValue;


}


void LdpProblem::initVectorForMapping(std::vector<int>& vectorForMapping ){

    assert(vectorForMapping.size()==xLength);

    if(!isOutFlow){
        const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
        const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
        for(size_t i=0;i<numberOfInnerNodes;i++){
            vectorForMapping[i]=int(i);
        }
        size_t counter=numberOfInnerNodes;

        for(size_t i=0;i<numberOfInnerNodes;i++){
            auto it=baseGraph.backwardNeighborsBegin(i);
            size_t counter2=0;
            for(;it!=baseGraph.backwardNeighborsEnd(i);it++){
                size_t secondNode=it->first;
                size_t reverseNeighborIndex=it->reverse_neighbor_index;
                int index=int(baseGraph.getIndexForward(secondNode,0));
                index+=reverseNeighborIndex;
                assert(counter==baseGraph.getIndexBackward(i,counter2)+numberOfInnerNodes);
                assert(size_t(index)==baseGraph.getIndexForward(secondNode,reverseNeighborIndex));
                vectorForMapping[counter]=index+int(numberOfInnerNodes);
                if(it->first>numberOfInnerNodes){
                    std::cout<<"s or t node "<<it->first<<std::endl;
                }
                counter++;
                counter2++;
            }
        }

        assert(counter==numberOfInnerNodes+numberOfLocalBaseEdges);
        size_t addToIndex=2*numberOfInnerNodes+numberOfLocalBaseEdges; //Inner nodes must be added twice to compensate missing base edges to t in the inflow subproblem
        for(size_t i=0;i<numberOfInnerNodes;i++){
            auto it=liftedGraph.backwardNeighborsBegin(i);
            for(;it!=liftedGraph.backwardNeighborsEnd(i);it++){
                size_t secondNode=it->first;
                size_t reverseNeighborIndex=it->reverse_neighbor_index;
                int index=int(liftedGraph.getIndexForward(secondNode,0));
                index+=reverseNeighborIndex;
                assert(size_t(index)==liftedGraph.getIndexForward(secondNode,reverseNeighborIndex));
                assert(counter<xLength);
                assert(size_t(index)<xLength);
                vectorForMapping[counter]=index+int(addToIndex);
                counter++;
            }
        }
        assert(counter=numberOfInnerNodes+numberOfLocalBaseEdges+numberOfLiftedEdges);
    }
    else{
        //0...numberOfInnerNodes-1: correspond to node variable to node variable
        //numberOfInnerNodes...numberOfLocalBaseEdges+numberOfInnerNodes: correspond to variables of base edges from inner nodes,
        //edges from s (there are numberOfInnerNodes of them) are not present in outflow
        for (size_t i = 0; i < numberOfLocalBaseEdges+numberOfInnerNodes; ++i) {
            vectorForMapping[i]=int(i);
        }
        size_t addToFirst=numberOfLocalBaseEdges+numberOfInnerNodes;
        size_t addToSecond=2*numberOfInnerNodes+numberOfLocalBaseEdges;
        for (size_t i = 0; i < numberOfLiftedEdges; ++i) {
            assert(i+addToFirst<xLength);
            assert(i+addToSecond<xLength+numberOfInnerNodes);
            vectorForMapping[i+addToFirst]=int(i+addToSecond);
        }
    }




}



double minInOutFlowLDP(double* wi, FWMAP::YPtr _y, FWMAP::TermData term_data) // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
{
  LdpProblem* ldpProblem = (LdpProblem*) term_data;
  size_t* y = (size_t*) _y;
  memset(y, 0, ldpProblem->getYLength()*sizeof(size_t));
  double optValue=0;
  for (size_t i = 0; i < ldpProblem->getNumberOfNodes(); ++i) {
      optValue+=ldpProblem->topDownMethod(i,wi,y);

  }


  return optValue;
}




}
