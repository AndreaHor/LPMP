#ifndef LDP_FW_HXX
#define LDP_FW_HXX


#include "FW-MAP.h"
#include <cassert>
#include <array>
#include <iostream>
#include <cmath>
#include <vector>
#include "ldp_instance.hxx"
#include <algorithm>

namespace LPMP{


class LdpProblem{
public:

    LdpProblem(const lifted_disjoint_paths::LdpInstance* _pInstance,bool _isOutFlow);
    void copyYToXi(double* x,size_t* y) const;
    double dotProduct(double* wi,size_t* y) const;

    size_t getNumberOfNodes() const{
        return numberOfInnerNodes;
    }

    bool getIsOutFlow()const{
        return isOutFlow;
    }

    size_t getYLength()const{
        return yLength;
    }

    size_t getXLength()const{
        return xLength;
    }

    double topDownMethod(size_t centralNodeID,double* wi,size_t* y);

    void initVectorForMapping(std::vector<int>& vectorForMapping );
     bool debugOutputs=false;


private:
    size_t getVertexToReach()const{
        if(isOutFlow) return numberOfInnerNodes+1;
        else return numberOfInnerNodes;
    }

    size_t getIndexInYLifted(size_t centralNodeID)const;

    size_t getIndexInYBase(size_t centralNodeID) const;

    //Can be called just to get the index of the first neighbor and then iterate without calling this function always
    size_t getIndexInWILifted(size_t centralNodeID)const;

    //Can be called just to get the index of the first neighbor and then iterate without calling this function always
    size_t getIndexInWIBase(size_t centralNodeID)const;

    size_t getIndexInWINode(size_t centralNodeID);


    size_t nodeIDtoNeighborIndex(const size_t& centralNodeID, const size_t& neighborID);





    const lifted_disjoint_paths::LdpInstance* pInstance;
    size_t numberOfInnerNodes; //without s and t
    size_t numberOfAllNodes;
    size_t numberOfLiftedEdges;
    size_t numberOfLocalBaseEdges;
    bool isOutFlow;
    std::vector<std::map<size_t,size_t>> nodeIDToIndex;
    size_t maxTimeGap;
    size_t yLength;
    size_t xLength;
    std::vector<std::vector<size_t>> traverseOrders; //Obtained from reachability structure


};

double minInOutFlowLDP(double* wi, FWMAP::YPtr _y, FWMAP::TermData term_data);

static void copyYtoXLDP(double* xi, FWMAP::YPtr _y, FWMAP::TermData term_data)
{
    LdpProblem* ldpProblem = (LdpProblem*) term_data;
    size_t* y = (size_t*) _y;

    ldpProblem->copyYToXi(xi,y);
}

static double dotProductLDP(double* wi, FWMAP::YPtr _y, FWMAP::TermData term_data)
{
    LdpProblem* ldpProblem = (LdpProblem*) term_data;
    size_t* y = (size_t*) _y;

    // std::vector<double> zeroLambda(ldpProblem->getXLength());
    double dp=ldpProblem->dotProduct(wi,y);
     //double dp=ldpProblem->dotProduct(zeroLambda.data(),y);
     return dp;
}



} //End of namespace LPMP

#endif // LDP_FW_HXX
