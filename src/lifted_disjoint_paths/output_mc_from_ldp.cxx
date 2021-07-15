#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "tclap/CmdLine.h"
#include <sstream>

using namespace LPMP;

int main(int argc, char** argv)
{


    TCLAP::CmdLine cmd_(std::string("Command line options for transforming LDP to MC"), ' ', "0.0.1");

    TCLAP::ValueArg<std::string> inputFileArg_("i","inputFile","file from which to read problem instance",true,"","file name",cmd_);
    TCLAP::ValueArg<std::string> outputFileArg_("o","outputFile","file to save result",true,"","file name",cmd_);
    TCLAP::ValueArg<char> separatorArg_("s","separatingCharacter","separating character for output",false,' ',"separating character",cmd_);
    TCLAP::ValueArg<double> mustCutCostArg_("","mustCutCost","must cut cost",false,-100,"must cut cost",cmd_);

    cmd_.parse(argc,argv);
    std::string inputFileName = inputFileArg_.getValue();
    std::string outputFileName = outputFileArg_.getValue();


    LPMP::lifted_disjoint_paths::LdpParameters<> configParams(inputFileName);

    LPMP::CompleteStructure<> completeStructure(configParams);
    std::vector<std::array<size_t,2>> edgeList=completeStructure.getEdgeList();
    const auto& vg=completeStructure.getVertexGroups();

    const LdpDirectedGraph& completeGraph=completeStructure.myCompleteGraph;


    std::string fileForSave=outputFileName;
    std::ofstream file;
    file.open(fileForSave, std::ofstream::out);
    if(!file.is_open()) {
        throw std::runtime_error("could not open file " + fileForSave);
    }


    char separator=separatorArg_.getValue();
    double withinFramePenalty=mustCutCostArg_.getValue();

    assert(withinFramePenalty<=0);

    //All scorese with minus because multicut has an oposite meaning for positive and negative edges

    for (size_t i = 0; i < completeGraph.getNumberOfVertices(); ++i) {
        //TODO here come edges within layer with a high score
        size_t groupIndex=vg.getGroupIndex(i);
        const std::vector<size_t>& groupVertices=vg.getGroupVertices(groupIndex);
        for(auto& v:groupVertices){
            if(v>i){
                file<<i<<separator<<v<<withinFramePenalty<<"\n";
            }
        }

        for (auto iter=completeGraph.forwardNeighborsBegin(i);iter!=completeGraph.forwardNeighborsEnd(i);iter++) {
            file<<i<<separator<<(iter->first)<<separator<<-(iter->second)<<"\n";

        }
    }

    //TODO which parameters in the parameter file are used? Create an example param file.

     //file << solution.str();
     file.close();


    //ldpInstance.setOutputFileName(outputFileName);

}
