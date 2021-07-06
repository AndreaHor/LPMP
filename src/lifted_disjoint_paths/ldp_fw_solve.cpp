#include"lifted_disjoint_paths/ldp_fw.hxx"
#include "tclap/CmdLine.h"
#include <sstream>
#include "FW-MAP.h"

using namespace LPMP;

int main(int argc, char** argv)
{


    std::cout<<"start"<<std::endl;
    std::cout<<"start 2"<<std::endl;

    TCLAP::CmdLine cmd_(std::string("Command line options for LDP FW"), ' ', "0.0.1");

    TCLAP::ValueArg<std::string> inputFileArg_("i","inputFile","file from which to read problem instance",true,"","file name",cmd_);
    TCLAP::ValueArg<std::string> outputFileArg_("o","outputFile","file to save result",true,"","file name",cmd_);

    cmd_.parse(argc,argv);
    std::string inputFileName = inputFileArg_.getValue();
    std::string outputFileName = outputFileArg_.getValue();

    //ProblemConstructorRoundingSolver<Solver<LP<lifted_disjoint_paths_FMC>,StandardTighteningVisitor>> solver(argc,argv);
    //I need a way how to obtain the command line arguments without this constructor

    //std::string inputFileName=solver.get_input_file();
    //std::string outputFileName=solver.getOutputFileName();
    //std::cout<<"output file name "<<outputFileName<<std::endl;

    LPMP::lifted_disjoint_paths::LdpParameters<> configParams(inputFileName);

    LPMP::CompleteStructure<> completeStructure(configParams);

    LPMP::lifted_disjoint_paths::LdpInstance ldpInstance(configParams,completeStructure);
    ldpInstance.setOutputFileName(outputFileName);


    LdpProblem ldpProblemOut(&ldpInstance,true);
    LdpProblem ldpProblemIn(&ldpInstance,false);

    std::vector<int> mappingIn(ldpProblemIn.getXLength());
    ldpProblemIn.initVectorForMapping(mappingIn);
    int* mappingDataIn=mappingIn.data();

    std::vector<int> mappingOut(ldpProblemOut.getXLength());
    ldpProblemOut.initVectorForMapping(mappingOut);
    int* mappingDataOut=mappingOut.data();

    int globalXLength=int(ldpProblemOut.getXLength()+ldpProblemOut.getNumberOfNodes());

    FWMAP s (globalXLength, 2, minInOutFlowLDP, copyYtoXLDP, dotProductLDP);//int d, int n, MaxFn max_fn, CopyFn copy_fn, DotProductFn dot_product_fn;

//TODO mapping!
  //std::array<int,2> mapping = {0,1}; //TODO map in variables to out variables

  //mapping should map variables of problem In to variables of problem Out
  //TODO: check the sizes of x and y of both problems, they should be the same.
  s.SetTerm(0, &ldpProblemIn, ldpProblemIn.getXLength(), mappingDataIn, ldpProblemIn.getYLength()*sizeof(size_t));
  s.SetTerm(1, &ldpProblemOut, ldpProblemOut.getXLength(), mappingDataOut, ldpProblemOut.getYLength()*sizeof(size_t));

  s.init();
  double cost=0;
  for(size_t i=0; i<100; ++i) {
    cost = s.do_descent_step();
    std::cout<<std::endl;
    std::cout<<"cost in "<<i<<"th step "<<cost<<std::endl;
  }
  //double cost = s.Solve();

  double* lambda1 = s.GetLambda(0);
  double* lambda2 = s.GetLambda(1);

//  std::cout << "\n\nsolution = ";
//  for(int i=0; i<ldpProblemIn.getXLength(); ++i) { std::cout << lambda1[i] << ","; }
//  std::cout << "\n";
//  for(int i=0; i<ldpProblemOut.getXLength(); ++i) { std::cout << lambda2[i] << ","; }
//  std::cout << "\n";

  //assert(std::abs(cost - 10.0) < 10e-6);
  return 0;
}


