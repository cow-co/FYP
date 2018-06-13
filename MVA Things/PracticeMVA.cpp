#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

/**
 * Executes an adaptive BDT and a Neural Network on the given root files
 * @param  argc Number of command-line arguments
 * @param  argv The command-line arguments themselves
 * @return      0 for success, -1 for failure.
 */
int main(int argc, char* argv[])
{
  // if(argc != 4)
  // {
  //  std::cout << "Input form ./programname <sigfile> <bgfile> <outputfilename> <Variable1> <Variable2> ... <VariableN> " << std::endl;
  //  return -1;
  // }

  TMVA::Tools::Instance();  //This loads the library, singleton-style.
  std::map<std::string, int> toUse; //List of tools and whether we will use them.

  toUse["MLP"] = 1; //lol MLP.  This is a Neural network.
  toUse["BDT"] = 1; //Adaptive Boosted Decision Tree.

  std::cout << "Beginning Classification using MLP Neural Network and BDT." << std::endl;

  TString outFileName = "OutputFiles/";
  outFileName.Append(argv[3]);
  outFileName.Append(".root");
  TFile* outFile = TFile::Open(outFileName, "RECREATE");

  //Set up the tool Factory
  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification", outFile,
                                             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification");

  int i;
  std::string cutString = "";

  //for loop to add all variables to be analysed
  for(i = 4; i<argc; i+=2)
    {TString VariableName = argv[i];
     TString VariableString = argv[i];
     double max = atof(argv[i+1]);
     std::ostringstream strs;
     strs << max;
     std::string str = strs.str();
     factory->AddVariable(VariableName, VariableString,"Units", 'F');
     std::cout << "Variable " << VariableString << " is being analysed."<< std::endl;
     if(i==4){
       cutString.append(VariableName);
       cutString.append("<");
       cutString.append(str);
     }
     else{
     cutString.append(" && ");
     cutString.append(VariableName);
     cutString.append("<");
     cutString.append(str);}
}
  TString cutTString(cutString);
  TCut cut(cutTString);
  //Grab the names of the signal and background files from the cmd args.
  TString sigFilename = argv[1];
  TString bgFilename = argv[2];

  std::cout << sigFilename << ", " << bgFilename << std::endl;

  //Open the signal and BG files.
  TFile* sigFile = TFile::Open(sigFilename);
  TFile* bgFile = TFile::Open(bgFilename);

  TDirectoryFile* sigTupleDirectory;
  sigFile->GetObject("Incl_Tuple;1", sigTupleDirectory);

  TDirectoryFile* bgTupleDirectory;
  bgFile->GetObject("Incl_Tuple;1", bgTupleDirectory);

  if(!sigTupleDirectory || !bgTupleDirectory)
    {
      std::cout << "RIP" << std::endl;

      return -1;
    }

  TTree* sigTree = (TTree*)sigTupleDirectory->Get("DecayTree;1");
  TTree* bgTree = (TTree*)bgTupleDirectory->Get("DecayTree;1");

  Double_t sigEventWeight = 1.0;
  Double_t bgEventWeight = 1.0;

  factory->AddSignalTree(sigTree, sigEventWeight);
  factory->AddBackgroundTree(bgTree, bgEventWeight);

  factory->PrepareTrainingAndTestTree(cut, cut,
                                      "SplitMode=Random:NormMode=NumEvents:!V");

  factory->BookMethod(TMVA::Types::kMLP, "MLP",
                      "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator")\
    ; //book the Neural Network
  factory->BookMethod(TMVA::Types::kBDT, "BDT",
                      "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedB\
oost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"); //book the BDT

  //Time to train and test the MVAs!
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outFile->Close();

  std::cout << "Done!" << std::endl;

  delete factory;
}
