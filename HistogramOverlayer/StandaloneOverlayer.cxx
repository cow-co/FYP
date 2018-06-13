#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[])
{
  if(argc != 9 )
    {
      std::cout << "Usage: ./StandaloneOverlayer <signal file> <bg file> <variable name> <min> <max> <numBins> <outputname> <verbose>" << std::endl;

      if(argc!= 8)
        {
          return -1;
        }
    }

  TString sigFilename = argv[1]; //zeroth arg is the program name
  TString bgFilename = argv[2];

  TFile* sigFile = TFile::Open(sigFilename);
  TFile* bgFile = TFile::Open(bgFilename);

  TDirectoryFile* sigDir;
  TDirectoryFile* bgDir;
  sigFile->GetObject("Incl_Tuple;1", sigDir);
  bgFile->GetObject("Incl_Tuple;1", bgDir);

  if(!sigDir || !bgDir)
    {
      std::cout << "RIP dir" << std::endl;

      return -1;
    }

  TTree* sigTree;
  TTree* bgTree;
  sigTree = (TTree*)sigDir->Get("DecayTree;1");
  bgTree = (TTree*)bgDir->Get("DecayTree;1");

  if(!sigTree || !bgTree)
    {
      std::cout << "RIP tree" << std::endl;

      return -1;
    }

  double min = atof(argv[4]);
  double max = atof(argv[5]);
  int numBins = atoi(argv[6]);
  TH1D varSig("Signal"," " , numBins, min, max);
  TH1D varBG("Background"," " , numBins, min, max);
  TString verb = argv[8];
  if(verb == "-v")
    {
      std::string VarName;
      std::cout << "Please enter the variable name in math mode format e.g. B_{0} MM \n" << std::endl;
      std::getline(std::cin,VarName);
      TString varTString(VarName);
      std::string Title;
      std::cout << "Please enter a title for the histogram:" << std::endl;
      std::getline(std::cin, Title);
      TString TTitle(Title);
      varSig.SetTitle(TTitle);
      varSig.GetXaxis()->SetTitle(varTString);
      varSig.GetYaxis()->SetTitle("Candidates");
      varSig.SetStats(false);
      varBG.SetStats(false);
    }
  else{
    TString varName(argv[3]);
    varSig.GetXaxis()->SetTitle(varName);
    varSig.GetYaxis()->SetTitle("Candidates");
    varSig.SetStats(false);
    varBG.SetStats(false);

  }
  //Analyse the signal first
  Long64_t numEntries = sigTree->GetEntriesFast();
  TLeaf* varLeaf = sigTree->GetLeaf(argv[3]);

  if(!varLeaf)
    {
      std::cout << "RIP leaf" << std::endl;

      return -1;
    }

  for (Long64_t jEntry = 0; jEntry < numEntries; ++jEntry)
    {
      int numBytes = varLeaf->GetBranch()->GetEntry(jEntry);

      if(numBytes == 0)
        {
          std::cout << "RIP entry" << std::endl;

          return -1;
        }

      varSig.Fill(varLeaf->GetValue());
    }

  varSig.Scale(1.0/varSig.Integral(min, max)); //Normalise the area of the histogram
  varSig.SetLineColor(1);

  numEntries = bgTree->GetEntriesFast();
  varLeaf = bgTree->GetLeaf(argv[3]);
  
  //Grab the origin vertex positions for the two muons, in order to strip out signal within our background
  TLeaf* mmOrigVtxXLeaf = bgTree->GetLeaf("muminus_TRUEORIGINVERTEX_X");  
  TLeaf* mmOrigVtxYLeaf = bgTree->GetLeaf("muminus_TRUEORIGINVERTEX_Y");  
  TLeaf* mmOrigVtxZLeaf = bgTree->GetLeaf("muminus_TRUEORIGINVERTEX_Z");
  
  TLeaf* mpOrigVtxXLeaf = bgTree->GetLeaf("muplus_TRUEORIGINVERTEX_X");
  TLeaf* mpOrigVtxYLeaf = bgTree->GetLeaf("muplus_TRUEORIGINVERTEX_Y");
  TLeaf* mpOrigVtxZLeaf = bgTree->GetLeaf("muplus_TRUEORIGINVERTEX_Z");
  
  double mmOrigVtxX;
  double mmOrigVtxY;
  double mmOrigVtxZ;
  
  double mpOrigVtxX;
  double mpOrigVtxY;
  double mpOrigVtxZ;

  if(!varLeaf)
    {
      std::cout << "RIP leaf" << std::endl;

      return -1;
    }

  for(Long64_t jEntry = 0; jEntry < numEntries; ++jEntry)
    {
      int nBytes = varLeaf->GetBranch()->GetEntry(jEntry);

      if(nBytes <= 0)
        {
          std::cout << "RIP Entry" << std::endl;

          return -1;
        }
      
      //Grab the origin vertex positions and check if they are equal; if so, chuck the candidate out.
      mmOrigVtxX = mmOrigVtxXLeaf->GetValue();
      mmOrigVtxY = mmOrigVtxYLeaf->GetValue();
      mmOrigVtxZ = mmOrigVtxZLeaf->GetValue();
      
      mpOrigVtxX = mpOrigVtxXLeaf->GetValue();
      mpOrigVtxY = mpOrigVtxYLeaf->GetValue();
      mpOrigVtxZ = mpOrigVtxZLeaf->GetValue();
      
      if(mmOrigVtxX == mpOrigVtxX &&
         mmOrigVtxY == mpOrigVtxY &&
         mmOrigVtxZ == mpOrigVtxZ)
        {
          continue;
        }
      
      else
        {
          varBG.Fill(varLeaf->GetValue());
        }
    }

  varBG.Scale(1.0/varBG.Integral(min, max));
  varBG.SetLineColor(2);

  TString outFilename = "Plots/";
  outFilename.Append(argv[7]);
  outFilename.Append(".pdf");

  TCanvas c;
  varSig.Draw();
  varBG.Draw("SAME");
  c.SaveAs(outFilename);
  return 0;
}
