#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Config.h"

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "TMVAClassification.h"

int TMVAClassification(std::string inputSname, std::string inputBname, std::string mycuts, std::string mycutb, 
                       std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage)
#int TMVAClassification(std::string inputSname, std::string inputBname, std::string mycuts, std::string mycutb, 
                       std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage, int nTrees)
{
  std::vector<std::string> methods;
  std::vector<int> stages;
  std::string outfname = mytmva::mkname(outputname, ptmin, ptmax, mymethod, stage, methods, stages);
  std::string outputstr = xjjc::str_replaceallspecial(outfname);
  if(ptmax < 0) { ptmax = 1.e+10; }

  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the
  // corresponding lines from .rootrc

  // Methods to be processed can be given as an argument; use format:
  //
  //     mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  //
  // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // ------>>
  if(mymethod != "")
    {
      for(std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
      for(auto& m : methods)
        {
          if(Use.find(m) != Use.end())
            { Use[m] = 1; std::cout <<"==> " << __FUNCTION__ << ": Registered method " << m << std::endl; }
          else
            { std::cout << "==> Abort " << __FUNCTION__ << ": error: unknown method " << m << "." << std::endl; continue; }
        }
    }
  // <<------

  // Select methods (don't look at this code - not of interest)
  // if (myMethodList != "") {
  //   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

  //   // std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
  //   for (UInt_t i=0; i<mlist.size(); i++) {
  //     std::string regMethod(mlist[i]);

  //     if (Use.find(regMethod) == Use.end()) {
  //       std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
  //       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
  //       std::cout << std::endl;
  //       return 1;
  //     }
  //     Use[regMethod] = 1;
  //   }
  // }

  // --------------------------------------------------------------------------------------------------

  // Here the preparation phase begins

  // Read training and test data
  // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
  //// TString fname = "./tmva_class_example.root";

  //// if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
  ////    gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");

  //// TFile *input = TFile::Open( fname );

  TFile* inputS = TFile::Open(inputSname.c_str());
  TFile* inputB = TFile::Open(inputBname.c_str());

  //// std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
  std::cout << "--- TMVAClassification       : Using input file: " << inputS->GetName() << " & "<< inputB->GetName() <<std::endl;

  // Register the training and test trees

  //// TTree* signalTree     = (TTree*)input->Get("TreeS");
  //// TTree* background     = (TTree*)input->Get("TreeB");

  TTree* background = (TTree*)inputB->Get("Dfinder/ntDkpi");
  background->AddFriend("hltanalysis/HltTree");
  background->AddFriend("hiEvtAnalyzer/HiTree");
  background->AddFriend("skimanalysis/HltTree");
  background->AddFriend("ppTracks/trackTree");
  background->AddFriend("l1MetFilterRecoTree/MetFilterRecoTree");
  background->AddFriend("zdcanalyzer/zdcdigi");
  background->AddFriend("particleFlowAnalyser/pftree");

  TTree* signal = (TTree*)inputS->Get("Dfinder/ntDkpi");
  signal->AddFriend("hltanalysis/HltTree");
  signal->AddFriend("hiEvtAnalyzer/HiTree");
  signal->AddFriend("skimanalysis/HltTree");
  signal->AddFriend("ppTracks/trackTree");
  signal->AddFriend("l1MetFilterRecoTree/MetFilterRecoTree");
  signal->AddFriend("zdcanalyzer/zdcdigi");
  signal->AddFriend("particleFlowAnalyser/pftree");

  // auto hlt = (TTree*)inputB->Get("hltanalysis/HltTree");
  // hlt->Print();
 
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  //// TString outfileName( "TMVA.root" );
  //// TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  xjjroot::mkdir(outfname.c_str());
  TFile* outf = TFile::Open(outfname.c_str(), "RECREATE");

  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outf,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  gSystem->Exec(Form("mkdir -p dataset/weights/%s", outputstr.c_str()));
  (TMVA::gConfig().GetIONames()).fWeightFileDir = Form("weights/%s", outputstr.c_str());

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  //// dataloader->AddVariable( "myvar1 := var1+var2", 'F' );
  //// dataloader->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
  //// dataloader->AddVariable( "var3",                "Variable 3", "units", 'F' );
  //// dataloader->AddVariable( "var4",                "Variable 4", "units", 'F' );

  // Variable
  std::string varinfo = "";
  TString VarSet = "";
  int nvar = 0;
  for(auto& s : stages)
    {
      dataloader->AddVariable(mytmva::varlist[s].var);
      if(mytmva::varlist[s].cutsign != "")
        { VarSet += (Form(":VarProp[%d]=", nvar)+mytmva::varlist[s].cutsign); }
      varinfo += (";"+mytmva::varlist[s].varname);
      std::cout << "==> " << __FUNCTION__ << ": Registered variable " << mytmva::varlist[s].var << std::endl;
      nvar++;
    }
  if(!nvar) { std::cout << "==> Abort " << __FUNCTION__ << ": no variable registered." << std::endl; return 2; }
  std::cout << "==> " << __FUNCTION__ << ": VarSet = " << VarSet << std::endl;

  // Spectator
  dataloader->AddSpectator("Dmass");

  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables

  //// dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
  //// dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  // You can add an arbitrary number of signal or background trees
  //// dataloader->AddSignalTree    ( signalTree,     signalWeight );
  //// dataloader->AddBackgroundTree( background, backgroundWeight );
  dataloader->AddSignalTree    ( signal,     signalWeight );
  dataloader->AddBackgroundTree( background, backgroundWeight );

  // To give different trees for training and testing, do as follows:
  //
  //     dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
  //     dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

  // Use the following code instead of the above two or four lines to add signal and background
  // training and test events "by hand"
  // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
  //      variable definition, but simply compute the expression before adding the event
  // ```cpp
  // // --- begin ----------------------------------------------------------
  // std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
  // Float_t  treevars[4], weight;
  //
  // // Signal
  // for (UInt_t ivar=0; ivar<4; ivar++) signalTree->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
  // for (UInt_t i=0; i<signalTree->GetEntries(); i++) {
  //    signalTree->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is testing
  //    // note that the weight can also be event-wise
  //    if (i < signalTree->GetEntries()/2.0) dataloader->AddSignalTrainingEvent( vars, signalWeight );
  //    else                              dataloader->AddSignalTestEvent    ( vars, signalWeight );
  // }
  //
  // // Background (has event weights)
  // background->SetBranchAddress( "weight", &weight );
  // for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
  // for (UInt_t i=0; i<background->GetEntries(); i++) {
  //    background->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is testing
  //    // note that the weight can also be event-wise
  //    if (i < background->GetEntries()/2) dataloader->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
  //    else                                dataloader->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
  // }
  // // --- end ------------------------------------------------------------
  // ```
  // End of tree registration

  // Set individual event weights (the variables must exist in the original TTree)
  // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
  // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
  // dataloader->SetBackgroundWeightExpression( "weight" );
  // dataloader->SetSignalWeightExpression("pthatweight*Ncoll");

  // Apply additional cuts on the signal and background samples (can be different)
  //// TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  //// TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  TString cuts = Form("(%s) && Dpt>%f && Dpt<%f", mycuts.c_str(), ptmin, ptmax);
  TString cutb = Form("(%s) && Dpt>%f && Dpt<%f", mycutb.c_str(), ptmin, ptmax);

  TCut mycutS = (TCut)cuts;
  TCut mycutB = (TCut)cutb;

  // Tell the dataloader how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used
  // for training, and the other half for testing:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  //
  // To also specify the number of testing events, use:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut,
  //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
  //// dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
  ////                                      "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
  dataloader->PrepareTrainingAndTestTree( mycutS, mycutB,
                                          // "nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=10000:nTest_Background=10000:SplitMode=Random:NormMode=NumEvents:!V" );
                                          // "nTrain_Signal=0:nTrain_Background=100000:nTest_Signal=0:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V" );
                                          "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  // ### Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                         "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
  // "!H:!V:NTrees=%d:MinNodeSize=15.5:BoostType=Grad:Shrinkage=0.275:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                         "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
//   "!H:!V:NTrees=%d:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  if (Use["BDTB"]) // Bagging
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
                         "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
//  "!H:!V:NTrees=%d:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                         "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
  // "!H:!V:NTrees=%d:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );


  if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
                         "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );
// "!H:!V:NTrees=$d:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );


  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                         "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

  // For an example of the category classifier usage, see: TMVAClassificationCategory
  //
  // --------------------------------------------------------------------------------------------------
  //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
  // STILL EXPERIMENTAL and only implemented for BDT's !
  //
  //     factory->OptimizeAllMethods("SigEffAt001","Scan");
  //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
  //
  // --------------------------------------------------------------------------------------------------

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  outf->cd("dataset");
  TTree* info = new TTree("tmvainfo", "TMVA info");
  info->Branch("cuts", &cuts);
  info->Branch("cutb", &cutb);
  info->Branch("var", &varinfo);
  info->Fill();
  info->Write();

  // Save the output
  outf->Close();

  std::cout << "==> Wrote root file: " << outf->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
  delete dataloader;
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfname.c_str() );

  // gSystem->Exec(Form("mv dataset/weights/*.* dataset/weights/%s/", outputstr.c_str()));

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==8)
    { 
      for(int i=0; i<mytmva::nptbins; i++)
        {
          TMVAClassification(argv[1], argv[2], argv[3], argv[4], argv[5], mytmva::ptbins[i], mytmva::ptbins[i+1], argv[6], argv[7]); 
        }
      return 0;
    }

  return 1;
}
