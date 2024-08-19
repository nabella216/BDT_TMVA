// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/efficiencies_8cxx_source.html
#include <iostream>
#include <string>

#include <TKey.h>
#include <TList.h>
#include "TMVA/efficiencies.h"
#include "TMVA/mvas.h"
#include "TMVA/correlations.h"

#include "mvaeffs.h"
#include "TMVAClassification.h"
#include "xjjcuti.h"

namespace mytmva
{
  void guiefficiencies(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage = "0,1,2,3");
  void efficiencies(std::string outfname);
}

void mytmva::guiefficiencies(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage/* = "0,1,2,3"*/)
{
  std::string outfname = mytmva::mkname(outputname, ptmin, ptmax, mymethod, stage);
  mytmva::efficiencies(outfname);
}

void mytmva::efficiencies(std::string outfname)
{
  TString dataset("");
  std::string outputstr = xjjc::str_replaceallspecial(outfname);

  // set up dataset
  TFile* file = TFile::Open( outfname.c_str() );
  if(!file)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify filename." << std::endl; return; }
  if(file->GetListOfKeys()->GetEntries()<=0)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify if dataset exist." << std::endl; return; }
  // --->>
  if( (dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()==1))
    { dataset = ((TKey*)file->GetListOfKeys()->At(0))->GetName(); }
  // <<---
  else if((dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()>=1))
    {
      std::cout << "==> Warning " << __FUNCTION__ << ", more than 1 dataset." << std::endl;
      file->ls(); return;
    }
  else { return; }

  TMVA::efficiencies(dataset.Data(), outfname.c_str(), 2); // 1,2,3
  TMVA::mvas(dataset.Data(), outfname.c_str(), TMVA::kCompareType);
  TMVA::correlations(dataset.Data(), outfname.c_str());
  mytmva::mvaeffs(dataset.Data(), outfname.c_str());
  mytmva::mvaeffs(dataset.Data(), outfname.c_str(), 4.e+2, 1.e+5);

  // gSystem->Exec(Form("rm %s/plots/*.png", dataset.Data()));
  gSystem->Exec(Form("mkdir -p %s/plots/%s", dataset.Data(), outputstr.c_str()));
  gSystem->Exec(Form("mv %s/plots/*.png %s/plots/%s/", dataset.Data(), dataset.Data(), outputstr.c_str()));
}

int main(int argc, char* argv[])
{
  if(argc==4)
    { 
      for(int i=0; i<mytmva::nptbins; i++)
        {
          mytmva::guiefficiencies(argv[1], mytmva::ptbins[i], mytmva::ptbins[i+1], argv[2], argv[3]); 
        }
      return 0; 
    }
  return 1;
}
