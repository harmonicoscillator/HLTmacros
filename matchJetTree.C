#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <TMath.h>

#include "EventMatchingCMS.h"

const int MAXJETS = 500;
const int nBins = 100;
const double maxPt = 300;
const float offlineEtaCut = 2.0;
const bool offlineIso = false;

void matchJetTree(TString HiForestName, TString bitfileName, TString outFileName)
{
  TFile *HLTFile = TFile::Open(bitfileName);
  TTree *HLTTree = (TTree*)HLTFile->Get("hltbitanalysis/HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  const int NTRIG = 18;
  TString trigname[NTRIG] = {"L1_SingleS1Jet4_BptxAND",
			     "L1_SingleS1Jet8_BptxAND",
			     "L1_SingleS1Jet16_BptxAND",
			     "L1_SingleS1Jet28_BptxAND",
			     "L1_SingleJet36_BptxAND",
			     "L1_SingleS1Jet40_BptxAND",
			     "L1_SingleJet44_BptxAND",
			     "L1_SingleS1Jet56_BptxAND",
			     "L1_SingleJet68_BptxAND",
			     "L1_SingleJet80_BptxAND",
			     "L1_SingleJet92_BptxAND",
			     "L1_SingleJet128_BptxAND",
			     "HLT_PuAK4CaloJet40_v1",
			     "HLT_PuAK4CaloJet60_v1",
			     "HLT_PuAK4CaloJet80_v1",
			     "HLT_PuAK4CaloJet100_v1",
			     "HLT_PuAK4CaloJet110_v1",
			     "HLT_PuAK4CaloJet120_v1"};

  Int_t triggers[NTRIG];

  HLTTree->SetBranchAddress("Event",&hlt_event);
  HLTTree->SetBranchAddress("Run",&hlt_run);
  HLTTree->SetBranchAddress("LumiBlock",&hlt_lumi);

  for(int i = 0; i < NTRIG; i++)
  {
    HLTTree->SetBranchAddress(trigname[i], &(triggers[i]));
  }

  TFile *inFile = TFile::Open(HiForestName);
  TTree *f1Tree = (TTree*)inFile->Get("akPu4CaloJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("vz",&vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t f_num;
  Float_t f_pt[MAXJETS];
  Float_t f_eta[MAXJETS];
  Float_t f_phi[MAXJETS];
  Float_t f_rawpt[MAXJETS];

  f1Tree->SetBranchAddress("nref",&f_num);
  f1Tree->SetBranchAddress("jtpt",f_pt);
  f1Tree->SetBranchAddress("jteta",f_eta);
  f1Tree->SetBranchAddress("jtphi",f_phi);
  f1Tree->SetBranchAddress("rawpt",f_rawpt);

  TH1D *fPt = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  TH1D *accepted[NTRIG];

  for(int i = 0; i < NTRIG; ++i)
  {
    accepted[i] = new TH1D(Form("accepted_pt%d",i),";offline p_{T}",nBins,0,maxPt);
  }

  // Make the event-matching map ************
  EventMatchingCMS *matcher = new EventMatchingCMS();
  int duplicates = 0;
  std::cout << "Begin making map." << std::endl;
  Long64_t l_entries = HLTTree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    HLTTree->GetEntry(j);
    bool status = matcher->addEvent(hlt_event, hlt_lumi, hlt_run, j);
    if(status == false)
      duplicates++;
  }
  std::cout << "Finished making map." << std::endl;
  std::cout << "Duplicate entries: " << duplicates << std::endl;
  // **********************

  // analysis loop
  int matched = 0;
  for(Long64_t entry = 0; entry < fEvtTree->GetEntries(); ++entry)
  {
    fEvtTree->GetEntry(entry);
    long long hlt_entry = matcher->retrieveEvent(f_evt, f_lumi, f_run);
    if(hlt_entry == -1) continue;
    matched++;

    fSkimTree->GetEntry(entry);

    bool goodEvent = false;
    if((pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15))
    {
      goodEvent = true;
    }
    if(!goodEvent) continue;

    f1Tree->GetEntry(entry);
    HLTTree->GetEntry(hlt_entry);

    double maxfpt = 0;
    for(int i = 0; i < f_num; ++i)
    {
      if(TMath::Abs(f_eta[i]) > offlineEtaCut) continue;
      if(f_pt[i] > maxfpt) {
	maxfpt = f_pt[i];
      }
    }

    fPt->Fill(maxfpt);
    for(int i = 0; i < NTRIG; ++i)
    {
      if(triggers[i])
      {
	accepted[i]->Fill(maxfpt);
      }
    }
  }
  std::cout << "Events matched: " << matched << std::endl;

  //make turn-on curves
  TGraphAsymmErrors *a_pt[NTRIG];
  for(int i = 0; i < NTRIG; ++i){
    a_pt[i] = new TGraphAsymmErrors();
    a_pt[i]->BayesDivide(accepted[i],fPt);
    a_pt[i]->SetName(trigname[i]+"_asymm");
  }

  //save output
  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  outFile->cd();
  fPt->Write();
  for(int i = 0; i < NTRIG; ++i)
  {
    accepted[i]->Write();
    a_pt[i]->Write();
  }

  HLTFile->Close();
  inFile->Close();
  outFile->Close();
}

int main(int argc, char **argv)
{
  if(argc == 4)
  {
    matchJetTree(argv[1], argv[2], argv[3]);
    return 0;
  }
  return 1;
}
