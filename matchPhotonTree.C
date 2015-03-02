#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>

const TString AnaFilename = "/export/d00/scratch/luck/hlt_Ian_lowlumi_V2_neutrinogun_photonanalyzer.root";
// eta = 3.0 triggers
//const TString HLTFilename = "/export/d00/scratch/luck/hlt_jetsphotons_lowlumi_V6_neutrinogun_openHLT.root";
// eta = 1.479 triggers
const TString HLTFilename = "/export/d00/scratch/luck/hlt_Ian_lowlumi_V2_neutrinogun_openHLT.root";

const int nBins = 100;
const double maxpt = 100;

void matchPhotonTree()
{
  TFile *HLTFile = TFile::Open(HLTFilename);
  TTree *HLTTree = (TTree*)HLTFile->Get("HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;
  Int_t HLT_HISinglePhoton10_v1;
  Int_t HLT_HISinglePhoton15_v1;
  Int_t HLT_HISinglePhoton20_v1;
  Int_t HLT_HISinglePhoton40_v1;
  Int_t HLT_HISinglePhoton60_v1;

  TString trigname[5] = {"HLT_HISinglePhoton10_v1",
			 "HLT_HISinglePhoton15_v1",
			 "HLT_HISinglePhoton20_v1",
			 "HLT_HISinglePhoton40_v1",
			 "HLT_HISinglePhoton60_v1"};

  HLTTree->SetBranchAddress("Event",&hlt_event);
  HLTTree->SetBranchAddress("Run",&hlt_run);
  HLTTree->SetBranchAddress("LumiBlock",&hlt_lumi);

  HLTTree->SetBranchAddress(trigname[0],&HLT_HISinglePhoton10_v1);
  HLTTree->SetBranchAddress(trigname[1],&HLT_HISinglePhoton15_v1);
  HLTTree->SetBranchAddress(trigname[2],&HLT_HISinglePhoton20_v1);
  HLTTree->SetBranchAddress(trigname[3],&HLT_HISinglePhoton40_v1);
  HLTTree->SetBranchAddress(trigname[4],&HLT_HISinglePhoton60_v1);

  TFile *AnaFile = TFile::Open(AnaFilename);
  TTree *AnaTree = (TTree*)AnaFile->Get("SimpleGedPhotonAnalyzer/PhotonTree"); // other option is SimplePhotonAnalyzer

  Int_t ana_event, ana_run, ana_lumi;
  Int_t nPhotons;
  Double_t pt[500], eta[500], phi[500];

  AnaTree->SetBranchAddress("event",&ana_event);
  AnaTree->SetBranchAddress("run",&ana_run);
  AnaTree->SetBranchAddress("lumi",&ana_lumi);

  AnaTree->SetBranchAddress("nPhotons",&nPhotons);
  AnaTree->SetBranchAddress("pt",pt);
  AnaTree->SetBranchAddress("eta",eta);
  AnaTree->SetBranchAddress("phi",phi);

  //book histos
  TH1D *hists_pt[6], *hists_eta[6];
  hists_pt[0] = new TH1D("leading_photon_pt",";p_{T}^{#gamma}",nBins,0,maxpt);
  hists_eta[0] = new TH1D("leading_photon_eta",";#eta^{#gamma}",nBins,-5,5);
  for(int i = 0; i < 5; ++i)
  {
    hists_pt[i+1] = (TH1D*)hists_pt[0]->Clone(trigname[i]);
    hists_eta[i+1] = (TH1D*)hists_eta[0]->Clone(trigname[i]+"eta");
  }

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << AnaTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();
  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry)
  {
    HLTTree->GetEntry(entry);
    matcher->addEvent(hlt_event, hlt_lumi, hlt_run, entry);
  }

  // analysis loop
  int matched = 0;
  for(Long64_t entry = 0; entry < AnaTree->GetEntries(); ++entry)
  {
    AnaTree->GetEntry(entry);
    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, ana_run);
    if(hlt_entry == -1) continue;
    HLTTree->GetEntry(hlt_entry);
    matched++;

    Double_t maxAnaPt = -1;
    Double_t maxAnaEta = -100;
    for(int i = 0; i < nPhotons; ++i)
    {
      if(fabs(eta[i]) > 1.44) continue;
      if(pt[i] > maxAnaPt)
      {
	maxAnaPt = pt[i];
	maxAnaEta = eta[i];
      }
    }

    if(maxAnaPt == -1) continue;

    hists_pt[0]->Fill(maxAnaPt);
    hists_eta[0]->Fill(maxAnaEta);
    if(HLT_HISinglePhoton10_v1){
      hists_pt[1]->Fill(maxAnaPt);
      hists_eta[1]->Fill(maxAnaEta);
    }
    if(HLT_HISinglePhoton15_v1){
      hists_pt[2]->Fill(maxAnaPt);
      hists_eta[2]->Fill(maxAnaEta);
    }
    if(HLT_HISinglePhoton20_v1){
      hists_pt[3]->Fill(maxAnaPt);
      hists_eta[3]->Fill(maxAnaEta);
    }
    if(HLT_HISinglePhoton40_v1){
      hists_pt[4]->Fill(maxAnaPt);
      hists_eta[4]->Fill(maxAnaEta);
    }
    if(HLT_HISinglePhoton60_v1){
      hists_pt[5]->Fill(maxAnaPt);
      hists_eta[5]->Fill(maxAnaEta);
    }
  }
  std::cout << "Events matched: " << matched << std::endl;

  //make turn-on curves
  TGraphAsymmErrors *a_pt[5], *a_eta[5];
  for(int i = 0; i < 5; ++i){
    a_pt[i] = new TGraphAsymmErrors();
    a_pt[i]->BayesDivide(hists_pt[i+1],hists_pt[0]);
    a_pt[i]->SetName(trigname[i]+"_asymm");
    a_eta[i] = new TGraphAsymmErrors();
    a_eta[i]->BayesDivide(hists_eta[i+1],hists_eta[0]);
    a_eta[i]->SetName(trigname[i]+"_eta_asymm");
  }

  //save output
  TFile *outFile = TFile::Open("photonTurnOn_barrel.root","RECREATE");
  outFile->cd();
  hists_pt[0]->Write();
  hists_eta[0]->Write();
  for(int i = 0; i < 5; ++i)
  {
    hists_pt[i+1]->Write();
    hists_eta[i+1]->Write();
    a_pt[i]->Write();
    a_eta[i]->Write();
  }

  HLTFile->Close();
  AnaFile->Close();
  outFile->Close();
}

int main()
{
  matchPhotonTree();
  return 0;
}
