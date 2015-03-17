#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>

void prettyPlots()
{
  const int NTRIG = 5;
  TGraphAsymmErrors *turnOns[NTRIG];

  // TString trigname[NTRIG] = {"HLT_HISinglePhoton10_v1",
  // 			     "HLT_HISinglePhoton15_v1",
  // 			     "HLT_HISinglePhoton20_v1",
  // 			     "HLT_HISinglePhoton40_v1",
  // 			     "HLT_HISinglePhoton60_v1"};

  TString trigname[NTRIG] = {"HLT_HISinglePhoton10_barrel_v1",
  			     "HLT_HISinglePhoton15_barrel_v1",
  			     "HLT_HISinglePhoton20_barrel_v1",
  			     "HLT_HISinglePhoton40_barrel_v1",
  			     "HLT_HISinglePhoton60_barrel_v1"};

  Int_t trigColors[NTRIG] = {1, kBlue, kRed, 90, kMagenta};

  TFile *inFile = TFile::Open("photonTurnOn_barrel.root");
  for(int i = 0; i < NTRIG; i++)
  {
    turnOns[i] = (TGraphAsymmErrors *)inFile->Get(trigname[i]+"_asymm");
    turnOns[i]->SetMarkerColor(trigColors[i]);
    turnOns[i]->SetLineColor(trigColors[i]);
  }

  TH1D *hEmpty = new TH1D("hEmpty",";p_{T}^{#gamma};Efficiency",20,0,100);
  hEmpty->SetMaximum(1.2);
  hEmpty->SetMinimum(0.0);

  TCanvas *c1 = new TCanvas();//("c1","c1",700,700);
  hEmpty->Draw();
  TLine *oneLine = new TLine(0,1,100,1);
  oneLine->SetLineStyle(2);
  oneLine->Draw();

  TLegend *leg = new TLegend(0.6, 0.2, 0.9, 0.6);
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  leg->SetTextSize(14);

  for(int i = 0; i < NTRIG; ++i)
  {
    turnOns[i]->Draw("p e");
    leg->AddEntry(turnOns[i], trigname[i], "p l");
  }
  leg->Draw();

  c1->SaveAs("photon_turnOn_barrel.pdf");

}
