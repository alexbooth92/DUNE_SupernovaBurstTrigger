#include <iostream>
#include <vector>
#include <fstream>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>
#include <TROOT.h>
#include "/Users/alexanderbooth/Documents/Work/UsefulCode/colors.h"


int main()
{
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_"+s_Filename+".txt");

  std::map<int,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int Config;
  double Eff, Bkgd;
  while(inFile >> Config >> Eff >> Bkgd)
  {
    map_ConfigToEffAndBkgd[Config] = {Eff,Bkgd};
  }

  gROOT->ForceStyle();
  std::vector<int> vec_Colors = getColors(2);
  TFile *f_SNTheoryDistributions = new TFile("SNTheoryDistributions.root","READ");
  TH1D  *h_SNProbabilityVDistance = (TH1D*)f_SNTheoryDistributions->Get("h_SNProbabilityVDistance");
  h_SNProbabilityVDistance->SetLineWidth(3);
  h_SNProbabilityVDistance->SetLineColor(46);

  TFile *f_Input = new TFile("Analyse_SNBurst_GH_SNMC.root", "READ");

  std::vector<TH1D*>   vec_h_FakeRateVNClusters;
  std::vector<TH1D*>   vec_h_EfficiencyVEvents;
  std::vector<TH1D*>   vec_h_EfficiencyVDistance;
  std::vector<TH1D*>   vec_h_EffGalaxy;
  std::vector<TGraph*> vec_ROC;
  for(unsigned int config = 0; config < map_ConfigToEffAndBkgd.size(); config++)
  {
    TString s_FakeRateVNClusters  = Form("h_FakeRateVNClusters_Config%i",    config);
    TString s_EfficiencyVEvents   = Form("h_EfficiencyVEvents_Config%i",     config);
    TString s_EfficiencyVDistance = Form("h_EfficiencyVDistance_Config%i",   config);
    TString s_EffGalaxy           = Form("h_NeighbourhoodEffiency_Config%i", config);
    TString s_ROC                 = Form("g_ROC_Config%i",                   config);
    
    TH1D *h_FakeRateVNClusters = (TH1D*)f_Input->Get(s_FakeRateVNClusters); 

    TH1D *h_EfficiencyVEvents = (TH1D*)f_Input->Get(s_EfficiencyVEvents);
    int   maxBin_Events = h_EfficiencyVEvents->GetMaximumBin();
    for(int i = maxBin_Events; i < h_EfficiencyVEvents->GetSize()-1; i++)
    {
      h_EfficiencyVEvents->SetBinContent(i,1);
    }
    TH1D *h_EfficiencyVDistance = (TH1D*)f_Input->Get(s_EfficiencyVDistance);
    int   maxBin_Distance = h_EfficiencyVDistance->GetMaximumBin();
    for(int i = maxBin_Distance; i > 0; i--)
    {
      h_EfficiencyVDistance->SetBinContent(i,1);
    }

    TH1D   *h_EffGalaxy           = (TH1D*)f_Input->Get(s_EffGalaxy);
    TGraph *g_ROC                 = (TGraph*)f_Input->Get(s_ROC);

    h_FakeRateVNClusters->SetLineColor(vec_Colors.at(config));
    h_EfficiencyVEvents->SetLineColor(vec_Colors.at(config));
    h_EfficiencyVDistance->SetLineColor(vec_Colors.at(config));
    h_EffGalaxy->SetLineColor(vec_Colors.at(config));
    g_ROC->SetMarkerColor(vec_Colors.at(config));
    g_ROC->SetMarkerStyle(3);

    vec_h_FakeRateVNClusters.push_back(h_FakeRateVNClusters);
    vec_h_EfficiencyVEvents.push_back(h_EfficiencyVEvents);
    vec_h_EfficiencyVDistance.push_back(h_EfficiencyVDistance);
    vec_h_EffGalaxy.push_back(h_EffGalaxy);
    vec_ROC.push_back(g_ROC);
  }

  TCanvas *c_FakeRateVNClusters   = new TCanvas("c_FakeRateVNClusters", "c_FakeRateVNClusters", 800, 500);
  THStack *stk_FakeRateVNClusters = new THStack("stk_FakeRateVNClusters", "Number of Clusters in Time Window Required to Trigger vs. Trigger Rate");
  TLegend *leg_FakeRateVNClusters = new TLegend(0.5, 0.105, 0.8, 0.3);
  leg_FakeRateVNClusters->SetTextSize(0.023);

  leg_FakeRateVNClusters->SetHeader("Individual Marley Efficiency & 10kt Background Rate");
  for(unsigned int config = 0; config < vec_h_FakeRateVNClusters.size(); config++)
  {
    stk_FakeRateVNClusters->Add(vec_h_FakeRateVNClusters.at(config));
    leg_FakeRateVNClusters->AddEntry(vec_h_FakeRateVNClusters.at(config), 
                                    Form("Eff: %.2f, Bkgd: %.2fHz", map_ConfigToEffAndBkgd[config].first, map_ConfigToEffAndBkgd[config].second), "L");
  }
  c_FakeRateVNClusters->Draw();
  c_FakeRateVNClusters->SetLogy();
  stk_FakeRateVNClusters->Draw("NOSTACK");
  stk_FakeRateVNClusters->GetXaxis()->SetTitle("Number of Clusters/Time Window");
  stk_FakeRateVNClusters->GetYaxis()->SetTitle("Trigger Rate, (Hz)");
  leg_FakeRateVNClusters->Draw();
  TText *t_perMonth  = new TText(vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax()-15, 2e-7, "1/Month");
  TText *t_perWeek   = new TText(vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax()-15, 1.65e-6, "1/Week");
  TText *t_perDay    = new TText(vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax()-15, 1.16e-5, "1/Day");
  TLine *l_perMonth  = new TLine(vec_h_FakeRateVNClusters.at(vec_h_FakeRateVNClusters.size()-1)->GetXaxis()->GetXmin(), 4.13e-7, vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax(), 4.13e-7);
  TLine *l_perWeek   = new TLine(vec_h_FakeRateVNClusters.at(vec_h_FakeRateVNClusters.size()-1)->GetXaxis()->GetXmin(), 1.65e-6, vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax(), 1.65e-6);
  TLine *l_perDay    = new TLine(vec_h_FakeRateVNClusters.at(vec_h_FakeRateVNClusters.size()-1)->GetXaxis()->GetXmin(), 1.16e-5, vec_h_FakeRateVNClusters.at(0)->GetXaxis()->GetXmax(), 1.16e-5);
  l_perMonth->SetLineColor(4);
  l_perWeek ->SetLineColor(4);
  l_perDay  ->SetLineColor(4);
  l_perMonth->SetLineWidth(3);
  l_perWeek ->SetLineWidth(3);
  l_perDay  ->SetLineWidth(3);
  t_perMonth->Draw();
  t_perWeek->Draw();
  t_perDay->Draw();
  l_perMonth->Draw();
  l_perWeek->Draw();
  l_perDay->Draw();
  c_FakeRateVNClusters->SaveAs("FakeRateVNClusters.pdf");


  TCanvas *c_EfficiencyVEvents   = new TCanvas("c_EfficiencyVEvents", "c_EfficiencyVEvents", 800, 500);
  THStack *stk_EfficiencyVEvents = new THStack("stk_EfficiencyVEvents", "Efficiency vs. Number of Events in SN Burst, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVEvents = new TLegend(0.53, 0.65, 0.85, 0.85);
  leg_EfficiencyVEvents->SetTextSize(0.023);

  leg_EfficiencyVEvents->SetHeader("Individual Marley Efficiency & 10kt Background Rate");
  for(unsigned int config = 0; config < vec_h_EfficiencyVEvents.size(); config++)
  {
    stk_EfficiencyVEvents->Add(vec_h_EfficiencyVEvents.at(config));
    leg_EfficiencyVEvents->AddEntry(vec_h_EfficiencyVEvents.at(config), 
                                    Form("Eff: %.2f, Bkgd: %.2fHz", map_ConfigToEffAndBkgd[config].first, map_ConfigToEffAndBkgd[config].second), "L");
  }
  c_EfficiencyVEvents->Draw();
  c_EfficiencyVEvents->SetLogx();
  stk_EfficiencyVEvents->Draw("NOSTACK");
  stk_EfficiencyVEvents->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_EfficiencyVEvents->GetYaxis()->SetTitle("Efficiency");
  leg_EfficiencyVEvents->Draw();
  c_EfficiencyVEvents->SaveAs("EfficiencyVEvents.pdf");


  TCanvas *c_EfficiencyVDistance   = new TCanvas("c_EfficiencyVDistance", "c_EfficiencyVDistance", 800, 500);
  THStack *stk_EfficiencyVDistance = new THStack("stk_EfficiencyVDistance", "Efficiency vs. Distance to SN, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVDistance = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EfficiencyVDistance->SetTextSize(0.023);

  leg_EfficiencyVDistance->SetHeader("Individual Marley Efficiency & 10kt Background Rate");
  for(unsigned int config = 0; config < vec_h_EfficiencyVDistance.size(); config++)
  {
    stk_EfficiencyVDistance->Add(vec_h_EfficiencyVDistance.at(config));
    leg_EfficiencyVDistance->AddEntry(vec_h_EfficiencyVDistance.at(config), 
                                      Form("Eff: %.2f, Bkgd: %.2fHz", map_ConfigToEffAndBkgd[config].first, map_ConfigToEffAndBkgd[config].second), "L");
  }
  c_EfficiencyVDistance->Draw();
  c_EfficiencyVDistance->SetLogx();
  stk_EfficiencyVDistance->Draw("NOSTACK");
  stk_EfficiencyVDistance->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EfficiencyVDistance->GetYaxis()->SetTitle("Efficiency");
  leg_EfficiencyVDistance->Draw();
  c_EfficiencyVDistance->SaveAs("EfficiencyVDistance.pdf");


  TCanvas *c_EffGalaxy   = new TCanvas("c_EffGalaxy", "c_EffGalaxy", 800, 500);
  gStyle->SetOptStat(0);
  THStack *stk_EffGalaxy = new THStack("stk_EffGalaxy", "Galactic Neighbourhood Coverage, Fake Trigger Rate 1/Month");
  TLegend *leg_EffGalaxy = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EffGalaxy->SetTextSize(0.023);

  leg_EffGalaxy->SetHeader("Individual Marley Efficiency & 10kt Background Rate");
  stk_EffGalaxy->Add(h_SNProbabilityVDistance);
  for(unsigned int config = 0; config < vec_h_EffGalaxy.size(); config++)
  {
    stk_EffGalaxy->Add(vec_h_EffGalaxy.at(config));
    leg_EffGalaxy->AddEntry(vec_h_EffGalaxy.at(config), Form("Eff: %.2f, Bkgd: %.2fHz", 
                            map_ConfigToEffAndBkgd[config].first, map_ConfigToEffAndBkgd[config].second), "L");
  }
  leg_EffGalaxy->AddEntry(h_SNProbabilityVDistance, "SN Probability", "L");
  c_EffGalaxy->Draw();
  c_EffGalaxy->SetLogy();
  stk_EffGalaxy->Draw("NOSTACK");
  stk_EffGalaxy->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EffGalaxy->GetYaxis()->SetTitle("Efficiency x SN Probability");
  stk_EffGalaxy->Draw("NOSTACK");
  leg_EffGalaxy->Draw();
  c_EffGalaxy->SaveAs("EffGalaxy.pdf");


  TCanvas *c_ROC      = new TCanvas("c_ROC", "c_ROC", 800, 500);
  TLegend *leg_ROC    = new TLegend(0.15, 0.68, 0.48, 0.88);
  leg_ROC->SetTextSize(0.023);

  leg_ROC->SetHeader("Individual Marley Efficiency & 10kt Background Rate");
  leg_ROC->AddEntry(vec_ROC.at(0), Form("Eff: %.2f, Bkgd: %.2fHz", map_ConfigToEffAndBkgd[0].first, map_ConfigToEffAndBkgd[0].second), "p");
  c_ROC->Draw();
  c_ROC->SetLogy();
  vec_ROC.at(0)->Draw("AP");
  for(unsigned int config = 1; config < vec_ROC.size(); config++)
  {
    leg_ROC->AddEntry(vec_ROC.at(config), Form("Eff: %.2f, Bkgd: %.2fHz", map_ConfigToEffAndBkgd[config].first, map_ConfigToEffAndBkgd[config].second), "p");
    vec_ROC.at(config)->Draw("P");
  }
  vec_ROC.at(0)->GetXaxis()->SetRangeUser(0.8, 1);
  vec_ROC.at(0)->GetYaxis()->SetRangeUser(10e-15, 20e2);
  vec_ROC.at(0)->SetTitle("Fake Trigger Rate vs. Galactic Neighbourhood Coverage");
  vec_ROC.at(0)->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  vec_ROC.at(0)->GetYaxis()->SetTitle("Fake Trigger Rate, (Hz)");
  TLine *l_perMonth_2 = new TLine(0.8, 4.13e-7, 1, 4.13e-7);
  l_perMonth_2->SetLineColor(1);
  l_perMonth_2->SetLineWidth(3);
  TText *t_perMonth_2 = new TText(0.815, 8e-7, "1/Month");

  leg_ROC->Draw();
  l_perMonth_2->Draw();
  t_perMonth_2->Draw();
  c_ROC->SaveAs("ROC.pdf");

  return 0;
}
