#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <TText.h>
#include <TFrame.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include "Module_SNClustering_Config.h"
#include "/Users/alexanderbooth/Documents/Work/UsefulCode/colors.h"


std::vector<int> vec_Colors = getColors(2);
std::vector<std::pair<double,double>> vec_EffAndBkgd;

bool makePlot1(TFile *f_Input_Module, TFile *f_Input_Analyse, TFile *f_Input_Reco, TCanvas *c)
{
  std::vector<TGraph*> vec_Module;
  double vupp(0);
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TGraph *g_Module = (TGraph*)f_Input_Module->Get(Form("g_ConfigDefinitions%i",i));
    g_Module->SetMarkerStyle(3);
    g_Module->SetMarkerColor(vec_Colors.at(i));
    g_Module->SetLineColor(vec_Colors.at(i));
    for(int j = 0; j < g_Module->GetN(); j++)
    {
      double x, y;
      g_Module->GetPoint(j,x,y);
      if(y > vupp)
      {
        vupp = y;
      }
    }
    vec_Module.push_back(g_Module);
  }
  TCanvas *c_Module = new TCanvas("c_Module", "c_Module", 800, 500);
  c_Module->SetLogy();
  vec_Module.at(0)->Draw("APL");
  for(unsigned int i = 1; i < vec_Module.size(); i++)
  {
    vec_Module.at(i)->Draw("PL");
  }
  double yAxis = c_Module->GetFrame()->GetY1();
  TText *t_1 = new TText(1,yAxis,"Adjacent Channel Tolerance");
  t_1->SetTextSize(0.023);
  t_1->SetTextAngle(10);
  TText *t_2 = new TText(2,yAxis,"Adjacent Time Hits");
  t_2->SetTextSize(0.023);
  t_2->SetTextAngle(10);
  TText *t_3 = new TText(3,yAxis,"Minimum Channels");
  t_3->SetTextSize(0.023);
  t_3->SetTextAngle(10);
  TText *t_4 = new TText(4,yAxis,"Cluster Width");
  t_4->SetTextSize(0.023);
  t_4->SetTextAngle(10);
  TText *t_5 = new TText(5,yAxis,"Adjacent Hit Time Separation");
  t_5->SetTextSize(0.023);
  t_5->SetTextAngle(10);
  TText *t_6 = new TText(6,yAxis,"Total ADC");
  t_6->SetTextSize(0.023);
  t_6->SetTextAngle(10);
  vec_Module.at(0)->GetXaxis()->SetLimits(0.5,NCuts+0.5);
  vec_Module.at(0)->GetYaxis()->SetRangeUser(0.9,vupp+vupp*0.05);
  vec_Module.at(0)->GetXaxis()->SetLabelSize(0);
  vec_Module.at(0)->GetXaxis()->SetNdivisions(NCuts);
  vec_Module.at(0)->SetTitle("Parameters Defining Cluster Configurations");
  vec_Module.at(0)->GetYaxis()->SetTitle("Value");
  t_1->Draw();
  t_2->Draw();
  t_3->Draw();
  t_4->Draw();
  t_5->Draw();
  t_6->Draw();
  c_Module->SaveAs("ConfigurationDefinitions.pdf");

  c->cd();
  c->SetLogy();
  TGraph *g_Analyse = (TGraph*)f_Input_Analyse->Get("g_BkgdVsEff");
  TGraph *g_Module0 = new TGraph(1);
  double x0, y0;
  g_Analyse->GetPoint(0,x0,y0);
  g_Module0->SetPoint(0,x0,y0);
  vec_EffAndBkgd.push_back({x0,y0});
  g_Module0->SetMarkerStyle(3);
  g_Module0->SetMarkerColor(vec_Colors.at(0));
  g_Module0->SetMarkerSize(2);
  g_Module0->Draw("AP");
  double ylow(y0), yupp(y0);
  for(unsigned int i = 1; i < NConfigs; i++)
  {
    TGraph *g_Module = new TGraph(1);
    double x, y;
    g_Analyse->GetPoint(i,x,y);
    g_Module->SetPoint(0,x,y);
    vec_EffAndBkgd.push_back({x,y});
    g_Module->SetMarkerStyle(3);
    g_Module->SetMarkerColor(vec_Colors.at(i));
    g_Module->SetMarkerSize(2);
    g_Module->Draw("P");
    if(y > yupp)
    {
      yupp = y;
    }
    if(y < ylow)
    {
      ylow = y;
    }
  }
  g_Module0->GetXaxis()->SetLimits(0.5,1);
  g_Module0->GetYaxis()->SetRangeUser(ylow, yupp);
  g_Module0->SetTitle("Background Rate vs. Individual SN Neutrino Triggering Efficiency 10kt");
  g_Module0->GetXaxis()->SetTitle("Efficiency");
  g_Module0->GetYaxis()->SetTitle("Background Rate, (Hz)");
  g_Module0->SetMarkerStyle(3);
  c->SaveAs("BkgdVsIndividualEff.pdf");
  c->Clear();

  return true;
}


bool makePlot2(TFile *f_Input_Module, TFile *f_Input_Analyse, TCanvas *c)
{
  c->SetLogy(0);
  std::vector<TH1D*> vec_h_ENu;
  TH1D *h_ENu_MC = (TH1D*)f_Input_Module->Get("h_ENu_MC");
  THStack *stk_ENu = new THStack("stk_ENu","Individual SN Neutrino Triggering Efficiency 10kt vs. True Neutrino Energy");
  TLegend *leg_ENu = new TLegend(0.65,0.15,0.895,0.5);
  leg_ENu->SetTextSize(0.023);
  leg_ENu->SetHeader("Overall Efficiency & 10kt Background Rate");
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TH1D *h_ENu = (TH1D*)f_Input_Analyse->Get(Form("h_EfficiencyVEnergy%i",i));
    h_ENu->Divide(h_ENu,h_ENu_MC,1.0,1.0,"B");
    h_ENu->SetLineColor(vec_Colors.at(i));
    vec_h_ENu.push_back(h_ENu);
    stk_ENu->Add(h_ENu);
    leg_ENu->AddEntry(h_ENu, Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
  }
  TLine *l_Threshold = new TLine(1.5,0,1.5,1);
  l_Threshold->SetLineWidth(3);
  stk_ENu->Draw("NOSTACK");
  stk_ENu->GetXaxis()->SetTitle("True Neutrino Energy, (MeV)");
  stk_ENu->GetYaxis()->SetTitle("Efficiency");
  leg_ENu->Draw();
  l_Threshold->Draw();
  c->SaveAs("IndividualEffVsEnergy.pdf");
  c->Clear();

  c->SetLogy();
  int startEnergy   = 7;
  int nEnergyCurves = 8;
  int curveCount    = 0;
  double bkgdMin    = 1000000;
  double bkgdMax    = 0;
    
  std::vector<TGraph*> vec_g_ROC_EnergySplit;
  std::map<std::pair<double,int>,std::pair<double,double>> map_Efficiency_EnergySplit;
  std::map<double,int> map_EnergyToCurve;
  
  std::map<double,TH1D*> map_BkgdRateToHist;
  for(unsigned int i = 0; i < vec_EffAndBkgd.size(); i++)
  {
    map_BkgdRateToHist[vec_EffAndBkgd.at(i).second] = vec_h_ENu.at(i);
  }

  int count = 0;
  std::map<double,TH1D*>::iterator it_BkgdRateToHist;
  TLegend *leg_ROC_EnergySplit = new TLegend(0.105,0.65,0.22,0.895);
  leg_ROC_EnergySplit->SetHeader("Energy Bin Centre");
  leg_ROC_EnergySplit->SetTextSize(0.023);
  for(it_BkgdRateToHist=map_BkgdRateToHist.begin(); it_BkgdRateToHist!=map_BkgdRateToHist.end(); it_BkgdRateToHist++)
  {
    int startBin = vec_h_ENu.at(count)->FindBin(startEnergy);
    for(int j = startBin; j < startBin+nEnergyCurves; j++)
    {
      double energy = vec_h_ENu.at(count)->GetBinCenter(j);
      double eff    = vec_h_ENu.at(count)->GetBinContent(j);
      std::pair<double,int>    pair_EnergyConfig = {energy,count};
      std::pair<double,double> pair_EffBkgd      = {vec_EffAndBkgd.at(count).second, eff};
      map_Efficiency_EnergySplit[pair_EnergyConfig] = pair_EffBkgd;
      if(count==0)
      {
        TGraph *g_ROC_EnergySplit = new TGraph(vec_EffAndBkgd.size());
        g_ROC_EnergySplit->SetMarkerStyle(3);
        vec_g_ROC_EnergySplit.push_back(g_ROC_EnergySplit);
        vec_g_ROC_EnergySplit.at(curveCount)->SetMarkerColor(vec_Colors.at(curveCount+11));
        vec_g_ROC_EnergySplit.at(curveCount)->SetLineColor(vec_Colors.at(curveCount+11));
        leg_ROC_EnergySplit->AddEntry(vec_g_ROC_EnergySplit.at(curveCount), Form("%.2fMeV",energy), "P");
        map_EnergyToCurve[energy] = curveCount;
        curveCount++;
      }
      if(vec_EffAndBkgd.at(count).second<bkgdMin)
      {
        bkgdMin = vec_EffAndBkgd.at(count).second;
      }
      if(vec_EffAndBkgd.at(count).second>bkgdMax)
      {
        bkgdMax = vec_EffAndBkgd.at(count).second;
      }
    }
    count++;
  }

  std::map<std::pair<double,int>,std::pair<double,double>>::iterator it_Efficiency_EnergySplit;
  for(it_Efficiency_EnergySplit = map_Efficiency_EnergySplit.begin(); it_Efficiency_EnergySplit!=map_Efficiency_EnergySplit.end();
      it_Efficiency_EnergySplit++)
  {
    vec_g_ROC_EnergySplit.at(map_EnergyToCurve[it_Efficiency_EnergySplit->first.first])
                        ->SetPoint(it_Efficiency_EnergySplit->first.second,it_Efficiency_EnergySplit->second.second,
                                   it_Efficiency_EnergySplit->second.first);    
  }

  TLine *l_1PBperYear = new TLine(0,5.09e-3,1,5.09e-3);
  l_1PBperYear->SetLineWidth(3);
  TText *t_1PBperYear = new TText(1,5.09e-3,"1PB/Year");
  t_1PBperYear->SetTextSize(0.038);
  TLine *l_2PBperYear = new TLine(0,0.01,1,0.01);
  l_2PBperYear->SetLineWidth(3);
  TText *t_2PBperYear = new TText(1,0.01,"2PB/Year");
  t_2PBperYear->SetTextSize(0.038);
  TLine *l_5PBperYear = new TLine(0,0.025,1,0.025);
  l_5PBperYear->SetLineWidth(3);
  TText *t_5PBperYear = new TText(1,0.025,"5PB/Year");
  t_5PBperYear->SetTextSize(0.038);

  c->SetLogy();
  vec_g_ROC_EnergySplit.at(0)->Draw("APL");
  for(unsigned int i = 1; i < vec_g_ROC_EnergySplit.size(); i++)
  {
    vec_g_ROC_EnergySplit.at(i)->Draw("PL");
  }
  vec_g_ROC_EnergySplit.at(0)->GetXaxis()->SetLimits(0,1);
  vec_g_ROC_EnergySplit.at(0)->GetYaxis()->SetRangeUser(3e-3, bkgdMax);
  vec_g_ROC_EnergySplit.at(0)->SetTitle("Background Rate vs. Efficiency for True ENu Bands");
  vec_g_ROC_EnergySplit.at(0)->GetXaxis()->SetTitle("Efficiency");
  vec_g_ROC_EnergySplit.at(0)->GetYaxis()->SetTitle("Background Rate, (Hz)");
  leg_ROC_EnergySplit->Draw();
  l_1PBperYear->Draw();
  l_2PBperYear->Draw();
  l_5PBperYear->Draw();
  t_1PBperYear->Draw();
  t_2PBperYear->Draw();
  t_5PBperYear->Draw();

  c->SaveAs("ROCEnergySplit.pdf");
  c->Clear();

  return true;
}


bool makePlot3(TFile *f_Input_Analyse, TCanvas *c)
{
  c->SetLogy(0);
  THStack *stk_RemBkgd = new THStack("stk_RemBkgd","Generator Makeup of Remaining Background");
  TLegend *leg_RemBkgd = new TLegend(0.105,0.65,0.355,0.895);
  leg_RemBkgd->SetTextSize(0.023);
  leg_RemBkgd->SetHeader("Overall Efficiency & 10kt Background Rate");
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TH1I *h_RemBkgd = (TH1I*)f_Input_Analyse->Get(Form("h_RemainingBkgdGenType%i",i));
    h_RemBkgd->SetLineColor(vec_Colors.at(i));
    stk_RemBkgd->Add(h_RemBkgd);
    leg_RemBkgd->AddEntry(h_RemBkgd, Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
  }
  stk_RemBkgd->Draw("NOSTACK");
  stk_RemBkgd->GetXaxis()->SetTitle("Generator Type");
  stk_RemBkgd->GetYaxis()->SetTitle("Number of Hits");
  leg_RemBkgd->Draw();
  c->SaveAs("RemainingBackgroundGenType.pdf");
  c->Clear();

  return true;
}


bool makePlot4(TFile *f_Input_Analyse, TCanvas *c)
{
  c->SetLogy(0);
  THStack *stk_Just1 = new THStack("stk_Just1","Number of Marley Hits in Clusters Tagged Marley");
  THStack *stk_Just2 = new THStack("stk_Just2","Number of Background Hits in Clusters Tagged Marley");
  THStack *stk_Just3 = new THStack("stk_Just3","Number of Marley Hits in Clusters Tagged Background");
  THStack *stk_Just4 = new THStack("stk_Just4","Fraction of Background Hits in Clusters Tagged Marley");
  THStack *stk_Just5 = new THStack("stk_Just5","Fraction of Marley Hits in Clusters Tagged Background");
  THStack *stk_Just6 = new THStack("stk_Just5","Fraction of Marley Hits in Clusters Tagged Marley");
  TLegend *leg_Just1 = new TLegend(0.6,0.65,0.895,0.895);
  TLegend *leg_Just2 = new TLegend(0.105,0.65,0.4,0.895);
  leg_Just1->SetTextSize(0.023);
  leg_Just1->SetHeader("Overall Efficiency & 10kt Background Rate");
  leg_Just2->SetTextSize(0.023);
  leg_Just2->SetHeader("Overall Efficiency & 10kt Background Rate");
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TH1I *h_Just1 = (TH1I*)f_Input_Analyse->Get(Form("h_NMarlHit_Marl%i",i));
    h_Just1->SetLineColor(vec_Colors.at(i));
    stk_Just1->Add(h_Just1);
    TH1I *h_Just2 = (TH1I*)f_Input_Analyse->Get(Form("h_NBkgdHit_Marl%i",i));
    h_Just2->SetLineColor(vec_Colors.at(i));
    stk_Just2->Add(h_Just2);
    TH1I *h_Just3 = (TH1I*)f_Input_Analyse->Get(Form("h_NMarlHit_Bkgd%i",i));
    h_Just3->SetLineColor(vec_Colors.at(i));
    stk_Just3->Add(h_Just3);
    TH1D *h_Just4 = (TH1D*)f_Input_Analyse->Get(Form("h_FracBkgdHit_Marl%i",i));
    h_Just4->SetLineColor(vec_Colors.at(i));
    stk_Just4->Add(h_Just4);
    TH1D *h_Just5 = (TH1D*)f_Input_Analyse->Get(Form("h_FracMarlHit_Bkgd%i",i));
    h_Just5->SetLineColor(vec_Colors.at(i));
    stk_Just5->Add(h_Just5);
    TH1D *h_Just6 = (TH1D*)f_Input_Analyse->Get(Form("h_FracMarlHit_Marl%i",i));
    h_Just6->SetLineColor(vec_Colors.at(i));
    stk_Just6->Add(h_Just6);
    leg_Just1->AddEntry(h_Just1, Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
    leg_Just2->AddEntry(h_Just1, Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
  }
  stk_Just1->Draw("NOSTACK");
  stk_Just1->GetXaxis()->SetTitle("Number of Clusters");
  stk_Just1->GetYaxis()->SetTitle("Number of Hits");
  leg_Just1->Draw();
  c->SaveAs("Justify_MarlInMarl.pdf");
  c->Clear();

  stk_Just2->Draw("NOSTACK");
  stk_Just2->GetXaxis()->SetTitle("Number of Clusters");
  stk_Just2->GetYaxis()->SetTitle("Number of Hits");
  leg_Just1->Draw();
  c->SaveAs("Justify_BkgdInMarl.pdf");
  c->Clear();

  stk_Just3->Draw("NOSTACK");
  stk_Just3->GetXaxis()->SetTitle("Number of Clusters");
  stk_Just3->GetYaxis()->SetTitle("Number of Hits");
  leg_Just1->Draw();
  c->SaveAs("Justify_MarlInBkgd.pdf");
  c->Clear();

  stk_Just4->Draw("NOSTACK");
  stk_Just4->GetXaxis()->SetTitle("Fraction of Hits");
  stk_Just4->GetYaxis()->SetTitle("Number of Clusters");
  leg_Just1->Draw();
  c->SaveAs("Justify_FracBkgdInMarl.pdf");
  c->Clear();

  stk_Just5->Draw("NOSTACK");
  stk_Just5->GetXaxis()->SetTitle("Fraction of Hits");
  stk_Just5->GetYaxis()->SetTitle("Number of Clusters");
  leg_Just1->Draw();
  c->SaveAs("Justify_FracMarlInBkgd.pdf");
  c->Clear();

  stk_Just6->Draw("NOSTACK");
  stk_Just6->GetXaxis()->SetTitle("Fraction of Hits");
  stk_Just6->GetYaxis()->SetTitle("Number of Clusters");
  leg_Just2->Draw();
  c->SaveAs("Justify_FracMarlInMarl.pdf");
  c->Clear();

  return true;
}


bool makePlot5(TFile *f_Input_Analyse, TCanvas *c)
{
  c->SetLogy(0);
  THStack *stk_NMarlClusters = new THStack("stk_NMarlClusters","Number of Clusters Tagged Marley per Event (Drift Window)");
  TLegend *leg_NMarlClusters = new TLegend(0.6,0.65,0.895,0.895);
  leg_NMarlClusters->SetTextSize(0.023);
  leg_NMarlClusters->SetHeader("Overall Efficiency & 10kt Background Rate");
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TH1I *h_NMarlClusters = (TH1I*)f_Input_Analyse->Get(Form("h_NMarlClustersPerEvent%i",i));
    h_NMarlClusters->SetLineColor(vec_Colors.at(i));
    stk_NMarlClusters->Add(h_NMarlClusters);
    leg_NMarlClusters->AddEntry(h_NMarlClusters, 
                                Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
  }
  stk_NMarlClusters->Draw("NOSTACK");
  stk_NMarlClusters->GetXaxis()->SetTitle("Number of Clusters/Event");
  stk_NMarlClusters->GetYaxis()->SetTitle("Number of Events");
  leg_NMarlClusters->Draw();
  c->SaveAs("NMarlClustersPerEvent.pdf");
  c->Clear();

  return true;
}


bool makePlot6(TFile *f_Input_Module, TCanvas *c)
{
  c->Clear();
  TH1D *h_TimeElapsed = (TH1D*)f_Input_Module->Get("h_TimeElapsed");
  h_TimeElapsed->Draw();
  h_TimeElapsed->SetTitle("Time Elapsed Between Receiving a Channel Ordered List of Hits to an Issued Trigger per Event");
  h_TimeElapsed->GetXaxis()->SetTitle("Time/Event, (ms)");
  h_TimeElapsed->GetYaxis()->SetTitle("Number of Events");
  c->SaveAs("TimeElapsed.pdf");

  return true;
}


bool makePlot7(TFile *f_Input_Analyse, TCanvas *c)
{
  c->SetLogy(1);
  THStack *stk_BkgdVEnergy = new THStack("stk_BkgdVEnergy","Background Rate vs. Background Cluster ADC Total");
  TLegend *leg_BkgdVEnergy = new TLegend(0.65,0.7,0.895,0.895);
  leg_BkgdVEnergy->SetTextSize(0.023);
  leg_BkgdVEnergy->SetHeader("Overall Efficiency & 10kt Background Rate");
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TH1I *h_BkgdVEnergy = (TH1I*)f_Input_Analyse->Get(Form("h_BkgdVEnergy%i",i));
    h_BkgdVEnergy->SetLineColor(vec_Colors.at(i));
    stk_BkgdVEnergy->Add(h_BkgdVEnergy);
    leg_BkgdVEnergy->AddEntry(h_BkgdVEnergy, 
                                Form("Eff: %.2f, Bkgd: %.2fHz",vec_EffAndBkgd.at(i).first, vec_EffAndBkgd.at(i).second), "L");
  }
  stk_BkgdVEnergy->Draw("NOSTACK");
  stk_BkgdVEnergy->GetXaxis()->SetTitle("ADC Total of Background Cluster, (ADC)");
  stk_BkgdVEnergy->GetYaxis()->SetTitle("Rate, (Hz)");
  leg_BkgdVEnergy->Draw();
  c->SaveAs("BkgdVEnergy.pdf");
  c->Clear();

  return true;
}


int main()
{
  TString s_FileName = "GH_SNMC";
  TFile *f_Input_Module  = new TFile("Module_"+s_FileName+".root","READ");
  TFile *f_Input_Analyse = new TFile("Analyse_"+s_FileName+".root","READ");
  TFile *f_Input_Reco    = new TFile("/Users/alexanderbooth/Documents/Work/Year1/SNTrigger/Samples/"+s_FileName+".root", "READ");
  TCanvas *c = new TCanvas("c","c", 800, 500);

  //CONFIGURATION DEFINITIONS AND BACKGROUND RATE & EFFICIENCY. ONE POINT PER CONFIGURATION.
  bool plot1 = makePlot1(f_Input_Module, f_Input_Analyse, f_Input_Reco, c);

  //EFFIFICIENCY VS ENERGY AND ROC ENERGY SPLIT.
  bool plot2 = makePlot2(f_Input_Module, f_Input_Analyse, c);

  //REMAINING BACKGROUND.
  bool plot3 = makePlot3(f_Input_Analyse, c);

  //PLOTS TO JUSTIFY CHOICE OF MARLEY CLUSTER CLASSIFICATION.
  bool plot4 = makePlot4(f_Input_Analyse, c);

  //NUMBER OF MARLEY CLUSTERS PER EVENT.
  bool plot5 = makePlot5(f_Input_Analyse, c);

  //COMPUTING TIME IN MILLISECONDS.
  bool plot6 = makePlot6(f_Input_Module, c);

  //BACKGROUND RATE AT GIVEN ADC CLUSTER SUM.
  bool plot7 = makePlot7(f_Input_Analyse, c);

  return 0;
}
