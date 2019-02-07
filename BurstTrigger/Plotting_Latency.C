#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include "/Users/alexanderbooth/Documents/Work/UsefulCode/colors.h"


int main()
{
  TString s_Filename = "GH_SNMC";
  std::vector<int> vec_Colors = getColors(2);
  ifstream inputFile;
  inputFile.open("Analyse_SNBurst_"+s_Filename+".txt");
  std::map<int,std::vector<double>> map_ConfigToQuantities;
  int config;
  double eff, bkgd, perMonth;
  while(inputFile >> config >> eff >> bkgd >> perMonth)
  {
    std::vector<double> vec_Quantities;
    vec_Quantities.push_back(eff);
    vec_Quantities.push_back(bkgd);
    vec_Quantities.push_back(perMonth);
    map_ConfigToQuantities[config] = vec_Quantities;
    vec_Quantities.clear();
  }

  TFile *f_Input  = new TFile("Analyse_Latency.root","READ");
  TFile *f_Sample = new TFile("TimeSample.root","READ");

  int burstMin = 1;
  int burstMax = 30e4;
  std::vector<TH1D*> vec_h_Latency;
  for(unsigned int i = 0; i < map_ConfigToQuantities.size(); i++)
  {
    TH1D *hTemp = (TH1D*)f_Input->Get(Form("h_Latency_Config%i",i));
    hTemp->SetLineColor(vec_Colors.at(i));
    int    lastNonZeroBin;
    double lastTriggerTime;
    for(int j = hTemp->GetSize()-2; j > 0; j--)
    {
      if(hTemp->GetBinContent(j)>0)
      {
        lastNonZeroBin = j;
        lastTriggerTime = hTemp->GetBinContent(j);
        break;
      }
    }
    for(int j = lastNonZeroBin; j < hTemp->GetSize()-1; j++)
    {
      hTemp->SetBinContent(j, lastTriggerTime);
    }
    vec_h_Latency.push_back(hTemp);
  }

  TCanvas *c_Latency   = new TCanvas("c_Latency","c_Latency", 800, 500);
  THStack *stk_Latency = new THStack("stk_Latency", "Time Taken to Accumulate the Required Number of Events to Burst Trigger. 1/Month Fake");
  TLegend *leg_Latency = new TLegend(0.55,0.65,0.88,0.88);
  leg_Latency->SetTextSize(0.023);
  leg_Latency->SetHeader("Individual Marley Efficiency, 10kt Background Rate & Cut");
  for(unsigned int i = 0; i < vec_h_Latency.size(); i++)
  {
    stk_Latency->Add(vec_h_Latency.at(i));
    leg_Latency->AddEntry(vec_h_Latency.at(i), 
                          Form("Eff: %.2f, Bkgd: %.2fHz, Cut: %.2f",
                                map_ConfigToQuantities[i].at(0), map_ConfigToQuantities[i].at(1), 
                                map_ConfigToQuantities[i].at(2)/map_ConfigToQuantities[i].at(0)), "L");
  }
  c_Latency->Draw();
  c_Latency->SetLogx();
  c_Latency->SetLogy();
  stk_Latency->Draw("NOSTACK");
  stk_Latency->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_Latency->GetYaxis()->SetTitle("Latency, (s)");
  leg_Latency->Draw();
  c_Latency->SaveAs("LatencyVEvents.pdf");


  TFile *f_Theory = new TFile("SNTheoryDistributions.root","READ"); 
  TF1   *f_EventsVSNDistance_10kt = (TF1*)f_Theory->Get("f_EventsVSNDistance_10kt_100kpc");
  double gradient  = f_EventsVSNDistance_10kt->GetParameter(0);
  double intercept = f_EventsVSNDistance_10kt->GetParameter(1); 
  TF1 *f_Inverse   = new TF1("f_Inverse", "TMath::Power(x/(TMath::Power(10,[0])),1/[1])", 1,30);
  f_Inverse->SetParameter(0, intercept);
  f_Inverse->SetParameter(1, gradient);
  double distanceLow = f_Inverse->Eval((double)burstMax);
  double distanceUp  = f_Inverse->Eval((double)burstMin);
  
  std::vector<TH1D*> vec_h_LatencyVDistance;
  for(unsigned int i = 0; i < map_ConfigToQuantities.size(); i++)
  {
    TH1D *h_LatencyVDistance = new TH1D(Form("h_LatencyVDistance_Config%i",i),"", burstMax-burstMin+1, distanceLow, distanceUp);   
    h_LatencyVDistance->SetLineColor(vec_Colors.at(i));
    for(int j = 1; j <= burstMax-burstMin+1; j++)
    {
      double distance   = h_LatencyVDistance->GetBinCenter(j);
      double nEvents    = f_EventsVSNDistance_10kt->Eval(distance);
      double nEventsBin = vec_h_Latency.at(i)->FindBin(nEvents);
      double latency    = vec_h_Latency.at(i)->GetBinContent(nEventsBin);

      h_LatencyVDistance->SetBinContent(j, latency);
    }
    vec_h_LatencyVDistance.push_back(h_LatencyVDistance);
  }

  TCanvas *c_LatencyVDistance = new TCanvas("c_LatencyVDistance","c_LatencyVDistance", 800, 500);
  THStack *stk_LatencyVDistance = new THStack("stk_LatencyVDistance", 
                                              "Time Taken to Accumulate the Required Number of Events to Burst Trigger. 1/Month Fake");
  TLegend *leg_LatencyVDistance = new TLegend(0.12,0.65,0.45,0.88);
  leg_LatencyVDistance->SetTextSize(0.023);
  leg_LatencyVDistance->SetHeader("Individual Marley Efficiency, 10kt Background Rate & Cut");
  for(unsigned int i = 0; i < vec_h_LatencyVDistance.size(); i++)
  {
    stk_LatencyVDistance->Add(vec_h_LatencyVDistance.at(i));
    leg_LatencyVDistance->AddEntry(vec_h_LatencyVDistance.at(i),
                          Form("Eff: %.2f, Bkgd: %.2fHz, Cut: %.2f",
                                map_ConfigToQuantities[i].at(0), map_ConfigToQuantities[i].at(1), 
                                map_ConfigToQuantities[i].at(2)/map_ConfigToQuantities[i].at(0)), "L");
  }
  c_LatencyVDistance->Draw();
  c_LatencyVDistance->SetLogx();
  c_LatencyVDistance->SetLogy();
  stk_LatencyVDistance->Draw("NOSTACK");
  stk_LatencyVDistance->GetXaxis()->SetTitle("Supernova Distance, (kpc)");
  stk_LatencyVDistance->GetYaxis()->SetTitle("Latency, (s)");
  leg_LatencyVDistance->Draw();
  c_LatencyVDistance->SaveAs("LatencyVDistance.pdf");

  TH1D *h_SNProbabilityVDistance = (TH1D*)f_Theory->Get("h_SNProbabilityVDistance");
  double minDist = h_SNProbabilityVDistance->GetXaxis()->GetXmin();
  double maxDist = h_SNProbabilityVDistance->GetXaxis()->GetXmax();

  std::map<int,std::vector<std::pair<double,double>>> map_ConfigToCoverageAndLatency;
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++)
  {
    double fraction = 1 - h_SNProbabilityVDistance->Integral(i,-1);
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    for(unsigned int j = 0; j < map_ConfigToQuantities.size(); j++)
    {
      double latency = vec_h_LatencyVDistance.at(j)->GetBinContent(vec_h_LatencyVDistance.at(j)->FindBin(distance));
      map_ConfigToCoverageAndLatency[j].push_back({fraction,latency});
    }
  }

  double minLatency = 100;
  double maxLatency = 0;
  std::vector<TGraph*> vec_g_LatencyVCoverage;
  for(unsigned int i = 0; i < map_ConfigToQuantities.size(); i++)
  {
    TGraph *g_LatencyVCoverage = new TGraph(map_ConfigToCoverageAndLatency[i].size());
    g_LatencyVCoverage->SetMarkerStyle(3);
    g_LatencyVCoverage->SetMarkerColor(vec_Colors.at(i));
    for(unsigned int j = 0;  j < map_ConfigToCoverageAndLatency[i].size(); j++)
    {
      g_LatencyVCoverage->SetPoint(j, map_ConfigToCoverageAndLatency[i].at(j).first, 
                                      map_ConfigToCoverageAndLatency[i].at(j).second);
      if(map_ConfigToCoverageAndLatency[i].at(j).second>maxLatency)
      {
        maxLatency = map_ConfigToCoverageAndLatency[i].at(j).second;
      }
      if(map_ConfigToCoverageAndLatency[i].at(j).second<minLatency)
      {
        minLatency = map_ConfigToCoverageAndLatency[i].at(j).second;
      }
    }
    vec_g_LatencyVCoverage.push_back(g_LatencyVCoverage);
  }

  TCanvas *c_LatencyVCoverage = new TCanvas("c_LatencyVCoverage","c_LatencyVCoverage", 800, 500);
  TLegend *leg_LatencyVCoverage = new TLegend(0.12,0.695,0.47,0.895);
  leg_LatencyVCoverage->SetTextSize(0.023);
  leg_LatencyVCoverage->SetHeader("Individual Marley Efficiency, 10kt Background Rate & Cut");
  vec_g_LatencyVCoverage.at(0)->SetTitle("Burst Trigger Latency for a Given Galactic Neighbourhood Coverage");
  vec_g_LatencyVCoverage.at(0)->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  vec_g_LatencyVCoverage.at(0)->GetYaxis()->SetTitle("Latency, (s)");
  vec_g_LatencyVCoverage.at(0)->Draw("AP");
  leg_LatencyVCoverage->AddEntry(vec_g_LatencyVCoverage.at(0),
                                 Form("Eff: %.2f, Bkgd: %.2fHz, Cut: %.2f",
                                 map_ConfigToQuantities[0].at(0), map_ConfigToQuantities[0].at(1), 
                                 map_ConfigToQuantities[0].at(2)/map_ConfigToQuantities[0].at(0)), "P");
  for(unsigned int i = 1; i < vec_g_LatencyVCoverage.size(); i++)
  {
    leg_LatencyVCoverage->AddEntry(vec_g_LatencyVCoverage.at(i),
                          Form("Eff: %.2f, Bkgd: %.2fHz, Cut: %.2f",
                                map_ConfigToQuantities[i].at(0), map_ConfigToQuantities[i].at(1), 
                                map_ConfigToQuantities[i].at(2)/map_ConfigToQuantities[i].at(0)), "P");
    vec_g_LatencyVCoverage.at(i)->Draw("P"); 
  }
  c_LatencyVCoverage->SetLogy();
  vec_g_LatencyVCoverage.at(0)->GetXaxis()->SetRangeUser(-0.05,1.05);
  vec_g_LatencyVCoverage.at(0)->GetYaxis()->SetRangeUser(minLatency, maxLatency+1);
  leg_LatencyVCoverage->Draw();
  c_LatencyVCoverage->SaveAs("LatencyVCoverage.pdf");

  TCanvas *c_TimeProfile = new TCanvas("c_TimeProfile","c_TimeProfile", 800, 500);
  gStyle->SetOptStat(0);
  TH1D *hFine = (TH1D*)f_Sample->Get("h_MarlTime");
  TH1D *hFine_Reweighted = (TH1D*)f_Sample->Get("h_MarlTime_Zero20ms");
  TH1D *hFine_Reweighted_Extended = (TH1D*)f_Sample->Get("h_MarlTime_Zero20ms_Extrap");
  TH1D *hFine_Reweighted_Extended_Fixed = (TH1D*)f_Sample->Get("h_MarlTime_Zero20ms_Extrap_3secs");

  hFine->SetTitle("Marley Time Profile, Event Normalised");
  hFine->GetXaxis()->SetTitle("Time, (s)");
  hFine->GetYaxis()->SetTitle("Events");
  hFine->GetYaxis()->SetTitleOffset(1.2);
  hFine->GetFunction("f_Cooling_Fit")->SetBit(TF1::kNotDraw);
  hFine->GetYaxis()->SetRangeUser(0,hFine->GetMaximum());
  hFine->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley.pdf");
  c_TimeProfile->Clear();

  hFine_Reweighted->SetTitle("Marley Time Profile, Event Normalised. No Events < 20ms");
  hFine_Reweighted->GetXaxis()->SetTitle("Time, (s)");
  hFine_Reweighted->GetYaxis()->SetTitle("Events");
  hFine_Reweighted->GetYaxis()->SetTitleOffset(1.2);
  hFine_Reweighted->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley_ZeroTil20ms.pdf");
  c_TimeProfile->Clear();

  TH1D *hFine_First50ms = (TH1D*)hFine->Clone();
  hFine_First50ms->GetXaxis()->SetRange(hFine->FindBin(-0.001),hFine->FindBin(0.05));
  hFine_First50ms->SetTitle("Marley Time Profile, Event Normalised. First 50ms");
  hFine_First50ms->GetXaxis()->SetTitle("Time, (s)");
  hFine_First50ms->GetYaxis()->SetTitle("Events");
  hFine_First50ms->GetYaxis()->SetTitleOffset(1.2);
  hFine_First50ms->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley_First50ms.pdf");
  c_TimeProfile->Clear();

  TH1D *hFine_ReweightedFirst50ms = (TH1D*)hFine_Reweighted->Clone();
  hFine_ReweightedFirst50ms->GetXaxis()->SetRange(hFine_Reweighted->FindBin(-0.001),hFine_Reweighted->FindBin(0.05));
  hFine_ReweightedFirst50ms->SetTitle("Marley Time Profile, Event Normalised. First 50ms");
  hFine_ReweightedFirst50ms->GetXaxis()->SetTitle("Time, (s)");
  hFine_ReweightedFirst50ms->GetYaxis()->SetTitle("Events");
  hFine_ReweightedFirst50ms->GetYaxis()->SetTitleOffset(1.2);
  hFine_ReweightedFirst50ms->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley_ZeroTil20ms_First50ms.pdf");

  hFine_Reweighted_Extended->SetTitle("Marley Time Profile, Event Normalised. No Events < 20ms, Exponential Cooling to 120s. Fitted");
  hFine_Reweighted_Extended->GetXaxis()->SetTitle("Time, (s)");
  hFine_Reweighted_Extended->GetYaxis()->SetTitle("Events");
  hFine_Reweighted_Extended->GetYaxis()->SetTitleOffset(1.2);
  hFine_Reweighted_Extended->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley_ZeroTil20ms_Extrapolated120s.pdf");
  c_TimeProfile->Clear();

  hFine_Reweighted_Extended_Fixed->SetTitle(
                                  "Marley Time Profile, Event Normalised. No Events < 20ms, Exponential Cooling to 120s. Decay Time = 3secs");
  hFine_Reweighted_Extended_Fixed->GetXaxis()->SetTitle("Time, (s)");
  hFine_Reweighted_Extended_Fixed->GetYaxis()->SetTitle("Events");
  hFine_Reweighted_Extended_Fixed->GetYaxis()->SetTitleOffset(1.2);
  hFine_Reweighted_Extended_Fixed->Draw();
  c_TimeProfile->SaveAs("TimeProfile_Marley_ZeroTil20ms_Extrapolated120s_Fixed.pdf");
  c_TimeProfile->Clear();

  c_TimeProfile->Clear();
  hFine_First50ms->Divide(hFine_First50ms, hFine_ReweightedFirst50ms, 1, 1,"B");
  hFine_First50ms->SetTitle("Ratio of Time Profiles Before and After Weighting. First 50ms");
  hFine_First50ms->GetXaxis()->SetTitle("Time, (s)");
  hFine_First50ms->GetYaxis()->SetTitle("Ratio");
  hFine_First50ms->Draw();
  c_TimeProfile->SaveAs("RatioBeforeAndAfterWeighting_First50ms.pdf");


  //COUNTDOWN PLOTS.
  TGraph *g_TargetConfig = (TGraph*)f_Input->Get("g_TargetConfig");
  std::vector<double> vec_TargetConfig;
  for(unsigned int i = 0; i < g_TargetConfig->GetN(); i++)
  {
    double x, y;
    g_TargetConfig->GetPoint(i, x, y);
    vec_TargetConfig.push_back(y);
  }
  std::vector<TH1D*> vec_h_Latency_CountDown;
  for(unsigned int i = 0; i < vec_TargetConfig.size(); i++)
  {
    TH1D *hTemp = (TH1D*)f_Input->Get(Form("h_Latency_CountDown_Config%i",i));
    hTemp->SetLineColor(vec_Colors.at(i));
    int    lastNonZeroBin;
    double lastTriggerTime;
    for(int j = 1; j < hTemp->GetSize()-1; j++)
    {
      if(hTemp->GetBinContent(j)==0)
      {
        lastNonZeroBin = j-1;
        lastTriggerTime = hTemp->GetBinContent(j-1);
        break;
      }
    }
    for(int j = lastNonZeroBin; j < hTemp->GetSize()-1; j++)
    {
      if(lastTriggerTime!=0)
      {
        hTemp->SetBinContent(j, lastTriggerTime);
      }
    }
    vec_h_Latency_CountDown.push_back(hTemp);
  }

  TCanvas *c_Latency_CountDown   = new TCanvas("c_Latency_CountDown","c_Latency_CountDown", 800, 500);
  THStack *stk_Latency_CountDown = new THStack("stk_Latency_CountDown", "Time Taken to Collect all but a Small Fraction of Events for a Given Burst Size");
  TLegend *leg_Latency_CountDown = new TLegend(0.55,0.12,0.88,0.4);
  leg_Latency_CountDown->SetTextSize(0.023);
  for(unsigned int i = 0; i < vec_h_Latency_CountDown.size(); i++)
  {
    stk_Latency_CountDown->Add(vec_h_Latency_CountDown.at(i));
    leg_Latency_CountDown->AddEntry(vec_h_Latency_CountDown.at(i),Form("< %.2f Events Remaining",vec_TargetConfig.at(i)), "L");
  }
  c_Latency_CountDown->Draw();
  c_Latency_CountDown->SetLogx();
  c_Latency_CountDown->SetLogy();
  stk_Latency_CountDown->Draw("NOSTACK");
  stk_Latency_CountDown->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_Latency_CountDown->GetYaxis()->SetTitle("Time Since First Neutrino Passed Through Detector, (s)");
  leg_Latency_CountDown->Draw();
  c_Latency_CountDown->SaveAs("EventsRemainingVEvents.pdf");

  std::vector<TH1D*> vec_h_Latency_CountDownVDistance;
  for(unsigned int i = 0; i < vec_h_Latency_CountDown.size(); i++)
  {
    TH1D *h_Latency_CountDownVDistance = new TH1D(Form("h_Latency_CountDownVDistance_Config%i",i),"", burstMax-burstMin+1, distanceLow, distanceUp);   
    h_Latency_CountDownVDistance->SetLineColor(vec_Colors.at(i));
    for(int j = 1; j <= burstMax-burstMin+1; j++)
    {
      double distance   = h_Latency_CountDownVDistance->GetBinCenter(j);
      double nEvents    = f_EventsVSNDistance_10kt->Eval(distance);
      double nEventsBin = vec_h_Latency_CountDown.at(i)->FindBin(nEvents);
      double latency    = vec_h_Latency_CountDown.at(i)->GetBinContent(nEventsBin);

      h_Latency_CountDownVDistance->SetBinContent(j, latency);
    }
    vec_h_Latency_CountDownVDistance.push_back(h_Latency_CountDownVDistance);
  }

  TCanvas *c_Latency_CountDownVDistance = new TCanvas("c_Latency_CountDownVDistance","c_Latency_CountDownVDistance", 800, 500);
  THStack *stk_Latency_CountDownVDistance = new THStack("stk_Latency_CountDownVDistance", 
                                              "Time Taken to Collect all but a Small Fraction of Events for a Given SN Distance");
  TLegend *leg_Latency_CountDownVDistance = new TLegend(0.12,0.12,0.45,0.4);
  leg_Latency_CountDownVDistance->SetTextSize(0.023);
  for(unsigned int i = 0; i < vec_h_Latency_CountDownVDistance.size(); i++)
  {
    stk_Latency_CountDownVDistance->Add(vec_h_Latency_CountDownVDistance.at(i));
    leg_Latency_CountDownVDistance->AddEntry(vec_h_Latency_CountDownVDistance.at(i),
                                             Form("< %.2f Events Remaining",vec_TargetConfig.at(i)), "L");
  }
  c_Latency_CountDownVDistance->Draw();
  c_Latency_CountDownVDistance->SetLogx();
  c_Latency_CountDownVDistance->SetLogy();
  stk_Latency_CountDownVDistance->Draw("NOSTACK");
  stk_Latency_CountDownVDistance->GetXaxis()->SetTitle("Supernova Distance, (kpc)");
  stk_Latency_CountDownVDistance->GetYaxis()->SetTitle("Time Since First Neutrino Passed Through Detector, (s)");
  leg_Latency_CountDownVDistance->Draw();
  c_Latency_CountDownVDistance->SaveAs("EventsRemainingVDistance.pdf");

  std::map<int,std::vector<std::pair<double,double>>> map_ConfigToCoverageAndTime;
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++)
  {
    double fraction = 1 - h_SNProbabilityVDistance->Integral(i,-1);
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    for(unsigned int j = 0; j < vec_h_Latency_CountDownVDistance.size(); j++)
    {
      double time = vec_h_Latency_CountDownVDistance.at(j)->GetBinContent(vec_h_Latency_CountDownVDistance.at(j)->FindBin(distance));
      map_ConfigToCoverageAndTime[j].push_back({fraction,time});
    }
  }

  double minTime = 100;
  double maxTime = 0;
  std::vector<TGraph*> vec_g_Time_CountDownVCoverage;
  for(unsigned int i = 0; i < vec_TargetConfig.size(); i++)
  {
    TGraph *g_Time_CountDownVCoverage = new TGraph(map_ConfigToCoverageAndTime[i].size());
    g_Time_CountDownVCoverage->SetMarkerStyle(3);
    g_Time_CountDownVCoverage->SetMarkerColor(vec_Colors.at(i));
    for(unsigned int j = 0;  j < map_ConfigToCoverageAndTime[i].size(); j++)
    {
      if(map_ConfigToCoverageAndTime[i].at(j).first!=0)
      {
        g_Time_CountDownVCoverage->SetPoint(j, map_ConfigToCoverageAndTime[i].at(j).first, 
            map_ConfigToCoverageAndTime[i].at(j).second);
        if(map_ConfigToCoverageAndTime[i].at(j).second>maxTime)
        {
          maxTime = map_ConfigToCoverageAndTime[i].at(j).second;
        }
        if(map_ConfigToCoverageAndTime[i].at(j).second<minTime)
        {
          minTime = map_ConfigToCoverageAndTime[i].at(j).second;
        }
      }
    }
    vec_g_Time_CountDownVCoverage.push_back(g_Time_CountDownVCoverage);
  }

  TCanvas *c_Latency_CountDownVCoverage   = new TCanvas("c_Latency_CountDownVCoverage","c_Latency_CountDownVCoverage", 800, 500);
  TLegend *leg_Latency_CountDownVCoverage = new TLegend(0.12,0.12,0.5,0.42);
  leg_Latency_CountDownVCoverage->SetTextSize(0.023);
  vec_g_Time_CountDownVCoverage.at(0)->SetTitle("Time Taken to Collect all but a Small Fraction of Events for a Given Galactic Neighbourhood Coverage");
  vec_g_Time_CountDownVCoverage.at(0)->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  vec_g_Time_CountDownVCoverage.at(0)->GetYaxis()->SetTitle("Time Since First Neutrino Passed Through Detector, (s)");
  vec_g_Time_CountDownVCoverage.at(0)->Draw("AP");
  leg_Latency_CountDownVCoverage->AddEntry(vec_g_Time_CountDownVCoverage.at(0),Form("< %.2f Events Remaining",vec_TargetConfig.at(0)), "P");
  for(unsigned int i = 1; i < vec_g_Time_CountDownVCoverage.size(); i++)
  {
    leg_Latency_CountDownVCoverage->AddEntry(vec_g_Time_CountDownVCoverage.at(i),
                                             Form("< %.2f Events Remaining",vec_TargetConfig.at(i)), "P");
    vec_g_Time_CountDownVCoverage.at(i)->Draw("P"); 
  }
  vec_g_Time_CountDownVCoverage.at(0)->GetXaxis()->SetRangeUser(-0.05,1.05);
  vec_g_Time_CountDownVCoverage.at(0)->GetYaxis()->SetRangeUser(minTime-0.2, maxTime+0.2);
  leg_Latency_CountDownVCoverage->Draw();
  c_Latency_CountDownVCoverage->SaveAs("EventsRemainingVCoverage.pdf");

  return 0;
}
