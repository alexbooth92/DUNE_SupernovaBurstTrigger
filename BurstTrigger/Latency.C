#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TText.h>
#include <TRandom.h>


int main()
{
  TString s_Filename = "GH_SNMC";
  TFile *f_Input  = new TFile("TimeSample.root","READ");
  TFile *f_Output = new TFile("Analyse_Latency.root","RECREATE");
  TH1D  *h_Sample = (TH1D*)f_Input->Get("h_MarlTime_Zero20ms_Extrap");

  int burstMin = 1;
  int burstMax = 3e5; 

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

  int    sampleEvery   = 1;
  int    startBin      = 0;
  double startTime     = 0;
  int    endBin        = 0;
  double endTime       = 0;
  double snDuration    = 0;
  double secondsPerBin = 0;
  double sensitivity   = 0;
  for(int i = 1; i < h_Sample->GetSize()-1; i++)
  {
    if(h_Sample->GetBinContent(i)>0)
    {
      startBin  = i;
      startTime = h_Sample->GetBinCenter(i);
      break;
    }
  }
  for(int i = h_Sample->GetSize()-2; i > 0; i--)
  {
    if(h_Sample->GetBinContent(i)>0)
    {
      endBin  = i;
      endTime = h_Sample->GetBinCenter(i);
      break;
    }
  }
  snDuration    = endTime - startTime;
  secondsPerBin = h_Sample->GetBinWidth(0);
  sensitivity   = (double)sampleEvery*secondsPerBin;

  std::cout << "THE SNs START AT TIME t = " << startTime     << ", BIN " << startBin << std::endl;
  std::cout << "THE SNs END AT TIME t = "   << endTime       << ", BIN " << endBin << std::endl;
  std::cout << "THE SNs DURATION IS "       << snDuration    << "s"      << std::endl;
  std::cout << "THE SECONDS PER BIN IS "    << secondsPerBin << "s"      << std::endl;
  std::cout << "THE SENSITIVITY IS "        << sensitivity   << "s"      << std::endl;

  std::map<std::pair<int,int>,double> map_ConfigAndEventsToTime;
  //LOOP OVER THE CLUSTERING CONFIGURATIONS AND THEREFORE MONTHLY CUTS.
  //COUNTING UPWARDS.
  std::cout << "WORKING ON COUNTING UP" << std::endl;
  for(unsigned int i = 0; i < map_ConfigToQuantities.size(); i++)
  {
    bool   sensitivityReached = false;
    double cut_perMonth = map_ConfigToQuantities[i].at(2)/map_ConfigToQuantities[i].at(0);
    //LOOP OVER THE DIFFERENT SIZES OF SUPERNOVA BURSTS. ASSUMING THAT THERE WILL NEVER BE A TRIGGER IF THE TOTAL NUMBER OF EVENTS IN THE
    for(int j = std::floor(cut_perMonth); j <= burstMax; j++)
    {
      int    nEventsInSN      = j;
      bool   switch_passedCut = false;
      double triggerTime      = 0;
      for(int k = startBin; k < h_Sample->GetSize()-1; k++)
      {
        if(k % sampleEvery == 0)
        {
          double eventCount = (1 - h_Sample->Integral(k,-1))*nEventsInSN; 
          if(eventCount>cut_perMonth && switch_passedCut == false)
          {
            triggerTime = h_Sample->GetBinCenter(k)-startTime;
            map_ConfigAndEventsToTime[{i,nEventsInSN}] = triggerTime; 
            switch_passedCut = true;
            if(triggerTime<sensitivity)
            {
              std::cout << "SENSITIVITY REACHED! BURST "  << nEventsInSN 
                        << " TRIGGER TIME " << triggerTime << " NOW LOWER THAN THE SENSITIVITY " << sensitivity << std::endl; 
              sensitivityReached = true;
              break;
            }
            break;
          }
        }
      }
      if(sensitivityReached==true)
      {
        break;
      }
      if(j % 5000 == 0)
      {
        std::cout << "CONFIG: " << i << ", NUMBER OF EVENTS IN BURST: " << nEventsInSN << std::endl;  
        std::cout << "LAST TRIGGER TIME: " << triggerTime << " SENSITIVITY: " << sensitivity << std::endl;
      }
    }
  }

  gFile = f_Output;
  for(unsigned int i = 0; i < map_ConfigToQuantities.size(); i++)
  {
    TH1D *h_Latency = new TH1D(Form("h_Latency_Config%i",i),"",burstMax-burstMin+1,burstMin, burstMax);
    for(unsigned int j = burstMin; j < burstMax; j++)
    {
      int    bin = h_Latency->FindBin(j);
      h_Latency->SetBinContent(bin,map_ConfigAndEventsToTime[{i,j}]);
    }
    h_Latency->Write();
    delete h_Latency;
  }
  map_ConfigAndEventsToTime.clear();

  std::cout << "WORKING ON COUNTING DOWN" << std::endl;
  //CONFIG FROM NOW ON REFERS TO THE INT CORRESPODING TO THE TARGET NUMBER OF EVENTS.
  std::vector<double> vec_TargetEvents = {1, 0.5, 0.2, 0.05};
  TGraph *g_TargetConfig = new TGraph(vec_TargetEvents.size());
  g_TargetConfig->SetName("g_TargetConfig");
  for(unsigned int i = 0; i < vec_TargetEvents.size(); i++)
  {
    g_TargetConfig->SetPoint(i,i,vec_TargetEvents.at(i));
  }
  g_TargetConfig->Write();

  std::map<int,double> map_BinToProfileFraction_BinToEnd;
  for(int i = startBin; i <= endBin; i++)
  {
    double profileFraction = h_Sample->Integral(i,-1);
    map_BinToProfileFraction_BinToEnd[i] = profileFraction;
  }

  for(unsigned int i = 0; i < vec_TargetEvents.size(); i++)
  {
    double targetEvents = vec_TargetEvents.at(i);
    int    nEventsInSN  = burstMin;
    for(int j = startBin; j <= endBin; j++)
    {
      double nEventsRemaining = map_BinToProfileFraction_BinToEnd[j]*nEventsInSN;
      if(j % 5000 == 0)
      {
        std::cout << "CONFIG: "       << i           << ", THE BIN IS: "        << j 
                  << ", BURST SIZE: " << nEventsInSN << ", NEVENTS REMAINING: " << nEventsRemaining << std::endl; 
      }
      if(nEventsInSN>burstMax)
      {
        break;
      }
      if(nEventsRemaining<=targetEvents)
      {
        double time = h_Sample->GetBinCenter(j) - startTime;
        map_ConfigAndEventsToTime[{i,nEventsInSN}] = time; 
        nEventsInSN++;
        j--;
      }
    }
  }

  for(unsigned int i = 0; i < vec_TargetEvents.size(); i++)
  {
    TH1D *h_Latency = new TH1D(Form("h_Latency_CountDown_Config%i",i),"",burstMax-burstMin+1,burstMin, burstMax);
    for(unsigned int j = burstMin; j < burstMax; j++)
    {
      int    bin = h_Latency->FindBin(j);
      h_Latency->SetBinContent(bin,map_ConfigAndEventsToTime[{i,j}]);
    }
    h_Latency->Write();
    delete h_Latency;
  }
  return 0;
}
