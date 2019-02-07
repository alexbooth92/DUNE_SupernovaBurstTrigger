#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TGraph.h>


std::map<double,int> makeDistanceToEventsNeighbourhoodMap(TH1D* const &h_SNProbabilityVDistance, TF1* const &f_EventsVSNDistance_10kt)
{
  std::map<double,int> map_DistanceToEvents_Neighbourhood;
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++)
  {
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    int    nEvents  = std::floor(f_EventsVSNDistance_10kt->Eval(distance));
    map_DistanceToEvents_Neighbourhood[distance] = nEvents;
  }

  return map_DistanceToEvents_Neighbourhood;
}


void makeFakeRateVNClusters(double const &mean, double const &rms, double const &bkgd, int const &config, 
                            std::map<double,double> &map_FakeRateToNClusters)
{
  int    nBins  = 200; 
  double rmsMax = 8*rms;
  TH1D  *h_FakeRateVNClusters = new TH1D(Form("h_FakeRateVNClusters_Config%i",config),
                                         Form("h_FakeRateVNClusters_Config%i",config), nBins, mean, mean+rmsMax);

  for(unsigned int i = 1; i < h_FakeRateVNClusters->GetSize()-1; i++)
  {
    double nClusters = h_FakeRateVNClusters->GetBinCenter(i);
    double nRMS      = (nClusters-mean)/rms;
    double fraction  = (1-TMath::Erf(nRMS/std::sqrt(2.)))/2.;
    double fakeRate  = fraction*bkgd;

    h_FakeRateVNClusters->SetBinContent(i, fakeRate);
    map_FakeRateToNClusters[fakeRate] = nClusters;
  }

  h_FakeRateVNClusters->Write();
  delete h_FakeRateVNClusters;

  return;
}


std::map<std::pair<double,int>,double>  makeEfficiencyVEvents(double const &eff,       double const &mean,     double const &rms, 
                                                              double const &nClusters, int    const &burstMin, int    const &burstMax,
                                                              std::map<double,int> &map_ClustersToMaxEffEvent) 
{
  std::map<std::pair<double,int>,double> map_NClustersAndEventsToBurstEff;
  bool isOne = false;

  for(int i = burstMin; i <= burstMax; i++)
  {
    double integral      = 0;
    int    totalClusters = (int)std::ceil(mean) + i;
    TF1 *f_Poisson = new TF1("f_Poisson", "TMath::PoissonI(x,[0])", 0, mean + 20*rms);
    f_Poisson->SetParameter(0,totalClusters);

    if(isOne == false)
    {
      integral = 1. - f_Poisson->Integral(0, nClusters/eff);

      if(integral>=0 && integral<=1)
      {
        map_NClustersAndEventsToBurstEff[{nClusters,i}] = integral; 
        if(integral==1)
        {
          map_ClustersToMaxEffEvent[nClusters] = i;
          break;
        }
      }
      else if(integral<0)
      {
        map_NClustersAndEventsToBurstEff[{nClusters,i}] = 0; 
      }
    }

    delete f_Poisson;
  }

  return map_NClustersAndEventsToBurstEff;
}


void makeNeighbourhoodEfficiency(TH1D* const &h_SNProbabilityVDistance, std::map<double,int> &map_DistanceToEvents_Neighbourhood, 
                                 std::map<std::pair<double,int>,double> &map_NClustersAndEventsToBurstEff, double const &nClusters,
                                 TH1D* &h_NeighbourhoodEffiency, std::map<double,int> &map_ClustersToMaxEffEvent, 
                                 std::map<double,double> &map_ClustersToCoverage)
{
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++)
  {
    double prob     = h_SNProbabilityVDistance->GetBinContent(i);
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    int    nEvents  = map_DistanceToEvents_Neighbourhood[distance];
    double probXEff, burstEff;
    if(nEvents>=map_ClustersToMaxEffEvent[nClusters])
    {
      probXEff = prob*1; 
    }
    else
    {
      burstEff = map_NClustersAndEventsToBurstEff[{nClusters,nEvents}];
      probXEff = prob*burstEff;
    }
    h_NeighbourhoodEffiency->SetBinContent(i,probXEff);
  }

  double integral = h_NeighbourhoodEffiency->Integral();
  map_ClustersToCoverage[nClusters] = integral;

  return;
}


int main()
{
  //GRAB INFORMATION FROM THE CLUSTERING ALGORITHM AND DEFINE OUTPUT FILES.
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_"+s_Filename+".txt");
  ofstream outFile;
  outFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  //DEFINE PARAMETERS.
  int    burstMin(1), burstMax(30e4);
  double timeWindow   = 10; 
  double cut_PerMonth = 4.13e-7; 

  //GRAB THEORY PLOTS.
  TFile *f_Theory = new TFile("SNTheoryDistributions.root","READ");
  TFile *f_Output = new TFile("Analyse_SNBurst_"+s_Filename+".root", "RECREATE");

  TH1D *h_SNProbabilityVDistance = (TH1D*)f_Theory->Get("h_SNProbabilityVDistance");
  double min_GalacticNeighbourhood = h_SNProbabilityVDistance->GetXaxis()->GetXmin();
  double max_GalacticNeighbourhood = h_SNProbabilityVDistance->GetXaxis()->GetXmax();

  TF1  *f_EventsVSNDistance_10kt = (TF1*)f_Theory->Get("f_EventsVSNDistance_10kt_100kpc");
  double gradient  = f_EventsVSNDistance_10kt->GetParameter(0);
  double intercept = f_EventsVSNDistance_10kt->GetParameter(1); 
  TF1 *f_Inverse   = new TF1("f_Inverse", "TMath::Power(x/(TMath::Power(10,[0])),1/[1])", 1,40e4);
  f_Inverse->SetParameter(0, intercept);
  f_Inverse->SetParameter(1, gradient);
  double min_Distance = f_Inverse->Eval(burstMax);
  double max_Distance = f_Inverse->Eval(burstMin);
  std::map<double,int> map_DistanceToEvents_Neighbourhood = makeDistanceToEventsNeighbourhoodMap(h_SNProbabilityVDistance, 
                                                                                                 f_EventsVSNDistance_10kt);

  std::map<int,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int Config;
  double Eff, Bkgd;
  while(inFile >> Config >> Eff >> Bkgd)
  {
    map_ConfigToEffAndBkgd[Config] = {Eff,Bkgd};
  }

  //LOOP AROUND THE CLUSTERING CONFIGURATIONS.
  std::map<int,std::pair<double,double>>::iterator it_ConfigToEffAndBkgd;
  for(it_ConfigToEffAndBkgd=map_ConfigToEffAndBkgd.begin(); it_ConfigToEffAndBkgd!=map_ConfigToEffAndBkgd.end(); it_ConfigToEffAndBkgd++)
  {
    int    config = it_ConfigToEffAndBkgd->first;
    double eff    = it_ConfigToEffAndBkgd->second.first;
    double bkgd   = it_ConfigToEffAndBkgd->second.second;

    //ASSUME THE NUMBER OF CLUSTERS IN A GIVEN TIME WINDOW IS GAUSSIAN ABOUT THE GIVEN BACKGROUND RATE WITH AN RMS OF SQRT(MEAN).
    double mean   = bkgd*timeWindow;
    double rms    = std::sqrt(mean);

    std::map<double,double> map_FakeRateToNClusters;
    makeFakeRateVNClusters(mean, rms, bkgd, config, map_FakeRateToNClusters);

    //PULL OUT THE PER MONTH CLUSTER CUT FOR THIS CONFIGURATION.
    double cut_KeyNClustersForOnePerMonth = 0.;
    double cut_NClustersForOnePerMonth    = 0.;
    double lastKey, lastN;
    std::map<double,double>::iterator it_FakeRateToNClusters;
    for(it_FakeRateToNClusters=map_FakeRateToNClusters.begin();it_FakeRateToNClusters!=map_FakeRateToNClusters.end();it_FakeRateToNClusters++)
    {
      if(it_FakeRateToNClusters->first >= cut_PerMonth)
      {
        cut_KeyNClustersForOnePerMonth = lastKey;
        cut_NClustersForOnePerMonth    = lastN;
        break;
      }
      lastKey = it_FakeRateToNClusters->first;
      lastN   = it_FakeRateToNClusters->second;
    }

    std::cout << "CONFIG: " << config << ", EFF: " << eff << ", BKGD: " << bkgd << ", PERMONTH CUT:  " 
              << cut_NClustersForOnePerMonth << std::endl;
    outFile   << config << " " << eff << " " << bkgd << " " << cut_NClustersForOnePerMonth << std::endl;

    //LOOP OVER THE DIFFERENT FAKE RATES AND CUTS TO GET FAKE RATE AGAINST GALACTIC COVERAGE.
    std::map<std::pair<double,int>,double> map_NClustersAndEventsToBurstEff;
    std::map<double,double> map_ClustersToCoverage;
    std::map<double,int>    map_ClustersToMaxEffEvent;
    int count_Loop = 0;
    for(it_FakeRateToNClusters=map_FakeRateToNClusters.begin();it_FakeRateToNClusters!=map_FakeRateToNClusters.end();it_FakeRateToNClusters++)
    {
      if(count_Loop % 50 == 0)
      {
        std::cout << "WORKING ON CONFIG: " << config << ", ITERATION: " << count_Loop << std::endl;
      }
      TH1D *h_NeighbourhoodEffiency = (TH1D*)h_SNProbabilityVDistance->Clone();
      h_NeighbourhoodEffiency->Clear();
      h_NeighbourhoodEffiency->SetName(Form("h_NeighbourhoodEffiency_Config%i",config));

      //FOR EACH BURST SIZE
      if(count_Loop % 1 == 0 || (it_FakeRateToNClusters->first == cut_KeyNClustersForOnePerMonth 
      && it_FakeRateToNClusters->second == cut_NClustersForOnePerMonth))
      {
        map_NClustersAndEventsToBurstEff = makeEfficiencyVEvents(eff, mean, rms, it_FakeRateToNClusters->second, burstMin, burstMax,
                                                                 map_ClustersToMaxEffEvent);
        makeNeighbourhoodEfficiency(h_SNProbabilityVDistance, map_DistanceToEvents_Neighbourhood,
                                    map_NClustersAndEventsToBurstEff, it_FakeRateToNClusters->second, 
                                    h_NeighbourhoodEffiency, map_ClustersToMaxEffEvent, map_ClustersToCoverage);
      }

      //MAKE THE 1 PER MONTH HISTOGRAMS.
      if(it_FakeRateToNClusters->first == cut_KeyNClustersForOnePerMonth && it_FakeRateToNClusters->second == cut_NClustersForOnePerMonth)
      {
        TH1D *h_EfficiencyVEvents   = new TH1D(Form("h_EfficiencyVEvents_Config%i",config),Form("h_EfficiencyVEvents_Config%i",config),
                                                  burstMax-burstMin+1, burstMin, burstMax);
        for(int i = burstMin; i <= burstMax; i++)
        {
          int bin = h_EfficiencyVEvents->FindBin(i);
          h_EfficiencyVEvents->SetBinContent(bin,map_NClustersAndEventsToBurstEff[{cut_NClustersForOnePerMonth,i}]);
        }
        h_EfficiencyVEvents->Write();

        TH1D *h_EfficiencyVDistance = new TH1D(Form("h_EfficiencyVDistance_Config%i",config),Form("h_EfficiencyVDistance_Config%i",config),
                                               burstMax-burstMin+1, min_Distance, max_Distance); 
        for(int i = 1; i < h_EfficiencyVDistance->GetSize()-1; i++)
        {
          double distance   = h_EfficiencyVDistance->GetBinCenter(i);
          double nEvents    = f_EventsVSNDistance_10kt->Eval(distance);
          int    nEventsBin = h_EfficiencyVEvents->FindBin(nEvents);
          double burstEff   = h_EfficiencyVEvents->GetBinContent(nEventsBin);

          h_EfficiencyVDistance->SetBinContent(i,burstEff);
        }
        h_EfficiencyVDistance->Write();
        delete h_EfficiencyVEvents;
        delete h_EfficiencyVDistance;

        h_NeighbourhoodEffiency->Write();
      }
      delete h_NeighbourhoodEffiency;

      count_Loop++;
    }

    //EXTRACT THE INFORMATION WE HAVE PICKED UP AND MAKE THE ROC PLOT.
    TGraph *g_ROC = new TGraph(map_FakeRateToNClusters.size());
    g_ROC->SetName(Form("g_ROC_Config%i",config));
    count_Loop = 0;
    for(it_FakeRateToNClusters=map_FakeRateToNClusters.begin();it_FakeRateToNClusters!=map_FakeRateToNClusters.end();it_FakeRateToNClusters++)
    {
      if(map_ClustersToCoverage.count(it_FakeRateToNClusters->second)==1)
      {
        g_ROC->SetPoint(count_Loop,map_ClustersToCoverage[it_FakeRateToNClusters->second], it_FakeRateToNClusters->first);
        count_Loop++;
      }
    }
    g_ROC->Write();
  }

  return 0;
}
