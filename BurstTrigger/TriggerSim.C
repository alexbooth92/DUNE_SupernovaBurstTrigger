#include <iostream>
#include <vector>
#include <utility>
#include <TRandom.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <algorithm>
#include <TGraph.h>
#include <TString.h>
#include <TLine.h>
#include <TText.h>
#include <TMath.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>

int main()
{
  //BACKGROUND RATE IN HERTZ, TIME WINDOW AND SIMULATION TIME IN SECONDS.
  double detectorScaling = 0.12;
  double bkgdRate(2.66), timeWindow(10), simulationTime(1e4), sampleEvery(1);
  bkgdRate/=detectorScaling;

  TFile *f_Output = new TFile("TriggerSim.root", "UPDATE");
  gFile = f_Output;

  double nRandom(bkgdRate*simulationTime), startTime(0.), endTime(simulationTime), seed(1.);
  gRandom->SetSeed(seed);

  TH1D *h_Rand = new TH1D("h_Rand", "h_Rand", 100, startTime, endTime);
  TH1D *h_Rand_Zoomed = new TH1D("h_Rand_Zoomed", "h_Rand_Zoomed", 30, 0, 300);
  std::vector<double> vec_Times;

  std::cout << "GENERATING RANDOM CLUSTERS" << std::endl; 

  for(long long i = 0; i < nRandom; i++)
  {
    if(i % 500000 == 0)
    {
      std::cout << "GENERATING RANDOM NUMBER " << i << std::endl;
    }
    double randomTime = gRandom->Uniform(startTime, endTime);
    h_Rand->Fill(randomTime);
    h_Rand_Zoomed->Fill(randomTime);
    vec_Times.push_back(randomTime);
  }

  std::cout << "SORTING" << std::endl;
  std::sort(vec_Times.begin(), vec_Times.end());

  TH1D *h_ClustersInTimeWindow      = new TH1D("h_ClustersInTimeWindow","h_ClustersInTimeWindow",50,
                                               bkgdRate*timeWindow-3*std::sqrt(bkgdRate*timeWindow), 
                                               bkgdRate*timeWindow+3*std::sqrt(bkgdRate*timeWindow));
  TH1D *h_ClustersInTimeWindow_Full = new TH1D("h_ClustersInTimeWindow_Full","h_ClustersInTimeWindow_Full",40,
                                               bkgdRate*timeWindow-6*std::sqrt(bkgdRate*timeWindow), 
                                               bkgdRate*timeWindow+6*std::sqrt(bkgdRate*timeWindow));
  TH1D *h_StartTimes                = new TH1D("h_StartTimes", "h_StartTimes", 1001, -0.5, 1000.5); 
  long long nextTime = 0;

  long long second = 0;
  while(second < simulationTime)
  {
    h_StartTimes->Fill(second);
    double nEvents = 0;
    for(long long i = 0; i < vec_Times.size(); i++)
    {
      if(vec_Times.at(i)>=second && vec_Times.at(i)<second+timeWindow)
      {
        nEvents++;
      }
      else if(vec_Times.at(i)>=second+timeWindow)
      {
        break;
      }
    }
    h_ClustersInTimeWindow->Fill(nEvents);
    h_ClustersInTimeWindow_Full->Fill(nEvents);
    second+=sampleEvery;

    if(second % 50000 == 0)
    {
      std::cout << "WORKING ON SECOND: " << second << std::endl;
    }
  }

  h_ClustersInTimeWindow_Full->Fit("gaus");
  TF1 *f_Gaus = h_ClustersInTimeWindow_Full->GetFunction("gaus");

  TText *t_Fit = new TText(bkgdRate*timeWindow-5.5*std::sqrt(bkgdRate*timeWindow), h_ClustersInTimeWindow_Full->GetMaximum(),
                     Form("Gaussian Fit, Mean: %f, Sigma: %f", f_Gaus->GetParameter(1), f_Gaus->GetParameter(2))); 
  t_Fit->SetTextSize(0.03);

  h_Rand->Write();
  h_Rand_Zoomed->Write();
  h_ClustersInTimeWindow->Write();
  h_ClustersInTimeWindow_Full->Write();
  h_StartTimes->Write();

  std::cout << "FILES WRITTEN, SIMULATION OVER." << std::endl;

  TCanvas *c = new TCanvas("c","c", 800, 500);
  gStyle->SetOptStat(0);
  h_ClustersInTimeWindow_Full->SetTitle("Number of Clusters in Time Window");
  h_ClustersInTimeWindow_Full->GetXaxis()->SetTitle("Number of Clusters");
  h_ClustersInTimeWindow_Full->GetYaxis()->SetTitle("Number of Windows");
  h_ClustersInTimeWindow_Full->Draw();
  t_Fit->Draw();
  c->SaveAs("ClustersInTimeWindow.pdf");

  return 0;
}
