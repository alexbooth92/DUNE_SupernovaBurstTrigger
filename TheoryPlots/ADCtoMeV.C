#include <iostream>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TMath.h>


int main()
{
  TString s_FileName = "GH_SNMC";
  TFile *f_Sample = new TFile("/Users/alexanderbooth/Documents/Work/Year1/SNTrigger/Samples/"+s_FileName+".root","READ"); 
  TTree *t_Sample = (TTree*)f_Sample->Get("DAQSimTree");

  int Event_Sample;
  double Px_Sample;
  double Pz_Sample;
  double VertX_Sample;

  t_Sample->SetBranchAddress("Event", &Event_Sample);
  t_Sample->SetBranchAddress("Px",    &Px_Sample);
  t_Sample->SetBranchAddress("Pz",    &Pz_Sample);
  t_Sample->SetBranchAddress("VertX", &VertX_Sample);

  std::cout << "EXTRACTING ANGLES AND DISTANCES" << std::endl;
  std::map<int,std::pair<double,double>> map_EventToDistanceAndAngle;
  for(int i = 0; i < t_Sample->GetEntries(); i++)
  {
    t_Sample->GetEntry(i);
    double angle = TMath::Tan(Pz_Sample/Px_Sample);
    map_EventToDistanceAndAngle[Event_Sample] = {VertX_Sample,angle};
  }

  TFile *f_Input     = new TFile("/Users/alexanderbooth/Documents/Work/Year1/SNTrigger/WC_180212/Clustering_dunetpc/"+s_FileName+"/Module_"+s_FileName+".root","READ");
  TTree *t_Input  = (TTree*)f_Input->Get("t_Output");

  int Event;
  int Config;
  int Type;
  float SumADC;
  double ENu_Lep;

  int nClusters = t_Input->GetEntries();
  t_Input->SetBranchAddress("Event",  &Event  );
  t_Input->SetBranchAddress("Config", &Config );
  t_Input->SetBranchAddress("Type",   &Type   );
  t_Input->SetBranchAddress("SumADC", &SumADC );
  t_Input->SetBranchAddress("ENu_Lep",&ENu_Lep);

  std::cout << "EXTRACTING ENERGY INFORMATION" << std::endl;
  std::map<int,double> map_EventToTotalADC;
  std::map<int,double> map_EventToELep;
  for(int i = 0; i < nClusters; i++)
  {
    t_Input->GetEntry(i);
    if(Type==1 && Config == 4)
    {
      map_EventToTotalADC[Event] += SumADC;
      map_EventToELep[Event]      = ENu_Lep*1000; 
    }
  }

  TH2D *h_RatioVELep_D1 = new TH2D("h_RatioVELep_D1","h_RatioVELep_D1",30,0,60,30,20,1000);
  TH2D *h_RatioVELep_D2 = new TH2D("h_RatioVELep_D2","h_RatioVELep_D2",30,0,60,30,20,1000);
  TH2D *h_RatioVELep_D3 = new TH2D("h_RatioVELep_D3","h_RatioVELep_D3",30,0,60,30,20,1000);

  TH2D *h_RatioVELep_A1 = new TH2D("h_RatioVELep_A1","h_RatioVELep_A1",30,0,60,30,20,1000);
  TH2D *h_RatioVELep_A2 = new TH2D("h_RatioVELep_A2","h_RatioVELep_A2",30,0,60,30,20,1000);
  TH2D *h_RatioVELep_A3 = new TH2D("h_RatioVELep_A3","h_RatioVELep_A3",30,0,60,30,20,1000);

  std::cout << "FILLING HISTOGRAMS" << std::endl;
  std::map<int,double>::iterator it_EventToELep;
  for(it_EventToELep=map_EventToELep.begin(); it_EventToELep!=map_EventToELep.end(); it_EventToELep++)
  {
    if(std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].first) < 100)
    {
      h_RatioVELep_D1->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
    else if(std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].first) >= 100 
         && std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].first)  < 200)
    {
      h_RatioVELep_D2->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
    else    
    {
      h_RatioVELep_D3->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
    
    if(std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].second) < 10)
    {
      h_RatioVELep_A1->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
    else if(std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].second) >= 10 
         && std::abs(map_EventToDistanceAndAngle[it_EventToELep->first].second)  < 20)
    {
      h_RatioVELep_A2->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
    else
    {
      h_RatioVELep_A3->Fill(it_EventToELep->second, map_EventToTotalADC[it_EventToELep->first]/it_EventToELep->second);
    }
  }

  TFile *f_Output = new TFile("ADCtoMeV.root","RECREATE");

  h_RatioVELep_D1->Write();
  h_RatioVELep_D2->Write();
  h_RatioVELep_D3->Write();
  h_RatioVELep_A1->Write();
  h_RatioVELep_A2->Write();
  h_RatioVELep_A3->Write();

  return 0;
}
