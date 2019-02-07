#include <iostream>
#include <map>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TCut.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TBox.h>


TTree *t;
TCut cut_Cluster;
TCut cut_Event;
TCut cut_Config;
int Cluster;
int event;
int config;

int Config;
int Event;
int NChan;
int NHits;
int Type;
int StartChan;
int EndChan;
float FirstTimeHit;
float LastTimeHit;
float SumADC;

std::vector<int>   *HitView = 0;
std::vector<int>   *GenType = 0;
std::vector<int>   *HitChan = 0;
std::vector<float> *HitTime = 0;
std::vector<float> *HitSADC = 0;
std::vector<float> *HitRMS  = 0;

std::vector<int> vec_Event;
std::vector<int> vec_NChan;
std::vector<int> vec_NHits;
std::vector<int> vec_Type;
std::vector<int> vec_StartChan;
std::vector<int> vec_EndChan;
std::vector<float> vec_FirstTimeHit;
std::vector<float> vec_LastTimeHit;
std::vector<float> vec_SumADC;

std::vector<std::vector<int>>   vec_HitView;
std::vector<std::vector<int>>   vec_GenType;
std::vector<std::vector<int>>   vec_HitChan;
std::vector<std::vector<float>> vec_HitTime;
std::vector<std::vector<float>> vec_HitSADC;
std::vector<std::vector<float>> vec_HitRMS;


void drawHitChan()
{
  std::cout << "DRAWING CLUSTER: " << Cluster <<  " CHANNELS." << std::endl;
  t->Draw("HitChan", cut_Cluster);

  return;
}


void drawSingleCluster(bool marlOnly, bool makePDF)
{
  std::cout << "DRAWING CLUSTER: " << Cluster << std::endl;

  if(marlOnly == false)
  {
    TH2F *h_Cluster = new TH2F("h_Cluster", Form("Cluster %i of Type %i, Event %i, %i Total Hits.",Cluster,Type,Event, NHits), 
        std::abs(EndChan-StartChan)+11,StartChan-5, EndChan+5, 
        std::abs(LastTimeHit-FirstTimeHit)+201, FirstTimeHit-100, LastTimeHit+100);
    h_Cluster->GetXaxis()->SetTitle("Channel, (Channel No)");
    h_Cluster->GetYaxis()->SetTitle("Time, (tick)");

    for(int i = 0; i < NHits; i++)
    {
      int HitSpread      = std::ceil(2*HitRMS->at(i)); 
      float StartSpread  = HitTime->at(i)-(HitSpread/2.); 
      float EndSpread    = HitTime->at(i)+(HitSpread/2.);
      float ChargeSpread = HitSADC->at(i)/(float)HitSpread; 
      for(int j = 0; j < HitSpread; j++)
      {
        float TimePosition = StartSpread + j;
        h_Cluster->Fill(HitChan->at(i), TimePosition, ChargeSpread);
      }
    }

    gStyle->SetOptStat(0);
    TCanvas *c_Cluster = new TCanvas("c_Cluster", "c_Cluster", 800, 500);
    h_Cluster->SetContour(100);
    h_Cluster->Draw("COLZ");
    if(makePDF==true)
    {
      c_Cluster->SaveAs(Form("SingleCluster_%i_Event%i_AllHits.pdf", Cluster, Event));
    }
  }
  else
  {
    int   marl_StartChan = 10e5;
    int   marl_EndChan   = 0;
    float marl_FirstTimeHit = 5000;
    float marl_LastTimeHit  = 0; 

    std::vector<int> marleyHits;
    for(int i = 0; i < NHits; i++)
    {
      if(GenType->at(i)==1)
      {
        marleyHits.push_back(i);
        if(HitTime->at(i)<marl_FirstTimeHit)
        {
          marl_FirstTimeHit = HitTime->at(i);
        }
        if(HitTime->at(i)>marl_LastTimeHit)
        {
          marl_LastTimeHit = HitTime->at(i);
        }
        if(HitChan->at(i)<marl_StartChan)
        {
          marl_StartChan = HitChan->at(i);
        }
        if(HitChan->at(i)>marl_EndChan)
        {
          marl_EndChan = HitChan->at(i);
        }
      }
    }

    TH2F *h_Cluster = new TH2F("h_Cluster", Form("Cluster %i of Type %i, Event %i, %i Total Hits.",Cluster,Type,Event, marleyHits.size()), 
        std::abs(marl_EndChan-marl_StartChan)+5,marl_StartChan-2, marl_EndChan+2, 
        std::abs(marl_LastTimeHit-marl_FirstTimeHit)+61, marl_FirstTimeHit-30, marl_LastTimeHit+30);
    h_Cluster->GetXaxis()->SetTitle("Channel, (Channel No)");
    h_Cluster->GetYaxis()->SetTitle("Time, (tick)");

    for(int i = 0; i < marleyHits.size(); i++)
    {
      int HitSpread      = std::ceil(2*HitRMS->at(marleyHits.at(i))); 
      float StartSpread  = HitTime->at(marleyHits.at(i))-(HitSpread/2.); 
      float EndSpread    = HitTime->at(marleyHits.at(i))+(HitSpread/2.);
      float ChargeSpread = HitSADC->at(marleyHits.at(i))/(float)HitSpread; 
      for(int j = 0; j < HitSpread; j++)
      {
        float TimePosition = StartSpread + j;
        h_Cluster->Fill(HitChan->at(marleyHits.at(i)), TimePosition, ChargeSpread);
      }
    }

    gStyle->SetOptStat(0);
    TCanvas *c_Cluster = new TCanvas("c_Cluster", "c_Cluster", 800, 500);
    h_Cluster->SetContour(100);
    h_Cluster->Draw("COLZ");
    if(makePDF==true)
    {
      c_Cluster->SaveAs(Form("SingleCluster_%i_Event%i_MarlHits.pdf", Cluster, Event));
    }
  }

  return;
}


void drawAllClusters(bool marlOnly, bool makePDF)
{
  if(marlOnly == false)
  {
    std::vector<TBox*> vec_Box;
    int   all_StartChan = 10e5;
    int   all_EndChan   = 0;
    float all_FirstTimeHit = 5000;
    float all_LastTimeHit  = 0; 

    for(unsigned int i = 0; i < vec_HitTime.size(); i++)
    {
      for(unsigned int  j = 0; j < vec_HitTime.at(i).size(); j++)
      {
        if(vec_HitTime.at(i).at(j)<all_FirstTimeHit)
        {
          all_FirstTimeHit = vec_HitTime.at(i).at(j);
        }
        if(vec_HitTime.at(i).at(j)>all_LastTimeHit)
        {
          all_LastTimeHit = vec_HitTime.at(i).at(j);
        }
        if(vec_HitChan.at(i).at(j)<all_StartChan)
        {
          all_StartChan = vec_HitChan.at(i).at(j);
        }
        if(vec_HitChan.at(i).at(j)>all_EndChan)
        {
          all_EndChan = vec_HitChan.at(i).at(j);
        }
      }
      TBox *box = new TBox(vec_StartChan.at(i)-0.8, vec_FirstTimeHit.at(i)-25, vec_EndChan.at(i)+0.8, vec_LastTimeHit.at(i)+25);
      vec_Box.push_back(box);
    }
    TH2F *h_Cluster = new TH2F("h_Cluster_All", Form("Event %i with %lu Clusters", event, vec_HitTime.size()), 
        std::abs(all_EndChan-all_StartChan)+5,all_StartChan-2, all_EndChan+2, 
        std::abs(all_LastTimeHit-all_FirstTimeHit)+61, all_FirstTimeHit-30, all_LastTimeHit+30);
    h_Cluster->GetXaxis()->SetTitle("Channel, (Channel No)");
    h_Cluster->GetYaxis()->SetTitle("Time, (tick)");

    for(unsigned int i = 0; i < vec_HitTime.size(); i++)
    {
      for(unsigned int  j = 0; j < vec_HitTime.at(i).size(); j++)
      {
        int HitSpread      = std::ceil(2*vec_HitRMS.at(i).at(j)); 
        float StartSpread  = vec_HitTime.at(i).at(j)-(HitSpread/2.); 
        float EndSpread    = vec_HitTime.at(i).at(j)+(HitSpread/2.);
        float ChargeSpread = vec_HitSADC.at(i).at(j)/(float)HitSpread; 
        for(int k = 0; k < HitSpread; k++)
        {
          float TimePosition = StartSpread + k;
          h_Cluster->Fill(vec_HitChan.at(i).at(j), TimePosition, ChargeSpread);
        }
      }
    }

    gStyle->SetOptStat(0);
    TCanvas *c_Cluster = new TCanvas("c_Cluster", "c_Cluster", 800, 500);
    h_Cluster->SetContour(100);
    h_Cluster->Draw("COLZ");
    for(unsigned int i = 0; i < vec_Box.size(); i++)
    {
      vec_Box.at(i)->SetFillStyle(0);
      vec_Box.at(i)->SetLineColor(1);
      vec_Box.at(i)->Draw();
    }
    if(makePDF==true)
    {
      c_Cluster->SaveAs(Form("AllClusters_Event%i_AllHits.pdf", event));
    }
  }
  else
  {
    std::vector<TBox*> vec_Box;
    std::vector<std::vector<int>> marleyHits;
    int   all_StartChan = 10e5;
    int   all_EndChan   = 0;
    float all_FirstTimeHit = 5000;
    float all_LastTimeHit  = 0; 

    for(unsigned int i = 0; i < vec_HitTime.size(); i++)
    {
      int ClustersEdge_StartChan = 10e5;
      int ClustersEdge_EndChan = 0;
      float ClustersEdge_FirstTimeHit = 5000;
      float ClustersEdge_LastTimeHit = 0;
      std::vector<int> temp_marleyHits;
      for(unsigned int  j = 0; j < vec_HitTime.at(i).size(); j++)
      {
        if(vec_GenType.at(i).at(j)==1)
        {
          temp_marleyHits.push_back(j);
          if(vec_HitTime.at(i).at(j)<all_FirstTimeHit)
          {
            all_FirstTimeHit = vec_HitTime.at(i).at(j);
          }
          if(vec_HitTime.at(i).at(j)>all_LastTimeHit)
          {
            all_LastTimeHit = vec_HitTime.at(i).at(j);
          }
          if(vec_HitChan.at(i).at(j)<all_StartChan)
          {
            all_StartChan = vec_HitChan.at(i).at(j);
          }
          if(vec_HitChan.at(i).at(j)>all_EndChan)
          {
            all_EndChan = vec_HitChan.at(i).at(j);
          }
          if(vec_HitTime.at(i).at(j)<ClustersEdge_FirstTimeHit)
          {
            ClustersEdge_FirstTimeHit = vec_HitTime.at(i).at(j);
          }
          if(vec_HitTime.at(i).at(j)>ClustersEdge_LastTimeHit)
          {
            ClustersEdge_LastTimeHit = vec_HitTime.at(i).at(j);
          }
          if(vec_HitChan.at(i).at(j)<ClustersEdge_StartChan)
          {
            ClustersEdge_StartChan = vec_HitChan.at(i).at(j);
          }
          if(vec_HitChan.at(i).at(j)>ClustersEdge_EndChan)
          {
            ClustersEdge_EndChan = vec_HitChan.at(i).at(j);
          }
        }
      }
      TBox *box = new TBox(ClustersEdge_StartChan-0.8, ClustersEdge_FirstTimeHit-25, ClustersEdge_EndChan+0.8, ClustersEdge_LastTimeHit+25);
      vec_Box.push_back(box);
      marleyHits.push_back(temp_marleyHits);
    }
    TH2F *h_Cluster = new TH2F("h_Cluster_All", Form("Event %i with %lu Clusters",event, vec_HitTime.size()), 
        std::abs(all_EndChan-all_StartChan)+5,all_StartChan-2, all_EndChan+2, 
        std::abs(all_LastTimeHit-all_FirstTimeHit)+61, all_FirstTimeHit-30, all_LastTimeHit+30);
    h_Cluster->GetXaxis()->SetTitle("Channel, (Channel No)");
    h_Cluster->GetYaxis()->SetTitle("Time, (tick)");

    for(unsigned int i = 0; i < vec_HitTime.size(); i++)
    {
      for(unsigned int  j = 0; j < marleyHits.at(i).size(); j++)
      {
        int HitSpread      = std::ceil(2*vec_HitRMS.at(i).at(marleyHits.at(i).at(j))); 
        float StartSpread  = vec_HitTime.at(i).at(marleyHits.at(i).at(j))-(HitSpread/2.); 
        float EndSpread    = vec_HitTime.at(i).at(marleyHits.at(i).at(j))+(HitSpread/2.);
        float ChargeSpread = vec_HitSADC.at(i).at(marleyHits.at(i).at(j))/(float)HitSpread; 
        for(int k = 0; k < HitSpread; k++)
        {
          float TimePosition = StartSpread + k;
          h_Cluster->Fill(vec_HitChan.at(i).at(marleyHits.at(i).at(j)), TimePosition, ChargeSpread);
        }
      }
    }

    gStyle->SetOptStat(0);
    TCanvas *c_Cluster = new TCanvas("c_Cluster", "c_Cluster", 800, 500);
    h_Cluster->SetContour(100);
    h_Cluster->Draw("COLZ");
    for(unsigned int i = 0; i < vec_Box.size(); i++)
    {
      vec_Box.at(i)->SetFillStyle(0);
      vec_Box.at(i)->SetLineColor(1);
      vec_Box.at(i)->Draw();
    }
    if(makePDF==true)
    {
      c_Cluster->SaveAs(Form("AllClusters_Event%i_MarlHits.pdf", event));
    }
  }

  return;
}


void getClusterInformation(int cCluster)
{
  Cluster = cCluster;
  cut_Cluster = Form("Cluster==%i", Cluster);

  t->SetBranchAddress("Event", &Event);
  t->SetBranchAddress("NChan", &NChan);
  t->SetBranchAddress("NHits", &NHits);
  t->SetBranchAddress("Type", &Type);
  t->SetBranchAddress("StartChan", &StartChan);
  t->SetBranchAddress("EndChan", &EndChan);
  t->SetBranchAddress("FirstTimeHit", &FirstTimeHit);
  t->SetBranchAddress("LastTimeHit",  &LastTimeHit);
  t->SetBranchAddress("SumADC",  &SumADC);
  t->SetBranchAddress("HitView", &HitView);
  t->SetBranchAddress("GenType", &GenType);
  t->SetBranchAddress("HitChan", &HitChan);
  t->SetBranchAddress("HitTime", &HitTime);
  t->SetBranchAddress("HitSADC", &HitSADC);
  t->SetBranchAddress("HitRMS",  &HitRMS);

  t->GetEntry(cCluster);

  std::cout << "***********************************************************************************************************************" 
            << std::endl;
  std::cout << "CLUSTER " << cCluster  << " BELONGS TO EVENT: " << Event << ", CONTAINS " << NChan << " CHANNELS AND " << NHits
            << " TOTAL HITS. THE ADC SUM IS " << HitSADC << ". TYPE IS " << Type << ", THE HITS ARE:" << std::endl;
  for(int i = 0; i < NHits; i++)
  {
    std::cout << "Event: " << Event << " Plane: " << HitView->at(i) << " Generated by: " << GenType->at(i) 
              << " On channel: " << HitChan->at(i) << " At time: " << HitTime->at(i) << std::endl; 
  }
  std::cout << "***********************************************************************************************************************" 
            << std::endl;
}


void drawAllMarleyHits(TString s_File, TString s_Path, long long cEvent, bool makePDF)
{
  std::vector<int> v_HitChan;
  std::vector<float> v_HitTime;
  std::vector<float> v_HitRMS;
  std::vector<float> v_HitSADC;
  std::cout << "DRAWING ALL MARLEY HITS, EVENT: " << cEvent << std::endl;

  TFile *f_Original = new TFile(s_File, "READ");
  TTree *t_Original = (TTree*)f_Original->Get(s_Path);

  int Or_NColHits;
  int Or_HitChan[50000];
  int Or_GenType[50000];
  float Or_HitRMS[50000];
  float Or_HitTime[50000];
  float Or_HitSADC[50000];

  t_Original->SetBranchAddress("NColHits", &Or_NColHits);
  t_Original->SetBranchAddress("HitChan",  &Or_HitChan);
  t_Original->SetBranchAddress("GenType",  &Or_GenType);
  t_Original->SetBranchAddress("HitRMS",   &Or_HitRMS);
  t_Original->SetBranchAddress("HitTime",  &Or_HitTime);
  t_Original->SetBranchAddress("HitSADC",  &Or_HitSADC);

  int   Or_StartChan    = 10e5;
  int   Or_EndChan      = 0;
  float Or_FirstTimeHit = 5000;
  float Or_LastTimeHit  = 0; 

  t_Original->GetEntry(cEvent);

  for(int j = 0; j < Or_NColHits; j++)
  {
    if(Or_GenType[j]==1)
    {
      v_HitChan.push_back(Or_HitChan[j]);
      v_HitRMS.push_back(Or_HitRMS[j]);
      v_HitTime.push_back(Or_HitTime[j]);
      v_HitSADC.push_back(Or_HitSADC[j]);
      if(Or_HitChan[j]<Or_StartChan)
      {
        Or_StartChan = Or_HitChan[j];
      }
      if(Or_HitChan[j]>Or_EndChan)
      {
        Or_EndChan = Or_HitChan[j];
      }
      if(Or_HitTime[j]<Or_FirstTimeHit)
      {
        Or_FirstTimeHit = Or_HitTime[j];
      }
      if(Or_HitTime[j]>Or_LastTimeHit)
      {
        Or_LastTimeHit = Or_HitTime[j];
      }
    }
  }

  TH2F *h_Marley = new TH2F("h_Marley", Form("All Marley Hits, Event %lli", cEvent), 
                             std::abs(Or_EndChan-Or_StartChan)+5,Or_StartChan-2, Or_EndChan+2, 
                             std::abs(Or_LastTimeHit-Or_FirstTimeHit)+61, Or_FirstTimeHit-30, Or_LastTimeHit+30);
  h_Marley->GetXaxis()->SetTitle("Channel, (Channel No)");
  h_Marley->GetYaxis()->SetTitle("Time, (tick)");

  for(long long i = 0; i < v_HitChan.size(); i++)
  {
    int HitSpread      = std::ceil(2*v_HitRMS.at(i)); 
    float StartSpread  = v_HitTime.at(i)-(HitSpread/2.); 
    float EndSpread    = v_HitTime.at(i)+(HitSpread/2.);
    float ChargeSpread = v_HitSADC.at(i)/(float)HitSpread; 
    for(int j = 0; j < HitSpread; j++)
    {
      float TimePosition = StartSpread + j;
      h_Marley->Fill(v_HitChan.at(i), TimePosition, ChargeSpread);
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *c_Marley = new TCanvas("c_Marley", "c_Marley", 800, 500);
  h_Marley->SetContour(100);
  h_Marley->Draw("COLZ");
  if(makePDF==true)
  {
    c_Marley->SaveAs(Form("AllMarleyHits_Event%lli.pdf", cEvent));
  }

  return;
}

void getEventInformation(int cEvent, int cConfig)
{
  event  = cEvent;
  config = cConfig;
  cut_Event  = Form("Event==%i", event);
  cut_Config = Form("Config==%i", config);

  t->SetBranchAddress("Config", &Config);
  t->SetBranchAddress("Cluster", &Cluster);
  t->SetBranchAddress("Event", &Event);
  t->SetBranchAddress("NChan", &NChan);
  t->SetBranchAddress("NHits", &NHits);
  t->SetBranchAddress("Type", &Type);
  t->SetBranchAddress("StartChan", &StartChan);
  t->SetBranchAddress("EndChan", &EndChan);
  t->SetBranchAddress("FirstTimeHit", &FirstTimeHit);
  t->SetBranchAddress("LastTimeHit",  &LastTimeHit);
  t->SetBranchAddress("SumADC", &SumADC);
  t->SetBranchAddress("HitView", &HitView);
  t->SetBranchAddress("GenType", &GenType);
  t->SetBranchAddress("HitChan", &HitChan);
  t->SetBranchAddress("HitTime", &HitTime);
  t->SetBranchAddress("HitSADC", &HitSADC);
  t->SetBranchAddress("HitRMS", &HitRMS);

  std::vector<int>   vec_temp_hit_Event;
  std::vector<int>   vec_temp_HitView;
  std::vector<int>   vec_temp_GenType;
  std::vector<int>   vec_temp_HitChan;
  std::vector<float> vec_temp_HitTime;
  std::vector<float> vec_temp_HitSADC;
  std::vector<float> vec_temp_HitRMS;

  for(int i = 0; i < t->GetEntries(); i++)
  {
    vec_temp_hit_Event.clear();
    vec_temp_HitView.clear();
    vec_temp_GenType.clear();
    vec_temp_HitChan.clear();
    vec_temp_HitTime.clear();
    vec_temp_HitSADC.clear();
    vec_temp_HitRMS.clear();

    t->GetEntry(i);
    if(Event==cEvent && Config==cConfig)
    {
      vec_Event.push_back(Event);
      vec_NChan.push_back(NChan);
      vec_NHits.push_back(NHits);
      vec_Type.push_back(Type);
      vec_StartChan.push_back(StartChan);
      vec_EndChan.push_back(EndChan);
      vec_FirstTimeHit.push_back(FirstTimeHit);
      vec_LastTimeHit.push_back(LastTimeHit);
      vec_SumADC.push_back(SumADC);

      for(int j = 0; j < NHits; j++)
      {
        vec_temp_HitView.push_back(HitView->at(j));
        vec_temp_GenType.push_back(GenType->at(j));
        vec_temp_HitChan.push_back(HitChan->at(j));
        vec_temp_HitTime.push_back(HitTime->at(j));
        vec_temp_HitSADC.push_back(HitSADC->at(j));
        vec_temp_HitRMS.push_back(HitRMS->at(j));
      }
      vec_HitView.push_back(vec_temp_HitView);
      vec_GenType.push_back(vec_temp_GenType);
      vec_HitChan.push_back(vec_temp_HitChan);
      vec_HitTime.push_back(vec_temp_HitTime);
      vec_HitSADC.push_back(vec_temp_HitSADC);
      vec_HitRMS.push_back(vec_temp_HitRMS);

      std::cout << "***********************************************************************************************************************" 
        << std::endl;
      std::cout << "CLUSTER " << Cluster  << " BELONGS TO EVENT: " << Event << ", CONTAINS " << NChan << " CHANNELS AND " << NHits
        << " TOTAL HITS. THE ADC SUM IS " << SumADC << ". TYPE IS " << Type << ", THE HITS ARE:" << std::endl;
      for(int j = 0; j < NHits; j++)
      {
        std::cout << "Event: " << Event << " Plane: " << HitView->at(j) << " Generated by: " << GenType->at(j) 
          << " On channel: " << HitChan->at(j) << " At time: " << HitTime->at(j) << std::endl; 
      }
      std::cout << "***********************************************************************************************************************" 
        << std::endl;
    }
  }

  return;
}


void loadFile(TString s_File)
{
  TFile *f_Input = new TFile(s_File, "READ");

  t = (TTree*)f_Input->Get("t_Output");

  std::cout << "TTREE LOADED. THERE ARE " << t->GetEntries() << " CLUSTERS." << std::endl;

  return;
}


bool clearCluster()
{
  Cluster = 0;

  Event = 0;
  NChan = 0;
  NHits = 0;
  Type = 0;
  StartChan = 0;
  EndChan = 0;
  FirstTimeHit = 0;
  LastTimeHit = 0;
  SumADC = 0;

  HitView->clear();
  GenType->clear();
  HitChan->clear();
  HitTime->clear();
  HitSADC->clear();
  HitRMS->clear();

  vec_Event.clear();
  vec_NChan.clear();
  vec_NHits.clear();
  vec_Type.clear();
  vec_StartChan.clear();
  vec_EndChan.clear();
  vec_FirstTimeHit.clear();
  vec_LastTimeHit.clear();
  vec_SumADC.clear();

  vec_HitView.clear();
  vec_GenType.clear();
  vec_HitChan.clear();
  vec_HitTime.clear();
  vec_HitSADC.clear();
  vec_HitRMS.clear();

  return true;
}


void EventDisplay()
{
  std::cout << "***********************************************************************************************************************" 
            << std::endl;
  std::cout << std::endl;
  std::cout <<"                                          THIS IS THE CLUSTER EVENT DISPLAY                                             "
            << std::endl;
  std::cout <<"                                             Tuesday January 23rd, 2018                                                 "
            << std::endl;
  std::cout <<"                                                 Written By A. Booth                                                    "
            << std::endl;
  std::cout << std::endl;
  std::cout << "***********************************************************************************************************************" 
            << std::endl;
  std::cout << "TO LOOK AT INDIVIDUAL CLUSTERS" << std::endl;
  std::cout << "Step 1: Load the .root file -> loadFile(\"RootFile.root\")" << std::endl;
  std::cout << "Step 2: View cluster information and pull that information from the tree -> getClusterInformation(int ClusterNumber)" 
            << std::endl;
  std::cout << "Step 3: Analyse the cluster -> drawSingleCluster(bool MarleyOnly, bool MakePDF)" << std::endl;
  std::cout << "Step 4: Before moving on to next cluster, clear the previous cluster's information -> clearCluster()" << std::endl;
  std::cout << "Step 5: Go from Step 2." << std::endl;
  std::cout << std::endl;
  std::cout << "TO LOOK AT A WHOLE EVENT" << std::endl;
  std::cout << "Step 1: Load the .root file -> loadFile(\"RootFile.root\")" << std::endl;
  std::cout << "Step 2: View the event information and pull that information from the tree -> getEventInformation(int EventNumber, int Config)" 
            << std::endl;
  std::cout << "Step 3: Analyse the cluster -> drawAllClusters(bool MarleyOnly, bool MakePDF)" << std::endl;
  std::cout << "Step 4: Before moving on to next cluster, clear the previous cluster's information -> clearCluster()" << std::endl;
  std::cout << "Step 5: Go from Step 2." << std::endl;
  std::cout << std::endl;
  std::cout << "TO COMPARE TO ORIGINAL MARLEY HITS IN SNANA" << std::endl;
  std::cout << "Step 1: Can draw straight away -> drawAllMarleyHits(\"DAQSimFile.root\", \"daqsimtree/path\", int Event, bool makePDF)"
            << std::endl;
  std::cout << std::endl;
  return;
}

int main()
{
  return 0;
}
