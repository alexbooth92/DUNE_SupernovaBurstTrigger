#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TMath.h>


int main()
{
  TFile *f_Input = new TFile("/Users/alexanderbooth/Documents/Work/Year1/SNTrigger/Samples/GH_SNMC.root","READ");
  TFile *f_Output = new TFile("TimeSample.root", "RECREATE");

  TTree *t = (TTree*)f_Input->Get("DAQSimTree");
  TH1D *h_MarlTime = new TH1D("h_MarlTime","h_MarlTime", 10000, -1, 11);
  TH1D *h_MarlTime_Zero20ms = new TH1D("h_MarlTime_Zero20ms","h_MarlTime_Zero20ms", 10000, -1, 11);

  double extrapToTime = 120;
  int    nBins_Extrap = std::ceil((extrapToTime+1)/h_MarlTime_Zero20ms->GetBinWidth(1));
  TH1D *h_MarlTime_Extrap = new TH1D("h_MarlTime_Extrap","h_MarlTime_Extrap", nBins_Extrap, -1, extrapToTime);
  TH1D *h_MarlTime_Zero20ms_Extrap = new TH1D("h_MarlTime_Zero20ms_Extrap","h_MarlTime_Zero20ms_Extrap", nBins_Extrap, -1, extrapToTime);
  TH1D *h_MarlTime_Extrap_Fixed = new TH1D("h_MarlTime_Extrap_3secs","h_MarlTime_Extrap_3secs", nBins_Extrap, -1, extrapToTime);
  TH1D *h_MarlTime_Zero20ms_Extrap_Fixed = new TH1D("h_MarlTime_Zero20ms_Extrap_3secs",
                                                    "h_MarlTime_Zero20ms_Extrap_3secs", nBins_Extrap, -1, extrapToTime);

  double endSuppressionTime = 20e-3;
  double MarlTime;
  t->SetBranchAddress("MarlTime", &MarlTime);

  for(unsigned int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    h_MarlTime->Fill(MarlTime);
    h_MarlTime_Extrap->Fill(MarlTime);
    h_MarlTime_Extrap_Fixed->Fill(MarlTime);
    if(MarlTime<=endSuppressionTime)
    {
      h_MarlTime_Zero20ms->Fill(MarlTime,0);
      h_MarlTime_Zero20ms_Extrap->Fill(MarlTime,0);
      h_MarlTime_Zero20ms_Extrap_Fixed->Fill(MarlTime,0);
    }
    else
    {
      h_MarlTime_Zero20ms->Fill(MarlTime);
      h_MarlTime_Zero20ms_Extrap->Fill(MarlTime);
      h_MarlTime_Zero20ms_Extrap_Fixed->Fill(MarlTime);
    }
  }
  //GET GENERAL INFORMATION ABOUT THE TIME PROFILE.
  int    sampleEvery   = 1;
  int    startBin      = 0;
  double startTime     = 0;
  int    endBin        = 0;
  double endTime       = 0;
  double snDuration    = 0;
  double secondsPerBin = 0;
  double sensitivity   = 0;
  for(int i = 1; i < h_MarlTime->GetSize()-1; i++)
  {
    if(h_MarlTime->GetBinContent(i)>0)
    {
      startBin  = i;
      startTime = h_MarlTime->GetBinCenter(i);
      break;
    }
  }
  for(int i = h_MarlTime->GetSize()-2; i > 0; i--)
  {
    if(h_MarlTime->GetBinContent(i)>0)
    {
      endBin  = i;
      endTime = h_MarlTime->GetBinCenter(i);
      break;
    }
  }
  snDuration    = endTime - startTime;

  std::cout << "THE SNs START AT TIME t = " << startTime  << ", BIN " << startBin << std::endl;
  std::cout << "THE SNs END AT TIME t = "   << endTime    << ", BIN " << endBin << std::endl;
  std::cout << "THE SNs DURATION IS "       << snDuration << "s"      << std::endl;

  //FIT AN EXPONENTIAL TO THE COOLING PART OF THE SPECTRUM THEN EXTRAPOLATED TO 30SEC.
  TF1 *f_Cooling_Fit = new TF1("f_Cooling_Fit", "expo",3,10);
  f_Cooling_Fit->SetLineColor(1);
  h_MarlTime->Fit("f_Cooling_Fit","R");
  double constant = f_Cooling_Fit->GetParameter(0);
  double slope    = f_Cooling_Fit->GetParameter(1);
  TF1 *f_Cooling_Extrap = new TF1("f_Cooling_Extrap", "TMath::Exp([0]+[1]*x)", endTime, extrapToTime);
  f_Cooling_Extrap->SetParameter(0,constant);
  f_Cooling_Extrap->SetParameter(1,slope);

  //ADD AN EXPONONTIAL WITH A DECAY TIME OF 3 SECONDS. OBTAIN THE CONSTANT BY EVALUATING THE FIT ABOVE AT 10SECS.
  double slopeFixed    = -1./3.;
  double constantFixed = TMath::Log(f_Cooling_Extrap->Eval(10)) - slopeFixed*10; 
  TF1 *f_Cooling_Extrap_Fixed = new TF1("f_Cooling_Extrap_3secs", "TMath::Exp([0]+[1]*x)", endTime, extrapToTime);
  f_Cooling_Extrap_Fixed->SetParameter(0,constantFixed);
  f_Cooling_Extrap_Fixed->SetParameter(1,slopeFixed);

  for(unsigned int i = endBin; i < h_MarlTime_Extrap->GetSize()-1; i++)
  {
    double binCenter = h_MarlTime_Extrap->GetBinCenter(i);
    h_MarlTime_Extrap->SetBinContent(i,f_Cooling_Extrap->Eval(binCenter));
    h_MarlTime_Zero20ms_Extrap->SetBinContent(i,f_Cooling_Extrap->Eval(binCenter));
    h_MarlTime_Extrap_Fixed->SetBinContent(i,f_Cooling_Extrap_Fixed->Eval(binCenter));
    h_MarlTime_Zero20ms_Extrap_Fixed->SetBinContent(i,f_Cooling_Extrap_Fixed->Eval(binCenter));
  }

  TCanvas *c_MarlTime_Extrap = new TCanvas("c_MarlTime_Extrap", "c_MarlTime_Extrap", 800, 500);
  h_MarlTime_Extrap->Draw();
  f_Cooling_Fit->Draw("SAME");
  f_Cooling_Extrap->Draw("SAME");
  f_Cooling_Extrap_Fixed->Draw("SAME");
  c_MarlTime_Extrap->Write();

  h_MarlTime->Scale(1/(double)h_MarlTime->GetEntries());
  h_MarlTime_Zero20ms->Scale(1/(double)h_MarlTime_Zero20ms->GetEntries());

  h_MarlTime->Write();
  h_MarlTime_Zero20ms->Write();

  double fracEventsInFirst20ms = h_MarlTime_Zero20ms_Extrap->Integral(1,-1)/h_MarlTime_Extrap->Integral(1,-1);
  h_MarlTime_Extrap->Scale(1/h_MarlTime_Extrap->Integral(1,-1));
  h_MarlTime_Extrap_Fixed->Scale(1/h_MarlTime_Extrap_Fixed->Integral(1,-1));
  h_MarlTime_Zero20ms_Extrap->Scale(fracEventsInFirst20ms/h_MarlTime_Zero20ms_Extrap->Integral(1,-1));
  h_MarlTime_Zero20ms_Extrap_Fixed->Scale(fracEventsInFirst20ms/h_MarlTime_Zero20ms_Extrap_Fixed->Integral(1,-1));
  h_MarlTime_Extrap->Write();
  h_MarlTime_Extrap_Fixed->Write();
  h_MarlTime_Zero20ms_Extrap->Write();
  h_MarlTime_Zero20ms_Extrap_Fixed->Write();


  std::cout << "INTEGRALS:" << std::endl; 
  std::cout << "MARLTIME: "        << h_MarlTime->Integral(1,-1)        << ", MARLTIME ZERO 20MS: " << h_MarlTime_Zero20ms->Integral(1,-1) << std::endl;
  std::cout << "MARLTIME EXTRAP: " << h_MarlTime_Extrap->Integral(1,-1) << ", MARLTIME ZERO 20MS EXTRAP: " 
            << h_MarlTime_Zero20ms_Extrap->Integral(1,-1) << std::endl;
  std::cout << "MARLTIME EXTRAP FIXED: " << h_MarlTime_Extrap_Fixed->Integral(1,-1) << ", MARLTIME ZERO 20MS EXTRAP FIXED: " 
            << h_MarlTime_Zero20ms_Extrap_Fixed->Integral(1,-1) << std::endl;
  return 0;
}
