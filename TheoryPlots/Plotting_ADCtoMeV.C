#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>

int main()
{
  TFile *f_Input = new TFile("ADCtoMeV.root", "READ");

  TH2D *h_RatioVELep_D1 = (TH2D*)f_Input->Get("h_RatioVELep_D1");
  TH2D *h_RatioVELep_D2 = (TH2D*)f_Input->Get("h_RatioVELep_D2");
  TH2D *h_RatioVELep_D3 = (TH2D*)f_Input->Get("h_RatioVELep_D3");

  TH2D *h_RatioVELep_A1 = (TH2D*)f_Input->Get("h_RatioVELep_A1");
  TH2D *h_RatioVELep_A2 = (TH2D*)f_Input->Get("h_RatioVELep_A2");
  TH2D *h_RatioVELep_A3 = (TH2D*)f_Input->Get("h_RatioVELep_A3");

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c", 800, 500);
  h_RatioVELep_D1->SetContour(1000);
  h_RatioVELep_D1->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. D1");
  h_RatioVELep_D1->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_D1->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_D1->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_D1->Draw("COLZ");
  c->SaveAs("RatioVELep_D1.pdf");
  c->Clear();

  h_RatioVELep_D2->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. D2");
  h_RatioVELep_D2->SetContour(1000);
  h_RatioVELep_D2->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_D2->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_D2->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_D2->Draw("COLZ");
  c->SaveAs("RatioVELep_D2.pdf");
  c->Clear();

  h_RatioVELep_D3->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. D3");
  h_RatioVELep_D3->SetContour(1000);
  h_RatioVELep_D3->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_D3->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_D3->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_D3->Draw("COLZ");
  c->SaveAs("RatioVELep_D3.pdf");
  c->Clear();

  h_RatioVELep_A1->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. A1");
  h_RatioVELep_A1->SetContour(1000);
  h_RatioVELep_A1->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_A1->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_A1->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_A1->Draw("COLZ");
  c->SaveAs("RatioVELep_A1.pdf");
  c->Clear();

  h_RatioVELep_A2->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. A2");
  h_RatioVELep_A2->SetContour(1000);
  h_RatioVELep_A2->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_A2->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_A2->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_A2->Draw("COLZ");
  c->SaveAs("RatioVELep_A2.pdf");
  c->Clear();

  h_RatioVELep_A3->SetTitle("Total Marley Cluster ADC Sum/True Primary Lepton Energy vs. True Primary Lepton Energy. A3");
  h_RatioVELep_A3->SetContour(1000);
  h_RatioVELep_A3->GetXaxis()->SetTitle("True Primary Lepton Energy, (MeV)");
  h_RatioVELep_A3->GetYaxis()->SetTitle("Cluster Energy/True Primary Lepton Energy, (ADC/MeV)");
  h_RatioVELep_A3->GetYaxis()->SetTitleOffset(1.2);
  h_RatioVELep_A3->Draw("COLZ");
  c->SaveAs("RatioVELep_A3.pdf");
  c->Clear();

  return 0;
}
