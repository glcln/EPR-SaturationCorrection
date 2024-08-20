#include <iostream>
#include <TFile.h>
#include <TH2F.h>

void Return_K_C_fit(double &K, double &C, std::string filename)
{
  std::ifstream infile(filename);
  infile >> K >> C;
  infile.close();

  return;
}

TF1 *Fit_MassParametrization(double mass, double p_start, double p_end, double C, double K)
{
  TF1 *fit = new TF1("fit", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit->SetParameter(0, K);
  fit->FixParameter(1, C);

  return fit;
}

void Display_TH2_Fit(TH2F *h2, TF1 *fit_pion, TF1 *fit_kaon, TF1 *fit_proton, TF1 *fit_deuteron, std::string filename_PDF, std::string filename_ROOT)
{
  gROOT->SetBatch(kTRUE);

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Track momentum [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, 14);

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h2->Draw("colz");
  c->SetLogz();
  fit_pion->Draw("same");
  fit_kaon->Draw("same");
  fit_proton->Draw("same");
  fit_deuteron->Draw("same");

  c->SaveAs(filename_PDF.c_str());
  c->SaveAs(filename_ROOT.c_str());

  delete c;

  return;
}


void Display_fit_wjet()
{
    // Initialisation
    TFile *file = new TFile("ROOT_SVG/Wjets_checkCorrection.root", "READ");

    // Récupérer le TH2F
    TH2F *dEdX0stripVsP_nocorr = (TH2F*)file->Get("Ih_VS_p_SatNoCorr");
    TH2F *dEdX0stripVsP_oldcorr = (TH2F*)file->Get("Ih_VS_p_OldCorr");
    TH2F *dEdX0stripVsP_newcorr = (TH2F*)file->Get("Ih_VS_p_SatCorr");
    
    // Récupérer le Fit
    double K_nocorr = -1, K_oldcorr = -1, K_newcorr = -1;
    double C_nocorr = -1, C_oldcorr = -1, C_newcorr = -1;
    Return_K_C_fit(K_nocorr, C_nocorr, "ROOT_SVG/K_C_fit_Wjets/K_C_fit_C_file_nocorr.txt");
    Return_K_C_fit(K_oldcorr, C_oldcorr, "ROOT_SVG/K_C_fit_Wjets/K_C_fit_C_file_oldcorr.txt");
    Return_K_C_fit(K_newcorr, C_newcorr, "ROOT_SVG/K_C_fit_Wjets/K_C_fit_C_file_newcorr.txt");

    TF1* Fit_pion_nocorr = Fit_MassParametrization(0.13957, 0.5, 10, C_nocorr, K_nocorr);
    TF1* Fit_kaon_nocorr = Fit_MassParametrization(0.49368, 0.5, 10, C_nocorr, K_nocorr);
    TF1* Fit_proton_nocorr = Fit_MassParametrization(0.93827, 0.5, 10, C_nocorr, K_nocorr);
    TF1* Fit_deuteron_nocorr = Fit_MassParametrization(1.87561, 0.5, 10, C_nocorr, K_nocorr);

    TF1* Fit_pion_oldcorr = Fit_MassParametrization(0.13957, 0.5, 10, C_oldcorr, K_oldcorr);
    TF1* Fit_kaon_oldcorr = Fit_MassParametrization(0.49368, 0.5, 10, C_oldcorr, K_oldcorr);
    TF1* Fit_proton_oldcorr = Fit_MassParametrization(0.93827, 0.5, 10, C_oldcorr, K_oldcorr);
    TF1* Fit_deuteron_oldcorr = Fit_MassParametrization(1.87561, 0.5, 10, C_oldcorr, K_oldcorr);

    TF1* Fit_pion_newcorr = Fit_MassParametrization(0.13957, 0.5, 10, C_newcorr, K_newcorr);
    TF1* Fit_kaon_newcorr = Fit_MassParametrization(0.49368, 0.5, 10, C_newcorr, K_newcorr);
    TF1* Fit_proton_newcorr = Fit_MassParametrization(0.93827, 0.5, 10, C_newcorr, K_newcorr);
    TF1* Fit_deuteron_newcorr = Fit_MassParametrization(1.87561, 0.5, 10, C_newcorr, K_newcorr);

    // Affichage
    Display_TH2_Fit(dEdX0stripVsP_nocorr, Fit_pion_nocorr, Fit_kaon_nocorr, Fit_proton_nocorr, Fit_deuteron_nocorr, "ROOT_SVG/K_C_fit_Wjets/Fit_nocorr.pdf", "ROOT_SVG/K_C_fit_Wjets/Fit_nocorr.root");
    Display_TH2_Fit(dEdX0stripVsP_oldcorr, Fit_pion_oldcorr, Fit_kaon_oldcorr, Fit_proton_oldcorr, Fit_deuteron_oldcorr, "ROOT_SVG/K_C_fit_Wjets/Fit_oldcorr.pdf", "ROOT_SVG/K_C_fit_Wjets/Fit_oldcorr.root");
    Display_TH2_Fit(dEdX0stripVsP_newcorr, Fit_pion_newcorr, Fit_kaon_newcorr, Fit_proton_newcorr, Fit_deuteron_newcorr, "ROOT_SVG/K_C_fit_Wjets/Fit_newcorr.pdf", "ROOT_SVG/K_C_fit_Wjets/Fit_newcorr.root");

}
