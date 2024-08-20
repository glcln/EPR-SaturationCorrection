#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>

double langaufun(double *x, double *par)
{
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

void CompareLaugausGaussian()
{
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    TF1 *landau = new TF1("landau", "landau", -5, 10);
    landau->SetParameters(1, +0.22278298, 1); // width, MPV, integral
    TF1 *langaus = new TF1("langaus", langaufun, -5, 10, 4);
    langaus->SetParameters(1, -0.22278298, 1, 1); // landau_width, MPV, integral, sigma_Gauss

    landau->SetLineColor(kRed);
    langaus->SetLineColor(kBlue);

    landau->Draw();
    langaus->Draw("same");
    gStyle->SetOptStat(0);

    TLegend *leg = new TLegend(0.55, 0.7, 0.85, 0.85);
    leg->AddEntry(landau, "Landau", "l");
    leg->AddEntry(langaus, "Landau convoluted with a gaussian", "l");
    leg->SetLineColor(0);
    leg->Draw();
}
