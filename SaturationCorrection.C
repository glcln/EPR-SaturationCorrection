#define SaturationCorrection_cxx
#include "SaturationCorrection.h"
#include <TH2.h>
#include <TMath.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace TMath;
using namespace std;

void DisplayCluster(const std::vector <int> C, TFile* SaveData, int entry, int itr, int layer, double VmaxOverVmin, double MaxOverVmax, bool print)
{
  TH1F *hist=new TH1F(Form("Cluster_%d_%d_%d__%+4.3f-%+4.3f",entry,itr,layer,VmaxOverVmin,MaxOverVmax),Form("Cluster_%d_%d_%d__%+4.3f-%+4.3f",entry,itr,layer,VmaxOverVmin,MaxOverVmax),C.size()+2,0,C.size()+2);
  hist->SetLineColor(kBlue);
  hist->SetLineWidth(2);
  hist->SetFillColor(kBlue);
  hist->SetFillStyle(3001);
  hist->GetXaxis()->SetTitle("Strip");
  hist->GetYaxis()->SetTitle("ADC");

  for (int i=0; i<C.size(); i++)
  {
    hist->SetBinContent(i+2,C[i]);
  }

  if (print)
  {
    SaveData->cd();
    hist->Write();
    hist->Delete();
    cout << "Writing in " << SaveData->GetName() << " : " << Form("Cluster_%d_%d_%d__%4.3f-%4.3f",entry,itr,layer,VmaxOverVmin,MaxOverVmax) << endl;
  }
}

void DisplayCluster(const std::vector <int> C1, const std::vector <int> C2, const std::vector <int> C3, TFile* SaveData, int shift, int entry, int itr, int layer, double local_theta, bool print)
{
  if (local_theta>Pi()/2 && local_theta<=Pi()) local_theta-=Pi();

  TH1F *hist_reco=new TH1F("hist_reco", "hist_reco", C1.size()+4, 0, C1.size()+4);
  TH1F *hist_sum=new TH1F("hist_sum", "hist_sum", C3.size()+4, 0, C3.size()+4);
  TH1F *hist_sig=new TH1F(Form("Cluster_%d_%d_%d_%+4.3f",entry,itr,layer,local_theta), Form("Cluster_%d_%d_%d_%+4.3f",entry,itr,layer,local_theta), C2.size()+4, 0, C2.size()+4);


  hist_reco->SetLineColor(kGreen);
  hist_reco->SetLineWidth(2);
  hist_reco->SetFillColor(kGreen);
  hist_reco->SetFillStyle(3001);
  hist_reco->GetXaxis()->SetTitle("Strip");
  hist_reco->GetYaxis()->SetTitle("ADC");

  hist_sig->SetLineColor(kBlue);
  hist_sig->SetLineWidth(2);
  hist_sig->GetXaxis()->SetTitle("Strip");
  hist_sig->GetYaxis()->SetTitle("ADC");

  hist_sum->SetLineColor(kRed);
  hist_sum->SetLineWidth(2);
  hist_sum->GetXaxis()->SetTitle("Strip");
  hist_sum->GetYaxis()->SetTitle("ADC");

  for (int i=0; i<C1.size(); i++) hist_reco->SetBinContent(i+2+shift,C1[i]);
  for (int i=0; i<C2.size(); i++) hist_sig->SetBinContent(i+2,C2[i]);
  for (int i=0; i<C3.size(); i++) hist_sum->SetBinContent(i+2+shift,C3[i]);

  if (print)
  {
    TCanvas *c1 = new TCanvas(Form("Cluster_%d_%d_%d_%4.3f",entry,itr,layer,local_theta),Form("Cluster_%d_%d_%d_%4.3f",entry,itr,layer,local_theta),600,600);
    hist_sig->Draw();
    hist_reco->Draw("same");
    hist_sum->Draw("same");
    TLegend *leg = new TLegend(0.15,0.7,0.3,0.8);
    leg->AddEntry(hist_sig,"SigDigi","l");
    leg->AddEntry(hist_sum,"SumDigi","l");
    leg->AddEntry(hist_reco,"SigReco","f");
    leg->SetBorderSize(0);
    leg->Draw();
    
    SaveData->cd();
    c1->Write();
    
    c1->Close();
    delete c1;
    cout << "Writing in " << SaveData->GetName() << " : " << Form("Cluster_%d_%d_%d_%4.3f",entry,itr,layer,local_theta) << endl;
  }
  
  delete hist_sig;
  delete hist_sum;
  delete hist_reco;
}

void DisplayCluster(const std::vector <int> C1, const std::vector <int> C2, TFile* SaveData, int shift, int layer, double local_theta, bool print)
{
  if (local_theta>Pi()/2 && local_theta<=Pi()) local_theta-=Pi();

  TH1F *hist_reco=new TH1F("hist_reco", "hist_reco", C1.size()+4, 0, C1.size()+4);
  TH1F *hist_sig=new TH1F(Form("Cluster_%d_%+4.3f",layer,local_theta), Form("Cluster_%d_%+4.3f",layer,local_theta), C2.size()+4, 0, C2.size()+4);

  hist_reco->SetLineColor(kGreen);
  hist_reco->SetLineWidth(2);
  hist_reco->SetFillColor(kGreen);
  hist_reco->SetFillStyle(3001);
  hist_reco->GetXaxis()->SetTitle("Strip");
  hist_reco->GetYaxis()->SetTitle("ADC");

  hist_sig->SetLineColor(kBlue);
  hist_sig->SetLineWidth(2);
  hist_sig->GetXaxis()->SetTitle("Strip");
  hist_sig->GetYaxis()->SetTitle("ADC");

  for (int i=0; i<C1.size(); i++) hist_reco->SetBinContent(i+2+shift,C1[i]);
  for (int i=0; i<C2.size(); i++) hist_sig->SetBinContent(i+2,C2[i]);

  if (print)
  {
    TCanvas *c1 = new TCanvas(Form("Cluster_%d_%+4.3f",layer,local_theta),Form("Cluster_%d_%+4.3f",layer,local_theta),600,600);
    hist_sig->Draw();
    hist_reco->Draw("same");
    hist_sig->Draw("same");
    TLegend *leg = new TLegend(0.15,0.7,0.3,0.8);
    leg->AddEntry(hist_sig,"True sim","l");
    leg->AddEntry(hist_reco,"Reco","f");
    leg->SetBorderSize(0);
    leg->Draw();
    
    SaveData->cd();
    c1->Write();
    
    c1->Close();
    delete c1;
    cout << "Writing in " << SaveData->GetName() << " : " << Form("Cluster_%d_%+4.3f",layer,local_theta) << endl;
  }
  
  delete hist_sig;
  delete hist_reco;
}

Double_t modulo2Pi(Double_t nb)
{
   if (nb>0){
      while(nb>Pi()) nb-=TwoPi();
   }
   else{
      while(nb<-Pi()) nb+=TwoPi();
   }
   return nb;
}

bool compareVectors(const std::vector<int>& vec1, const std::vector<int>& vec2)
{
  if (vec1.size() != vec2.size()) return false; 

  for (size_t i = 0; i < vec1.size(); ++i) 
  {
    if (vec1[i] != vec2[i]) return false; 
  }

  return true;
}

int FindLayer(int subdetid, int detid)
{
  if(subdetid==3)  // TIB
  {
    if(((detid>>14)&0x7)==1) return 1;
    else if(((detid>>14)&0x7)==2) return 2;
    else if(((detid>>14)&0x7)==3) return 3;
    else if(((detid>>14)&0x7)==4) return 4;
  }    
  else if(subdetid==5) // TOB
  {       
    if(((detid>>14)&0x7)==1) return 5;
    else if(((detid>>14)&0x7)==2) return 6;
    else if(((detid>>14)&0x7)==3) return 7;
    else if(((detid>>14)&0x7)==4) return 8;
    else if(((detid>>14)&0x7)==5) return 9;
    else if(((detid>>14)&0x7)==6) return 10;
  }       
  /*else if(subdetid==4)  //TID by disk 
  {    
    if(((detid>>11)&0x3)==1) return 11;
    else if(((detid>>11)&0x3)==2) return 12;
    else if(((detid>>11)&0x3)==3) return 13;
  }*/
  else if(subdetid==4)  //TID by ring 
  {    
    if(((detid>>9)&0x3)==1) return 11;
    else if(((detid>>9)&0x3)==2) return 12;
    else if(((detid>>9)&0x3)==3) return 13;
  }
  /*else if(subdetid==6) // TEC by wheel 
  {
    if(((detid>>14)&0xF)==1) return 14;
    else if(((detid>>14)&0xF)==2) return 15;
    else if(((detid>>14)&0xF)==3) return 16;
    else if(((detid>>14)&0xF)==4) return 17;
    else if(((detid>>14)&0xF)==5) return 18;
    else if(((detid>>14)&0xF)==6) return 19;
    else if(((detid>>14)&0xF)==7) return 20;
    else if(((detid>>14)&0xF)==8) return 21;
    else if(((detid>>14)&0xF)==9) return 22;
  }*/
  else if(subdetid==6) // TEC by ring 
  {
    if(((detid>>5)&0x7)==1) return 14;
    else if(((detid>>5)&0x7)==2) return 15;
    else if(((detid>>5)&0x7)==3) return 16;
    else if(((detid>>5)&0x7)==4) return 17;
    else if(((detid>>5)&0x7)==5) return 18;
    else if(((detid>>5)&0x7)==6) return 19;
    else if(((detid>>5)&0x7)==7) return 20;
    else return -1;
  }
  return -1;
}

int getIndex(vector<int> v, int K) 
{ 
  auto it = find(v.begin(), v.end(), K); 
  int index = it - v.begin(); 

  return index; 
} 

int Correction_FL_FR_xtalk(const std::vector <int>&  Q, int layer, std::string TemplateFile)
{
  int MaxCorr = 0;
  int cpt_sat = 0;
  int ThresholdSat = -1;

  if ((layer>=5 && layer<=10) || layer>=18) ThresholdSat = 25;
  else ThresholdSat = 16;
  
  // CORRECTION TEMPLATES
  std::string line;
  std::vector<double> template_FL;
  std::vector<double> template_FR;

  std::ifstream Template(TemplateFile);
  while (std::getline(Template, line))
  {
    std::istringstream iss(line);
    int i_layer=-1;
    double FL_a=-1, FR_a=-1;
    if (!(iss >> i_layer >> FL_a >> FR_a)) break; // error
      
    template_FL.push_back(FL_a);
    template_FR.push_back(FR_a);
  }
  Template.close();

  // NUMBER OF MAX
  for (unsigned int i=0;i<Q.size();i++)
  {
    if (Q[i]>253) cpt_sat++;
  }
  cpt_sat=1;      //A ENLEVER SI ON VEUT REACTIVER LA SATURATION

  // NO CORRECTION IF TOO SMALL OR LARGE
  if (Q.size()<2 || Q.size()>8) return -1;

  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
  
  if (cpt_sat == 1)
  {
    // FULL LEFT
    if (*mQ==Q[0] && *(mQ+1)>=ThresholdSat && *(mQ+1)<254) return template_FL[layer-1] * (*(mQ+1));

    // FULL RIGHT
    if (*mQ==Q[Q.size()-1] && *(mQ-1)>=ThresholdSat && *(mQ-1)<254) return template_FR[layer-1] * (*(mQ-1)); 
    
    return -1;
  }

  // NO SATURATION --> no x-talk inversion
  else return -1;
}

int Correction_wNeighbours(const std::vector <int>&  Q, int layer, std::string TemplateFile, bool CENTER)
{
	// SETUP 
	int cpt_sat = 0;
	int Vmin = 0, Vmax = 0;
  for (unsigned int i=0;i<Q.size();i++)
  {
    if (Q[i]>=254) cpt_sat++;
  }
  
	vector<int>::const_iterator max_Q = max_element(Q.begin(), Q.end());
	int i_max = getIndex(Q,*(max_Q));

  //  A ENLEVER SI ON VEUT REACTIVER LA SATURATION

	/*if (Q.size()<=2 || cpt_sat!=1 || *(max_Q)<254 || i_max==0 || i_max==Q.size()-1)
  {
    cout << "Problem correction on LEFT, RIGHT or CENTER clusters" << endl;
    return -1;
  }*/
  if (false) return -1;
	else
	{
		if (*(max_Q-1) > *(max_Q+1))
		{
			Vmax = *(max_Q-1);
			Vmin = *(max_Q+1);
		}
		else
		{
			Vmax = *(max_Q+1);
			Vmin = *(max_Q-1);
		}
	}

  // CORRECTION TEMPLATES
  std::ifstream Template(TemplateFile);
  std::string line;
  std::vector<double> template_a1;
  std::vector<double> template_a2;
  
  if (CENTER)
  {
    while (std::getline(Template, line))
    {
      std::istringstream iss(line);
      int i_layer=-1;
      double a1=-1;
      if (!(iss >> i_layer >> a1)) break; // error
        
      template_a1.push_back(a1);
    }
    Template.close();

    return 1.*(Vmax+Vmin)/2 * template_a1[layer-1];
  }

  while (std::getline(Template, line))
  {
    std::istringstream iss(line);
    int i_layer=-1;
    double a1=-1, a2=-1;
    if (!(iss >> i_layer >> a1 >> a2 )) break; // error
      
    template_a1.push_back(a1);
    template_a2.push_back(a2);
  }
  Template.close();

  return Vmax * template_a1[layer-1] + Vmin * template_a2[layer-1];
}

int Correction_OLD(const std::vector <int>&  Q, std::vector <int>& cpt_cout, bool ok)
{
  int QII=-1;
  float thresholdSat=25;

  cpt_cout[8]++;

  if(Q.size()<2 || Q.size()>8)
  {
    cpt_cout[0]++;
    return QII;
  }

  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
  
  if (*mQ>253)
  {
    if(*mQ==255 && *(mQ-1)>253)
    {
      cpt_cout[1]++;
      return QII;
    }
    if(*mQ==255 && *(mQ+1)>253)
    {
      cpt_cout[2]++;
      return QII;
    }

    if (*(mQ-1)>thresholdSat)
    {
      if (*(mQ+1)>thresholdSat)
      { 
        if (*(mQ-1)<254)
        {
          if (*(mQ+1)<254)
          {      
            if (abs(*(mQ-1) - *(mQ+1))<40)
            {
              QII = (10*(*(mQ-1)) + 10*(*(mQ+1))) / 2;
              return QII;
            }
            else cpt_cout[7]++;
            
          }
          else cpt_cout[6]++;
        }
        else cpt_cout[5]++;
      }
      else cpt_cout[4]++;
    }
    if (ok) for (int i=0; i<Q.size(); i++) cout << Q[i] << " ";
    if (ok) cout << endl;
    if (ok) cout << *(mQ-2) << " " << *(mQ-1) << " " << *mQ << " " << *(mQ+1) << " " << *(mQ+2) << endl;
    else cpt_cout[3]++;
  }

  return QII; // no saturation --> no x-talk inversion
}

float ClusterBarycenter(const std::vector <int>&  Q)
{
  float barycenter=0;
  float sum=0;

  for (int i=0; i<Q.size(); i++)
  {
    barycenter+=Q[i]*i;
    sum+=Q[i];
  }
  barycenter/=sum;

  return barycenter;
}

std::vector <int> Cluster_Corrected(const std::vector <int>& Qold, int Maxcorr)
{
  std::vector <int> Qnew;
  int Max_old = *max_element(Qold.begin(), Qold.end());
  
  for (int i=0; i<Qold.size(); i++)
  {
    if (Qold[i] != Max_old) Qnew.push_back(Qold[i]);
    else Qnew.push_back(Maxcorr);
  }

  return Qnew;
}

void SaturationCorrection::Loop()
{
  gROOT->SetBatch(kTRUE);

  TFile* SaveData = new TFile("ROOT_SVG/Check_Correction.root", "RECREATE");
  TFile* SaveDisplay = new TFile("ROOT_SVG/Check_Display.root", "RECREATE");
  TFile* SaveDataRaph = new TFile("ROOT_SVG/PourRaph.root", "RECREATE");
  ofstream Detid_overcorrected("ROOT_SVG/Detid_overcorrected.txt", std::ofstream::out);
  //ofstream Detid_correction("ROOT_SVG/Detid_correction.txt", std::ofstream::out);

  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  double DeltaR=0;
  double DeltaR_val;
  int index;

  int cpt_FR_bad=0, cpt_FL_bad=0, cpt_FR_bad_id=0, cpt_FL_bad_id=0, cpt_LR_bad=0, cpt_center_bad=0;
  int cpt_FR=0, cpt_FL=0;
  int cpt_all=0, cpt_corrected=0;
  int cpt_size1=0;
  int cpt_clusters=0, cpt_clusters_1strip=0, cpt_clusters_2strip=0, cpt_clusters_3strip=0, cpt_clusters_4strip=0, cpt_clusters_MoreThan4strip=0;
  bool bool_LR_bad=false;
  double local_theta_LR_bad = -1;

  int loctheta_00_05=0, loctheta_05_10=0, loctheta_10_15=0, loctheta_15_20=0, loctheta_20_25=0, loctheta_25_30=0;

  int cpt_120=0, cpt_110=0, cpt_100=0, cpt_90=0, cpt_80=0, cpt_60=0, cpt_none=0;
  int cpt_120_bad=0, cpt_110_bad=0, cpt_100_bad=0, cpt_90_bad=0, cpt_80_bad=0, cpt_60_bad=0, cpt_none_bad=0;

  std::vector <int> cpt_cout;
  std::vector <int> cpt_cout_FL;
  std::vector <int> cpt_cout_FR;
  std::vector <int> cpt_cout_L;
  std::vector <int> cpt_cout_R;
  std::vector <int> cpt_cout_C;
  
  for (int i=0; i<=8; i++)
  {
    cpt_cout.push_back(0);
    cpt_cout_FL.push_back(0);
    cpt_cout_FR.push_back(0);
    cpt_cout_L.push_back(0);
    cpt_cout_R.push_back(0);
    cpt_cout_C.push_back(0);
  }

  TH1F *QsigWsat_OVER_Qsig=new TH1F("QsigWsat_OVER_Qsig","QsigWsat_OVER_Qsig",200,0,3);//0.2,1.2);
  TH1F *Qcorr_OVER_Q_OLD=new TH1F("Qcorr_OVER_Q_OLD","Qcorr_OVER_Q_OLD",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_all=new TH1F("Qcorr_OVER_Q_all","Qcorr_OVER_Q_all",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_recoIn220_240=new TH1F("Qcorr_OVER_Q_recoIn220_240","Qcorr_OVER_Q_recoIn220_240",200,0.5,1.5);

  TH1F *Qcorr_OVER_Q_FLFR=new TH1F("Qcorr_OVER_Q_FLFR","Qcorr_OVER_Q_FLFR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_FL=new TH1F("Qcorr_OVER_Q_FL","Qcorr_OVER_Q_FL",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_FR=new TH1F("Qcorr_OVER_Q_FR","Qcorr_OVER_Q_FR",200,0,3);//0.5,1.5);

  TH1F *Qcorr_OVER_Q_LR=new TH1F("Qcorr_OVER_Q_LR","Qcorr_OVER_Q_LR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TIB_LR=new TH1F("Qcorr_OVER_Q_TIB_LR","Qcorr_OVER_Q_TIB_LR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TOB_LR=new TH1F("Qcorr_OVER_Q_TOB_LR","Qcorr_OVER_Q_TOB_LR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TID_LR=new TH1F("Qcorr_OVER_Q_TID_LR","Qcorr_OVER_Q_TID_LR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TEC_LR=new TH1F("Qcorr_OVER_Q_TEC_LR","Qcorr_OVER_Q_TEC_LR",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_LR_corrected=new TH1F("Qcorr_OVER_Q_LR_corrected","Qcorr_OVER_Q_LR_corrected",200,0,3);//0.5,1.5);

  TH1F *Qcorr_OVER_Q_CENTER=new TH1F("Qcorr_OVER_Q_CENTER","Qcorr_OVER_Q_CENTER",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TIB_CENTER=new TH1F("Qcorr_OVER_Q_TIB_CENTER","Qcorr_OVER_Q_TIB_CENTER",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TOB_CENTER=new TH1F("Qcorr_OVER_Q_TOB_CENTER","Qcorr_OVER_Q_TOB_CENTER",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TID_CENTER=new TH1F("Qcorr_OVER_Q_TID_CENTER","Qcorr_OVER_Q_TID_CENTER",200,0,3);//0.5,1.5);
  TH1F *Qcorr_OVER_Q_TEC_CENTER=new TH1F("Qcorr_OVER_Q_TEC_CENTER","Qcorr_OVER_Q_TEC_CENTER",200,0,3);//0.5,1.5);


  TH2F *eta_VS_size_LR_good = new TH2F("eta_VS_size_LR_good","eta_VS_size_LR_good",100,-3,3, 20,0,20);
  TH2F *eta_VS_size_LR_bad = new TH2F("eta_VS_size_LR_bad","eta_VS_size_LR_bad",100,-3,3, 20,0,20);
  TH2F *MaxcorrOverMax_VS_size_LR = new TH2F("MaxcorrOverMax_VS_size_LR","MaxcorrOverMax_VS_size_LR",30,0,30, 50,0.5,1.5);
  TH2F *MaxcorrOverMax_VS_eta_LR = new TH2F("MaxcorrOverMax_VS_eta_LR","MaxcorrOverMax_VS_eta_LR",30,-3,3, 50,0.5,1.5);

  TH2F *QcorrOQsum_VS_V2OV1 = new TH2F("QcorrOQsum_VS_V2OV1","QcorrOQsum_VS_V2OV1",50,0,2, 50,0,5);
  TH2F *QcorrOQsum_VS_V1 = new TH2F("QcorrOQsum_VS_V1","QcorrOQsum_VS_V1",26,0,255, 50,0,5);
  TH2F *V2OV1_VS_V1 = new TH2F("V2OV1_VS_V1","V2OV1_VS_V1",26,0,255, 50,0,2);
  TH3F *QcorrOQsum_VS_V2OV1_VS_V1 = new TH3F("QcorrOQsum_VS_V2OV1_VS_V1","QcorrOQsum_VS_V2OV1_VS_V1",26,0,255, 80,0,0.8, 50,0,5);
  QcorrOQsum_VS_V2OV1_VS_V1->GetXaxis()->SetTitle("1^{er} voisin");
  QcorrOQsum_VS_V2OV1_VS_V1->GetYaxis()->SetTitle("2^{nd} voisin / 1^{er} voisin");
  QcorrOQsum_VS_V2OV1_VS_V1->GetZaxis()->SetTitle("Qcorr/Qsum");
  TH3F *V2OV1_VS_V1h_V1l_QcOQmore3 = new TH3F("V2OV1_VS_V1h_V1l_QcOQmore3","V2OV1_VS_V1h_V1l_QcOQmore3",26,0,255, 26,0,255, 80,0,0.8);
  V2OV1_VS_V1h_V1l_QcOQmore3->GetXaxis()->SetTitle("1^{er} voisin haut");
  V2OV1_VS_V1h_V1l_QcOQmore3->GetYaxis()->SetTitle("1^{er} voisin bas");
  V2OV1_VS_V1h_V1l_QcOQmore3->GetZaxis()->SetTitle("2^{nd} voisin / 1^{er} voisin");
  TH3F *V2OV1_VS_V1h_V1l_QcOQnear1 = new TH3F("V2OV1_VS_V1h_V1l_QcOQnear1","V2OV1_VS_V1h_V1l_QcOQnear1",26,0,255, 26,0,255, 80,0,0.8);
  V2OV1_VS_V1h_V1l_QcOQnear1->GetXaxis()->SetTitle("1^{er} voisin haut");
  V2OV1_VS_V1h_V1l_QcOQnear1->GetYaxis()->SetTitle("1^{er} voisin bas");
  V2OV1_VS_V1h_V1l_QcOQnear1->GetZaxis()->SetTitle("2^{nd} voisin / 1^{er} voisin");
  TH3F *QcorrOQsum_VS_V1p_VS_V1m = new TH3F("QcorrOQsum_VS_V1p_VS_V1m","QcorrOQsum_VS_V1p_VS_V1m",26,0,255, 26,0,255, 50,0,5);
  QcorrOQsum_VS_V1p_VS_V1m->GetXaxis()->SetTitle("1^{er} voisin haut");
  QcorrOQsum_VS_V1p_VS_V1m->GetYaxis()->SetTitle("1^{er} voisin bas");
  QcorrOQsum_VS_V1p_VS_V1m->GetZaxis()->SetTitle("Qcorr/Qsum");

  TH1F *localtheta = new TH1F("localtheta","localtheta",100,-Pi()/2-0.05,Pi()/2+0.05);
  TH2F *BadCorrection_layer = new TH2F("BadCorrection_layer","BadCorrection_layer",21,0,21, 10,0,10);

  TH1F *ClusterSat_layer = new TH1F("ClusterSat_layer","ClusterSat_layer",21,0,21);
  TH1F *Cluster_layer = new TH1F("Cluster_layer","Cluster_layer",21,0,21);
  TH1F *SatValueADC = new TH1F("SatValueADC","SatValueADC",250,10,1200);

  TH2F *barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue = new TH2F("barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue","barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue",100,0.7,1.3,100,0.7,1.3);
  barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->GetXaxis()->SetTitle("bar(Qtrue with sat) / bar(Qtrue)");
  barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->GetYaxis()->SetTitle("bar (Qcorr) / bar (Qtrue)");


  //nentries = 10000;
	for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

        // MUON MATCHING
    int gen_index=0;
    for (int igen=0; igen<ngenpart; igen++)
    {
      if (gen_pdg[igen]==13 && gen_status[igen]==1) gen_index=igen;
    }

    DeltaR_val=100;
    index=-1;
    for (int itr=0; itr<ntracks; itr++)
    {
      DeltaR=Sqrt(Power(modulo2Pi(track_phi[itr]-gen_phi[gen_index]),2)+Power(track_eta[itr]-gen_eta[gen_index],2));

      if(DeltaR<DeltaR_val)
      {
        DeltaR_val=DeltaR;
        index=itr;
      }
    }


        // FULL CLUSTERS ANALYSIS
    if (track_nvalidhits[index]>8 && track_qual[index])
    {
      for (int idedx=track_index_hit[index]; idedx<track_index_hit[index]+track_nhits[index]; idedx++)
      {
        if (dedx_isstrip[idedx])
        { 
              // LAYER
          int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

              // SETUP
          vector<int> index_sig;
          vector<int> value_sig;
          vector<int> value_sig_nosat;
          double value_max=-1;
          int i_max=-1, pos_max=-1;
          int sum_Qtrue = 0, sum_QclusterSAT = 0;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
          bool_LR_bad = false;
          int Nleft = 0, Nright = 0;

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              index_sig.push_back(isigdigi);
              value_sig.push_back(sigdigi_strip_adc[isigdigi]);
            
              if (sigdigi_strip_adc[isigdigi]<=253) value_sig_nosat.push_back(sigdigi_strip_adc[isigdigi]);
              else if (sigdigi_strip_adc[isigdigi]>253 && sigdigi_strip_adc[isigdigi]<=1023) value_sig_nosat.push_back(254);
              else if (sigdigi_strip_adc[isigdigi]>1023) value_sig_nosat.push_back(255);
            }
          }
          /*for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            value_sig_nosat.push_back(strip_ampl[istrip]);
            index_sig.push_back(strip_channel[istrip]);
          }*/
        

          if (!index_sig.empty() && !value_sig.empty() && !value_sig_nosat.empty())
          {
            value_max = *max_element(value_sig_nosat.begin(), value_sig_nosat.end());
            i_max = index_sig.front() + getIndex(value_sig_nosat,value_max);  // max index : sigdigi_strip_adc[i_max]
            pos_max = getIndex(value_sig_nosat,value_max);  // max index relative to the cluster : value_sig[pos_max], value_sig_nosat[pos_max]
            sum_Qtrue = accumulate(value_sig.begin(), value_sig.end(), 0);
            sum_QclusterSAT = accumulate(value_sig_nosat.begin(), value_sig_nosat.end(), 0);
            if (value_sig_nosat.size()>=3) {Nleft = value_sig_nosat[pos_max-1]; Nright = value_sig_nosat[pos_max+1];}

            if (sum_Qtrue==0) sum_Qtrue=-1;

            int check_OneMaxSat=0;
            for (int i=0; i<value_sig_nosat.size(); i++)
            {
              if (value_sig_nosat[i]>253) check_OneMaxSat++;
            }
            cpt_clusters++;
            if (value_sig_nosat.size()==1 && check_OneMaxSat==1) cpt_size1++;
            if (check_OneMaxSat==1) cpt_clusters_1strip++;
            if (check_OneMaxSat==2) cpt_clusters_2strip++;
            if (check_OneMaxSat==3) cpt_clusters_3strip++;
            if (check_OneMaxSat==4) cpt_clusters_4strip++;
            if (check_OneMaxSat>4) cpt_clusters_MoreThan4strip++;
            

            if (index_sig.size()>=3 && check_OneMaxSat==1 && pos_max>0 && pos_max<index_sig.size()-1 && Nleft>1.1*Nright) left=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && pos_max>0 && pos_max<index_sig.size()-1 && Nleft<0.9*Nright) right=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && i_max>0 && pos_max<index_sig.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && pos_max==0) FullLeft=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && pos_max==index_sig.size()-1) FullRight=true;
            

                // DISPLAY RAPH
            if (value_sig.size()>=3 && value_sig[i_max-index_sig.front()]<=250)
            {
              if (sclus_loctheta[idedx]>=0 && sclus_loctheta[idedx]<0.5 && loctheta_00_05<=2)
              { loctheta_00_05++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx], 0, false); }
              if (sclus_loctheta[idedx]>=0.5 && sclus_loctheta[idedx]<1 && loctheta_05_10<=2)
              { loctheta_05_10++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx], 0, false); }
              if (sclus_loctheta[idedx]>=1 && sclus_loctheta[idedx]<Pi()/2 && loctheta_10_15<=2)
              { loctheta_10_15++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx], 0, false); }
              if (sclus_loctheta[idedx]>=Pi()/2 && sclus_loctheta[idedx]<2 && loctheta_15_20<=2)
              { loctheta_15_20++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx]-Pi(), 0, false); }
              if (sclus_loctheta[idedx]>=2 && sclus_loctheta[idedx]<2.5 && loctheta_20_25<=2)
              { loctheta_20_25++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx]-Pi(), 0, false); }
              if (sclus_loctheta[idedx]>=2.5 && sclus_loctheta[idedx]<-Pi() && loctheta_25_30<=2)
              { loctheta_25_30++; DisplayCluster(value_sig, SaveDataRaph, jentry, index, layer, sclus_loctheta[idedx]-Pi(), 0, false); }
            }
            if (sclus_loctheta[idedx]<Pi()/2) localtheta->Fill(sclus_loctheta[idedx]);
            else localtheta->Fill(sclus_loctheta[idedx]-Pi());
           

                // Old correction
            if (check_OneMaxSat>0) ClusterSat_layer->Fill(layer);
            Cluster_layer->Fill(layer);
            if (check_OneMaxSat>0) SatValueADC->Fill(value_max);

            if (FullLeft || FullRight || left || right || center)
            {
              int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout, false);

              cpt_all++;
              if (Qcorr!=-1)
              {
                Qcorr_OVER_Q_OLD->Fill((double) Qcorr/sum_Qtrue);
                cpt_corrected++;
              }

              QsigWsat_OVER_Qsig->Fill((double) sum_QclusterSAT / sum_Qtrue);
            }
            if (FullLeft) int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout_FL, false);
            if (FullRight) int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout_FR, false);
            if (left) int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout_L, false);
            if (right) int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout_R, false);
            if (center) int Qcorr = Correction_OLD(value_sig_nosat, cpt_cout_C, false);


                // dE/dx in [220 ; 240]
            bool recoIn220_240 = false;
            if (value_sig_nosat[pos_max] >= 220 && value_sig_nosat[pos_max] <= 240) recoIn220_240 = true;

            if (index_sig.size()>=3 && recoIn220_240 && pos_max>0 && pos_max<index_sig.size()-1 && (Nleft>1.1*Nright || Nleft<0.9*Nright))
            {
              int Maxcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
              Qcorr_OVER_Q_recoIn220_240->Fill(1.*Maxcorr/value_sig_nosat[pos_max]);
            }
            if (index_sig.size()>=3 && recoIn220_240 && i_max>0 && pos_max<index_sig.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright)
            {
              int Maxcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
              Qcorr_OVER_Q_recoIn220_240->Fill(1.*Maxcorr/value_sig_nosat[pos_max]);
            }
            if (index_sig.size()>=2 && recoIn220_240 && (pos_max==0 || pos_max==index_sig.size()-1))
            {
              int Maxcorr = Correction_FL_FR_xtalk(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");
              Qcorr_OVER_Q_recoIn220_240->Fill(1.*Maxcorr/value_sig_nosat[pos_max]);
            }


                // FL and FR correction
            if (FullLeft)
            {
              cpt_FL++;
							int Qcorr = Correction_FL_FR_xtalk(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");

              if (Qcorr == -1) cpt_FL_bad_id++;
              if (Qcorr < 254) cpt_FL_bad++;
              else
              {
                barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->Fill(ClusterBarycenter(value_sig_nosat)/ClusterBarycenter(value_sig), ClusterBarycenter(Cluster_Corrected(value_sig_nosat, Qcorr))/ClusterBarycenter(value_sig));

                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                }
                Qcorr_OVER_Q_FL->Fill((double) Qcorr / sum_Qtrue);
                
                //Detid_correction << dedx_detid[idedx] << " " << Qcorr*1./sum_Qtrue << "\n";
              }
						}
						if (FullRight)
						{
              cpt_FR++;
							int Qcorr = Correction_FL_FR_xtalk(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");

              if (Qcorr == -1) cpt_FR_bad_id++;
							if (Qcorr < 254) cpt_FR_bad++;
							else
							{
                barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->Fill(ClusterBarycenter(value_sig_nosat)/ClusterBarycenter(value_sig), ClusterBarycenter(Cluster_Corrected(value_sig_nosat, Qcorr))/ClusterBarycenter(value_sig));

                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                }
                Qcorr_OVER_Q_FR->Fill((double) Qcorr / sum_Qtrue);
                
                //Detid_correction << dedx_detid[idedx] << " " << Qcorr*1./sum_Qtrue << "\n";
              }
            }
             
                // LEFT and RIGHT correction
            if (left || right)
            {
              int Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);

              if (Qcorr < 254) cpt_LR_bad++;
              else
              {
                barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->Fill(ClusterBarycenter(value_sig_nosat)/ClusterBarycenter(value_sig), ClusterBarycenter(Cluster_Corrected(value_sig_nosat, Qcorr))/ClusterBarycenter(value_sig));

                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                }
                Qcorr_OVER_Q_LR->Fill((double) Qcorr / sum_Qtrue);
                
                MaxcorrOverMax_VS_size_LR->Fill(value_sig_nosat.size(), (double) Qcorr / sum_Qtrue);
                MaxcorrOverMax_VS_eta_LR->Fill(track_eta[index], (double) Qcorr / sum_Qtrue);

                //Detid_correction << dedx_detid[idedx] << " " << Qcorr*1./sum_Qtrue << "\n";

                
                if (value_sig_nosat[pos_max-1] < 120 && value_sig_nosat[pos_max+1] < 120){ cpt_120++; if (Qcorr >= 1.1*sum_Qtrue) cpt_120_bad++; }
                if (value_sig_nosat[pos_max-1] < 110 && value_sig_nosat[pos_max+1] < 110){ cpt_110++; if (Qcorr >= 1.1*sum_Qtrue) cpt_110_bad++; }
                if (value_sig_nosat[pos_max-1] < 100 && value_sig_nosat[pos_max+1] < 100){ cpt_100++; if (Qcorr >= 1.1*sum_Qtrue) cpt_100_bad++; }
                if (value_sig_nosat[pos_max-1] < 90 && value_sig_nosat[pos_max+1] < 90){ cpt_90++; if (Qcorr >= 1.1*sum_Qtrue) cpt_90_bad++; }
                if (value_sig_nosat[pos_max-1] < 80 && value_sig_nosat[pos_max+1] < 80){ cpt_80++; if (Qcorr >= 1.1*sum_Qtrue) cpt_80_bad++; }
                if (value_sig_nosat[pos_max-1] < 60 && value_sig_nosat[pos_max+1] < 60){ cpt_60++; if (Qcorr >= 1.1*sum_Qtrue) cpt_60_bad++; }
                cpt_none++; if(Qcorr >= 1.1*sum_Qtrue) cpt_none_bad++;


                if (left && pos_max-2>=0)
                {
                  if (value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1]>0.1 && value_sig_nosat[pos_max+1]<100) Qcorr_OVER_Q_LR_corrected->Fill((double) Qcorr / sum_Qtrue);

                  QcorrOQsum_VS_V2OV1_VS_V1->Fill(value_sig_nosat[pos_max-1], value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1], Qcorr*1./sum_Qtrue);
                  QcorrOQsum_VS_V2OV1->Fill(value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1], Qcorr*1./sum_Qtrue);
                  QcorrOQsum_VS_V1->Fill(value_sig_nosat[pos_max-1], Qcorr*1./sum_Qtrue);
                  V2OV1_VS_V1->Fill(value_sig_nosat[pos_max-1], value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1]);
                  QcorrOQsum_VS_V1p_VS_V1m->Fill(value_sig_nosat[pos_max-1], value_sig_nosat[pos_max+1], Qcorr*1./sum_Qtrue);

                  if (Qcorr*1./sum_Qtrue > 3) V2OV1_VS_V1h_V1l_QcOQmore3->Fill(value_sig_nosat[pos_max-1], value_sig_nosat[pos_max+1], value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1]);
                  if (Qcorr*1./sum_Qtrue < 1.08 && Qcorr*1./sum_Qtrue > 0.92) V2OV1_VS_V1h_V1l_QcOQnear1->Fill(value_sig_nosat[pos_max-1], value_sig_nosat[pos_max+1], value_sig_nosat[pos_max-2]*1./value_sig_nosat[pos_max-1]);
                }
                if (right && pos_max+2<=value_sig_nosat.size()-1)
                { 
                  if (value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1]>0.1 && value_sig_nosat[pos_max-1]<100) Qcorr_OVER_Q_LR_corrected->Fill((double) Qcorr / sum_Qtrue);


                  QcorrOQsum_VS_V2OV1_VS_V1->Fill(value_sig_nosat[pos_max+1], value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1], Qcorr*1./sum_Qtrue);
                  QcorrOQsum_VS_V2OV1->Fill(value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1], Qcorr*1./sum_Qtrue);
                  QcorrOQsum_VS_V1->Fill(value_sig_nosat[pos_max+1], Qcorr*1./sum_Qtrue);
                  V2OV1_VS_V1->Fill(value_sig_nosat[pos_max+1], value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1]);
                  QcorrOQsum_VS_V1p_VS_V1m->Fill(value_sig_nosat[pos_max+1], value_sig_nosat[pos_max-1], Qcorr*1./sum_Qtrue);

                  if (Qcorr*1./sum_Qtrue > 3) V2OV1_VS_V1h_V1l_QcOQmore3->Fill(value_sig_nosat[pos_max+1], value_sig_nosat[pos_max-1], value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1]);
                  if (Qcorr*1./sum_Qtrue < 1.08 && Qcorr*1./sum_Qtrue > 0.92) V2OV1_VS_V1h_V1l_QcOQnear1->Fill(value_sig_nosat[pos_max+1], value_sig_nosat[pos_max-1], value_sig_nosat[pos_max+2]*1./value_sig_nosat[pos_max+1]);
                }
              }

              if (Qcorr >= 1.1*sum_Qtrue) 
              {
                eta_VS_size_LR_bad->Fill(track_eta[index],value_sig.size());
                BadCorrection_layer->Fill(layer, Qcorr/sum_Qtrue);
              }
              else eta_VS_size_LR_good->Fill(track_eta[index],value_sig.size());
              if (Qcorr >= 4*sum_Qtrue)
              { 
                //bool_LR_bad = true;
                local_theta_LR_bad = sclus_loctheta[idedx];

                Detid_overcorrected << dedx_detid[idedx] << " " << Qcorr*1./sum_Qtrue << "\n";
              }

                  // layer by layer
              if (dedx_subdetid[idedx]==3)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TIB_LR->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==5)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TOB_LR->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==4)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TID_LR->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==6)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
                  
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TEC_LR->Fill((double) Qcorr / sum_Qtrue);
                }
              }
            }
            
                // CENTER correction              
            if (center)
            {
              int Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);

              if (Qcorr < 254) cpt_center_bad++;
              else
              {
                barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->Fill(ClusterBarycenter(value_sig_nosat)/ClusterBarycenter(value_sig), ClusterBarycenter(Cluster_Corrected(value_sig_nosat, Qcorr))/ClusterBarycenter(value_sig));

                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                }
                Qcorr_OVER_Q_CENTER->Fill((double) Qcorr / sum_Qtrue);
                
                //Detid_correction << dedx_detid[idedx] << " " << Qcorr*1./sum_Qtrue << "\n";
              }

              if (false) //Qcorr >= 4*sum_Qtrue && jentry <= 10000)
              { 
                bool_LR_bad = true;

                std::vector <int> test_digi;
                cout << "layer: " << layer << endl;
                cout << "RECO: ";
                for (int i=0; i<value_sig_nosat.size(); i++) cout << value_sig_nosat[i] << " ";
                cout << endl;
                cout << "DIGI: ";
                for (int i=0; i<value_sig.size(); i++) cout << value_sig[i] << " ";
                cout << endl;
                cout << "Correction RECO: " << Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true) << endl;
                cout << "Correction DIGI: " << Correction_wNeighbours(value_sig, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true) << endl;
                cout << endl;
                cout << "----------" << endl;
              }
              
                  // layer by layer
              if (dedx_subdetid[idedx]==3)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TIB_CENTER->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==5)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TOB_CENTER->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==4)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TID_CENTER->Fill((double) Qcorr / sum_Qtrue);
                }
              }
              else if (dedx_subdetid[idedx]==6)
              {
                Qcorr = Correction_wNeighbours(value_sig_nosat, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
                
                if (Qcorr >= 254)
                {
                  for (int i=0; i<value_sig_nosat.size(); i++)
                  {
                    if (i != i_max-index_sig.front()) Qcorr += value_sig_nosat[i];
                  }
                  Qcorr_OVER_Q_TEC_CENTER->Fill((double) Qcorr / sum_Qtrue);
                }
              }
            }
 

					} //end if index_sig

          std::vector<int> sigreco;
          std::vector<int> sigreco_channel;
          std::vector<int> sigdigi;
          std::vector<int> sigdigi_channel;
          std::vector<int> sumdigi;
          std::vector<int> sumdigi_channel;

          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            sigreco.push_back(strip_ampl[istrip]);
            sigreco_channel.push_back(strip_channel[istrip]);
          }

          for (int isigdigi=0; isigdigi<ndigi_sig_strip; isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]]-3 && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1]+3)
            {
              sigdigi.push_back(sigdigi_strip_adc[isigdigi]);
              sigdigi_channel.push_back(sigdigi_strip_channel[isigdigi]);
            }
          }

          for (int isumdigi=0; isumdigi<ndigi_sum_strip; isumdigi++)
          {
            if (sumdigi_strip_id[isumdigi]==dedx_detid[idedx] && sumdigi_strip_channel[isumdigi]>= strip_channel[sclus_index_strip[idedx]]-3 && sumdigi_strip_channel[isumdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1]+3)
            {
              sumdigi.push_back(sumdigi_strip_adc[isumdigi]);
              sumdigi_channel.push_back(sumdigi_strip_channel[isumdigi]);
            }
          }

          
          if (!sigdigi.empty() && !sigreco.empty() && !sumdigi.empty())
          {
            bool channel_miss = false; 
            for (int i=0; i<sigdigi_channel.size()-1; i++)
            {
              if (sigdigi_channel[i]!=sigdigi_channel[i+1]-1) channel_miss = true;
            }
            if (channel_miss == false)
            {
              int shift = sigreco_channel.front() - sigdigi_channel.front();
              if (jentry <= 10000 && bool_LR_bad) DisplayCluster(sigreco, sigdigi, SaveDisplay, shift, layer, local_theta_LR_bad, true);
              //DisplayCluster(sigreco, sigdigi, sumdigi, SaveDisplay, shift, jentry, index, layer, local_theta_LR_bad, true);
            }
          }



        } //end if strip
      } //end hit loop
    } //end track loop
	  
		if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  } // end entry loop

  Qcorr_OVER_Q_all->Add(Qcorr_OVER_Q_LR);
  Qcorr_OVER_Q_all->Add(Qcorr_OVER_Q_CENTER);
  Qcorr_OVER_Q_all->Add(Qcorr_OVER_Q_FL);
  Qcorr_OVER_Q_all->Add(Qcorr_OVER_Q_FR);

  Qcorr_OVER_Q_FLFR->Add(Qcorr_OVER_Q_FL);
  Qcorr_OVER_Q_FLFR->Add(Qcorr_OVER_Q_FR);



  SaveData->cd();

  QsigWsat_OVER_Qsig->Write();
  Qcorr_OVER_Q_OLD->Write();
  Qcorr_OVER_Q_recoIn220_240->Write();

  Qcorr_OVER_Q_all->Write();
  Qcorr_OVER_Q_LR->Write();
  Qcorr_OVER_Q_LR_corrected->Write();
  Qcorr_OVER_Q_CENTER->Write();
  Qcorr_OVER_Q_FLFR->Write();

  Qcorr_OVER_Q_TIB_LR->Write();
  Qcorr_OVER_Q_TOB_LR->Write();
  Qcorr_OVER_Q_TID_LR->Write();
  Qcorr_OVER_Q_TEC_LR->Write();
  Qcorr_OVER_Q_TIB_CENTER->Write();
  Qcorr_OVER_Q_TOB_CENTER->Write();
  Qcorr_OVER_Q_TID_CENTER->Write();
  Qcorr_OVER_Q_TEC_CENTER->Write();
  Qcorr_OVER_Q_FL->Write();
  Qcorr_OVER_Q_FR->Write();

  eta_VS_size_LR_good->Write();
  eta_VS_size_LR_bad->Write();
  MaxcorrOverMax_VS_size_LR->Write();
  MaxcorrOverMax_VS_eta_LR->Write();

  QcorrOQsum_VS_V2OV1->Write();
  QcorrOQsum_VS_V1->Write();
  V2OV1_VS_V1->Write();
  QcorrOQsum_VS_V2OV1_VS_V1->Write();

  localtheta->Write();


  V2OV1_VS_V1h_V1l_QcOQmore3->Write();
  V2OV1_VS_V1h_V1l_QcOQnear1->Write();

  BadCorrection_layer->Write();
  QcorrOQsum_VS_V1p_VS_V1m->Write();

  ClusterSat_layer->Sumw2();
  Cluster_layer->Sumw2();
  ClusterSat_layer->Divide(Cluster_layer);
  ClusterSat_layer->Write();

  SatValueADC->Write();

  barycenter_QtrueWsatOVERQtrue_VS_QcorrOVERQtrue->Write();

  Detid_overcorrected.close();
  SaveData->Close();

  cout << endl;
  cout << "-----------------------------------" << endl;
  cout << "      CLUSTERS CHARACTERISTCS      " << endl;
  cout << "-----------------------------------" << endl;
  cout << "Number of clusters : " << cpt_clusters << endl;
  cout << "With size=1 : " << cpt_size1 << endl;
  cout << "With 1 saturated strip : " << cpt_clusters_1strip << endl;
  cout << "With 2 saturated strips : " << cpt_clusters_2strip << endl;
  cout << "With 3 saturated strips : " << cpt_clusters_3strip << endl;
  cout << "With 4 saturated strips : " << cpt_clusters_4strip << endl;
  cout << "With more than 4 saturated strips : " << cpt_clusters_MoreThan4strip << endl;
  cout << endl;

  cout << "-----------------------------------" << endl;
  cout << "            OLD CORRECTION         " << endl;
  cout << "-----------------------------------" << endl;
  cout << "Clusters corrected : " << cpt_corrected << " / " << cpt_all << " (" << 100*(double) cpt_corrected/cpt_all << " %)" << endl;
  cout << "The way it's rejected: " << endl;
  cout << "all: "; for (int i=0; i<cpt_cout.size(); i++) cout << cpt_cout[i] << " ";
  cout << endl;
  cout << "FL: "; for (int i=0; i<cpt_cout_FL.size(); i++) cout << cpt_cout_FL[i] << " ";
  cout << endl;
  cout << "FR: "; for (int i=0; i<cpt_cout_FR.size(); i++) cout << cpt_cout_FR[i] << " ";
  cout << endl;
  cout << "L: "; for (int i=0; i<cpt_cout_L.size(); i++) cout << cpt_cout_L[i] << " ";
  cout << endl;
  cout << "R: "; for (int i=0; i<cpt_cout_R.size(); i++) cout << cpt_cout_R[i] << " ";
  cout << endl;
  cout << "C: "; for (int i=0; i<cpt_cout_C.size(); i++) cout << cpt_cout_C[i] << " ";
  cout << endl;
  cout << endl;
  
  cout << "-----------------------------------" << endl;
  cout << "            NEW CORRECTION         " << endl;
  cout << "-----------------------------------" << endl;
  cout << "Clusters not corrected: " << endl;
  cout << "CENTER: " << cpt_center_bad << " (Correction<254) / " << Qcorr_OVER_Q_CENTER->GetEntries() << " (" << 100*(double) cpt_center_bad/Qcorr_OVER_Q_CENTER->GetEntries() << " %)" << endl;
  cout << "LEFT and RIGHT: " << cpt_LR_bad << " (Correction<254) / " << Qcorr_OVER_Q_LR->GetEntries() << " (" << 100*(double) cpt_LR_bad/Qcorr_OVER_Q_LR->GetEntries() << " %)" << endl;
  cout << "FL and FR: " << cpt_FL_bad_id + cpt_FR_bad_id << " (Neighbour<thld) + " << cpt_FL_bad-cpt_FL_bad_id + cpt_FR_bad-cpt_FR_bad_id << " (Correction<254) / " << Qcorr_OVER_Q_FLFR->GetEntries() << " (" << 100*(double) (cpt_FL_bad+cpt_FR_bad)/Qcorr_OVER_Q_FLFR->GetEntries() << " %)" << endl;
  cout << "Total: " << cpt_FL_bad+cpt_FR_bad+cpt_LR_bad+cpt_center_bad << " (Correction<254) / " << Qcorr_OVER_Q_FLFR->GetEntries()+Qcorr_OVER_Q_LR->GetEntries()+Qcorr_OVER_Q_CENTER->GetEntries() << " (" << 100*(double) (cpt_FL_bad+cpt_FR_bad+cpt_LR_bad+cpt_center_bad)/(Qcorr_OVER_Q_FLFR->GetEntries()+Qcorr_OVER_Q_LR->GetEntries()+Qcorr_OVER_Q_CENTER->GetEntries()) << " %)" << endl;
  cout << endl;
  cout << "Integral in [mu-sigma ; mu+sigma]: " << endl;
  cout << "CENTER: " << Qcorr_OVER_Q_CENTER->Integral(1, Qcorr_OVER_Q_CENTER->FindBin(Qcorr_OVER_Q_CENTER->GetMean() + Qcorr_OVER_Q_CENTER->GetRMS())) - Qcorr_OVER_Q_CENTER->Integral(1, Qcorr_OVER_Q_CENTER->FindBin(Qcorr_OVER_Q_CENTER->GetMean() - Qcorr_OVER_Q_CENTER->GetRMS())) << " / " << Qcorr_OVER_Q_CENTER->GetEntries() << " (" << 100.* (Qcorr_OVER_Q_CENTER->Integral(1, Qcorr_OVER_Q_CENTER->FindBin(Qcorr_OVER_Q_CENTER->GetMean() + Qcorr_OVER_Q_CENTER->GetRMS())) - Qcorr_OVER_Q_CENTER->Integral(1, Qcorr_OVER_Q_CENTER->FindBin(Qcorr_OVER_Q_CENTER->GetMean() - Qcorr_OVER_Q_CENTER->GetRMS())))/ Qcorr_OVER_Q_CENTER->GetEntries() << " %)" << endl;
  cout << "LEFT and RIGHT: " << Qcorr_OVER_Q_LR->Integral(1, Qcorr_OVER_Q_LR->FindBin(Qcorr_OVER_Q_LR->GetMean() + Qcorr_OVER_Q_LR->GetRMS())) - Qcorr_OVER_Q_LR->Integral(1, Qcorr_OVER_Q_LR->FindBin(Qcorr_OVER_Q_LR->GetMean() - Qcorr_OVER_Q_LR->GetRMS())) << " / " << Qcorr_OVER_Q_LR->GetEntries() << " (" << 100.* (Qcorr_OVER_Q_LR->Integral(1, Qcorr_OVER_Q_LR->FindBin(Qcorr_OVER_Q_LR->GetMean() + Qcorr_OVER_Q_LR->GetRMS())) - Qcorr_OVER_Q_LR->Integral(1, Qcorr_OVER_Q_LR->FindBin(Qcorr_OVER_Q_LR->GetMean() - Qcorr_OVER_Q_LR->GetRMS())))/ Qcorr_OVER_Q_LR->GetEntries() << " %)" << endl;
  cout << "FL and FR: " << Qcorr_OVER_Q_FLFR->Integral(1, Qcorr_OVER_Q_FLFR->FindBin(Qcorr_OVER_Q_FLFR->GetMean() + Qcorr_OVER_Q_FLFR->GetRMS())) - Qcorr_OVER_Q_FLFR->Integral(1, Qcorr_OVER_Q_FLFR->FindBin(Qcorr_OVER_Q_FLFR->GetMean() - Qcorr_OVER_Q_FLFR->GetRMS())) << " / " << Qcorr_OVER_Q_FLFR->GetEntries() << " (" << 100.* (Qcorr_OVER_Q_FLFR->Integral(1, Qcorr_OVER_Q_FLFR->FindBin(Qcorr_OVER_Q_FLFR->GetMean() + Qcorr_OVER_Q_FLFR->GetRMS())) - Qcorr_OVER_Q_FLFR->Integral(1, Qcorr_OVER_Q_FLFR->FindBin(Qcorr_OVER_Q_FLFR->GetMean() - Qcorr_OVER_Q_FLFR->GetRMS())))/ Qcorr_OVER_Q_FLFR->GetEntries() << " %)" << endl;
  cout << "Total: " << Qcorr_OVER_Q_all->Integral(1, Qcorr_OVER_Q_all->FindBin(Qcorr_OVER_Q_all->GetMean() + Qcorr_OVER_Q_all->GetRMS())) - Qcorr_OVER_Q_all->Integral(1, Qcorr_OVER_Q_all->FindBin(Qcorr_OVER_Q_all->GetMean() - Qcorr_OVER_Q_all->GetRMS())) << " / " << Qcorr_OVER_Q_all->GetEntries() << " (" << 100.* (Qcorr_OVER_Q_all->Integral(1, Qcorr_OVER_Q_all->FindBin(Qcorr_OVER_Q_all->GetMean() + Qcorr_OVER_Q_all->GetRMS())) - Qcorr_OVER_Q_all->Integral(1, Qcorr_OVER_Q_all->FindBin(Qcorr_OVER_Q_all->GetMean() - Qcorr_OVER_Q_all->GetRMS())))/ Qcorr_OVER_Q_all->GetEntries() << " %)" << endl;
  cout << endl;
  cout << "LEFT and RIGHT correction, fraction of clusters with Qcorr/Qsum>=1.1: " << endl;
  cout << "No V1 criteria: " << cpt_none_bad << " / " << cpt_none << " (" << 100.*cpt_none_bad/cpt_none << " %)" << endl;
  cout << "With V1<120: " << cpt_120_bad << " / " << cpt_120 << " (" << 100.*cpt_120_bad/cpt_120 << " %)" << endl;
  cout << "With V1<110: " << cpt_110_bad << " / " << cpt_110 << " (" << 100.*cpt_110_bad/cpt_110 << " %)" << endl;
  cout << "With V1<100: " << cpt_100_bad << " / " << cpt_100 << " (" << 100.*cpt_100_bad/cpt_100 << " %)" << endl;
  cout << "With V1<90: " << cpt_90_bad << " / " << cpt_90 << " (" << 100.*cpt_90_bad/cpt_90 << " %)" << endl;
  cout << "With V1<80: " << cpt_80_bad << " / " << cpt_80 << " (" << 100.*cpt_80_bad/cpt_80 << " %)" << endl;
  cout << "With V1<60: " << cpt_60_bad << " / " << cpt_60 << " (" << 100.*cpt_60_bad/cpt_60 << " %)" << endl;
  cout << endl;
  
}
