/*
    -----------------------------------------------

    Code de avec la fausse saturation ne prenant
    en compte que les clusters tq Max in [240,250]

    LES FONCTIONS SONT MODIFIEES POUR, ATTENTION

    -----------------------------------------------
*/


#define MuonPU_checkCorrection_cxx
#include "MuonPU_checkCorrection.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <stdio.h>
#include <string.h>
#include <vector>

using namespace TMath;
using namespace std;

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
    if (Q[i]>=240 && Q[i]<=250) cpt_sat++;
  }

  // NO CORRECTION IF TOO SMALL OR LARGE
  if (Q.size()<2 || Q.size()>8) return -1;

  
  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());

  // EXCLUSION CASES
  if (Q.size()<2 || Q.size()>8 || *(mQ-1) <= 15 || *(mQ+1) <= 15
      || *(mQ-1)>=240 || *(mQ+1)>=240 || abs(*(mQ-1) - *(mQ+1))>=40) return -1;

  
  if (cpt_sat == 1)
  {
    // FULL LEFT
    if (*mQ==Q[0] && *(mQ+1)>=ThresholdSat && *(mQ+1)<240) return template_FL[layer-1] * (*(mQ+1));

    // FULL RIGHT
    if (*mQ==Q[Q.size()-1] && *(mQ-1)>=ThresholdSat && *(mQ-1)<240) return template_FR[layer-1] * (*(mQ-1));
    
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
    if (Q[i]>=240 && Q[i]<=250) cpt_sat++;
  }
	vector<int>::const_iterator max_Q = max_element(Q.begin(), Q.end());
	int i_max = getIndex(Q,*(max_Q));


  // Exclusion cases
  if (Q.size()<2 || Q.size()>8 || *(max_Q-1) <= 15 || *(max_Q+1) <= 15
      || *(max_Q-1)>=240 || *(max_Q+1)>=240 || abs(*(max_Q-1) - *(max_Q+1))>=40) return -1;


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

int Correction_OLD(const std::vector <int>&  Q)
{
  int QII = -1;
  float thresholdSat=25;
  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());

  if (Q.size() < 2 || Q.size() > 8) return QII;

  if (*mQ>=240 && *mQ<=250)
  {
    if (*(mQ-1)>15 && *(mQ+1)>15 && *(mQ-1)<240 && *(mQ+1)<240 && abs(*(mQ-1) - *(mQ+1))<40)
    {
      QII = (10*(*(mQ-1)) + 10*(*(mQ+1))) / 2;
      return QII;
    }
  }

  return QII; // no saturation --> no x-talk inversion
}

int Correction_OLD_NeighboursChanged(const std::vector <int>&  Q)
{
int QII = -1;
  float thresholdSat=25;

  if (Q.size() < 2 || Q.size() > 8) return QII;

  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
  
  if (*mQ>=240 && *mQ<=250)
  {
    if (*(mQ-1)>15 && *(mQ+1)>15 && *(mQ-1)<240 && *(mQ+1)<240 && abs(*(mQ-1) - *(mQ+1))<40)
    {
      QII = (11.9*(*(mQ-1)) + 11.9*(*(mQ+1))) / 2;
      return QII;
    }
  }

  return QII; // no saturation --> no x-talk inversion
}

float ClusterBarycenter(const std::vector <int>&  Q)
{
  float barycenter=0;
  float sum=0;

  int maxQ = *max_element(Q.begin(), Q.end());
  int i_max = getIndex(Q, maxQ);

  for (int i=0; i<Q.size(); i++)
  {
    barycenter+=Q[i]*i;
    sum+=Q[i];
  }
  barycenter/=sum;

  return barycenter - i_max;
}

void CoutClusters(const std::vector <int>& Q)
{
  for (int i=0; i<Q.size(); i++)
  {
    cout<<Q[i]<<" ";
  }
  cout<<endl;
}

void MuonPU_checkCorrection::Loop()
{
  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  int index;
  double DeltaR=0, DeltaR_val=0;
  int cpt_FLFR_bad=0, cpt_FLFR_bad_id=0, cpt_LR_bad=0, cpt_center_bad=0;

  TFile* SaveData = new TFile("ROOT_SVG/MuonPU_checkCorrection_RecoFausseSat.root", "RECREATE");

  TH2F *IhTRUE_VS_IhCORR = new TH2F("IhTRUE_VS_IhCORR","IhTRUE_VS_IhCORR",300,0,15,200,0.9,1.1);
  TH2F *IhTRUE_VS_IhOLDCORR = new TH2F("IhTRUE_VS_IhOLDCORR","IhTRUE_VS_IhOLDCORR",300,0,15,200,0.9,1.1);
  TH2F *IhTRUE_VS_IhSAT = new TH2F("IhTRUE_VS_IhSAT","IhTRUE_VS_IhSAT",300,0,15,200,0.5,1.5);
  IhTRUE_VS_IhCORR->GetXaxis()->SetTitle("I_{h}^{true} [MeV/cm]");
  IhTRUE_VS_IhCORR->GetYaxis()->SetTitle("I_{h}^{corr new} / I_{h}^{true}");
  IhTRUE_VS_IhOLDCORR->GetXaxis()->SetTitle("I_{h}^{true} [MeV/cm]");
  IhTRUE_VS_IhOLDCORR->GetYaxis()->SetTitle("I_{h}^{corr old} / I_{h}^{true}");
  IhTRUE_VS_IhSAT->GetXaxis()->SetTitle("I_{h}^{true} [MeV/cm]");
  IhTRUE_VS_IhSAT->GetYaxis()->SetTitle("I_{h}^{sat} / I_{h}^{true}");

  TH2F *BarycenterCheck_newcorr = new TH2F("BarycenterCheck_newcorr","BarycenterCheck_newcorr",200,-5,5,200,-5,5);
  BarycenterCheck_newcorr->GetXaxis()->SetTitle("bar(Q_{true})");
  BarycenterCheck_newcorr->GetYaxis()->SetTitle("bar(Q_{corr new})");

  TH1F* QCorrNew_VS_Qtrue = new TH1F("QCorrNew_VS_Qtrue","QCorrNew_VS_Qtrue",250,-0.2,2);
  QCorrNew_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr new} / Q^{true} with Max#in[240,250]");
  TH1F* QOldCorr_VS_Qtrue = new TH1F("QOldCorr_VS_Qtrue","QOldCorr_VS_Qtrue",250,-0.2,2);
  QOldCorr_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr old} / Q^{true} with Max#in[240,250]");
  TH1F* QOldCorrNC_VS_Qtrue = new TH1F("QOldCorrNC_VS_Qtrue","QOldCorrNC_VS_Qtrue",250,-0.2,2);
  QOldCorrNC_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr old NC} / Q^{true} with Max#in[240,250]");


  //nentries=100000;
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
    if (track_nvalidhits[index]>8 && track_qual[index]) selection = true;
    {
      float Ih_true = 0, Ih_sat = 0, Ih_corr = 0, Ih_oldcorr = 0;
      int cpt_isstrip = 0;

      for (int idedx=track_index_hit[index]; idedx<track_index_hit[index]+track_nhits[index]; idedx++)
      {
        if (dedx_isstrip[idedx])
        {
          cpt_isstrip++;

              // LAYER
          int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

              // SETUP
          std::vector <int> ClusterCharge;
          std::vector <int> ClusterChargeSAT;
          std::vector <int> ClusterCorr;
          int cpt_sat = 0;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
          int Nleft = 0, Nright = 0;
          int Qcorr = -1, Qoldcorr = -1, Qoldcorr_neighbourschanged = -1;

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              ClusterCharge.push_back(sigdigi_strip_adc[isigdigi]);
          
              /*if (sigdigi_strip_adc[isigdigi]<=253) ClusterChargeSAT.push_back(sigdigi_strip_adc[isigdigi]);
              else if (sigdigi_strip_adc[isigdigi]>253 && sigdigi_strip_adc[isigdigi]<=1023) ClusterChargeSAT.push_back(254);
              else if (sigdigi_strip_adc[isigdigi]>1023) ClusterChargeSAT.push_back(255);

              if (sigdigi_strip_adc[isigdigi]>=240 && sigdigi_strip_adc[isigdigi]<=250) cpt_sat++;*/
            }
          }
          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            ClusterChargeSAT.push_back(strip_ampl[istrip]);
            if (strip_ampl[istrip]>=240 && strip_ampl[istrip]<=250) cpt_sat++;
          }

          
          if (!ClusterCharge.empty() && !ClusterChargeSAT.empty())
          {
            int value_max = *max_element(ClusterChargeSAT.begin(), ClusterChargeSAT.end());
            int i_max = getIndex(ClusterChargeSAT,value_max);
            int sum_true = accumulate(ClusterCharge.begin(), ClusterCharge.end(), 0);
            int sum_sat = accumulate(ClusterChargeSAT.begin(), ClusterChargeSAT.end(), 0);
            if (ClusterChargeSAT.size()>=3) {Nleft = ClusterChargeSAT[i_max-1]; Nright = ClusterChargeSAT[i_max+1];}
            bool value_max_InRange = false;
            if (value_max >=240 && value_max <=250) value_max_InRange = true; 

                // SHAPE
            if (value_max_InRange && ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft>1.1*Nright) left=true;
            if (value_max_InRange && ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft<0.9*Nright) right=true;
            if (value_max_InRange && ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
            if (value_max_InRange && ClusterChargeSAT.size()>=2 && cpt_sat==1 && i_max==0) FullLeft=true;
            if (value_max_InRange && ClusterChargeSAT.size()>=2 && cpt_sat==1 && i_max==ClusterChargeSAT.size()-1) FullRight=true;


                // OLD CORRECTION
            if (center || left || right || FullLeft || FullRight)
            {
              Qoldcorr = Correction_OLD(ClusterChargeSAT);
              Qoldcorr_neighbourschanged = Correction_OLD_NeighboursChanged(ClusterChargeSAT);
            }

            if (Qoldcorr >= 0) QOldCorr_VS_Qtrue->Fill(1.*Qoldcorr/sum_sat);
            if (Qoldcorr_neighbourschanged >= 0) QOldCorrNC_VS_Qtrue->Fill(1.*Qoldcorr_neighbourschanged/sum_sat);
            
            if (Qoldcorr >= 0) Ih_oldcorr += 1./Power(Qoldcorr*1./dedx_pathlength[idedx],2);
            else Ih_oldcorr += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);

            
                // NEW CORRECTION
            if (center) Qcorr = Correction_wNeighbours(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
            if (left || right) Qcorr = Correction_wNeighbours(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
            if (FullLeft || FullRight) Qcorr = Correction_FL_FR_xtalk(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");
            
            int Qtemp = Qcorr;
            if (Qcorr >= 0)
            {
              for (int i=0; i<ClusterChargeSAT.size(); i++)
              {
                if (i != i_max)
                {
                  ClusterCorr.push_back(ClusterChargeSAT[i]);
                  Qcorr += ClusterChargeSAT[i];
                }
                else ClusterCorr.push_back(Qtemp);
              }

              BarycenterCheck_newcorr->Fill(ClusterBarycenter(ClusterChargeSAT), ClusterBarycenter(ClusterCorr));
              Ih_corr += 1./Power(Qcorr*1./dedx_pathlength[idedx],2);
              QCorrNew_VS_Qtrue->Fill(1.* Qcorr/sum_sat);
            }

            if (Qcorr < 240 || (!left && !right && !center && !FullLeft && !FullRight)) Ih_corr += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);
            Ih_true += 1./Power(sum_true*1./dedx_pathlength[idedx],2);
            Ih_sat += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);
            
          } // end check ClusterCharge.empty()
        } // end check dedx_isstrip
      } // end dedx loop

      IhTRUE_VS_IhCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_corr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)));
      IhTRUE_VS_IhOLDCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_oldcorr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)));
      IhTRUE_VS_IhSAT->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true)));

    } // end track loop

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  }

  SaveData->cd();

  QCorrNew_VS_Qtrue->Write();
  QOldCorr_VS_Qtrue->Write();
  QOldCorrNC_VS_Qtrue->Write();

  IhTRUE_VS_IhCORR->Write();
  IhTRUE_VS_IhOLDCORR->Write();
  //IhTRUE_VS_IhSAT->Write();

  BarycenterCheck_newcorr->Write();

  SaveData->Close();

}