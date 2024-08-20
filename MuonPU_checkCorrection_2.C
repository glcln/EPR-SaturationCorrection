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
    if (Q[i]>253) cpt_sat++;
  }

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

	if (Q.size()<=2 || cpt_sat!=1 || *(max_Q)<254 || i_max==0 || i_max==Q.size()-1)
    {
        cout << "Problem correction on LEFT, RIGHT or CENTER clusters" << endl;
        return -1;
    }
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

int Correction_OLD(const std::vector <int>&  Q)
{
  //int QII = accumulate(Q.begin(), Q.end(), 0);
  int QII = -1;
  float thresholdSat=25;

  if(Q.size()<2 || Q.size()>8) return QII;

  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
  
  if (*mQ>253)
  {
    if(*mQ==255 && (*(mQ-1)>253 || *(mQ+1)>253)) return QII;

    if (*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 && abs(*(mQ-1) - *(mQ+1))<40)
    {
      QII = (10*(*(mQ-1)) + 10*(*(mQ+1))) / 2;
      return QII;
    }
  }

  return QII; // no saturation --> no x-talk inversion
}

std::vector<int> SaturationCorrection(const std::vector <int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat, std::vector <int>& cpt_cout)
{
  std::vector<int> QII;
  std::vector<float> QI(Q.size(),0);

  if(Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    cpt_cout[1]++;
    return QII;
  }

  if(way)
  {
    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
    if(*mQ>253)
    {
      if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253)
      {
        cpt_cout[2]++;
        return Q;
      }
      if(*(mQ-1)>thresholdSat) 
      {
        if (*(mQ+1)>thresholdSat)
        {
          if (*(mQ-1)<254)
          {
            if (*(mQ+1)<254)
            {
              if (abs(*(mQ-1) - *(mQ+1))<40)
              {
                QII.push_back((10*(*(mQ-1)) + 10*(*(mQ+1)))/2);
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
      else cpt_cout[3]++;

    }
    else return Q; // no saturation --> no x-talk inversion
  }
  
  cpt_cout[0]++; //total of 3->7
  return Q;
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

int ConsecutiveSaturatedStrips (const std::vector <int> & Q)
{
  int cpt_sat = 0;
  int n = 0, warning = 0;
  std::vector <int> satInCluster;
  for (int i=0; i<Q.size(); i++)
  {
    if (warning>0)
    {
      warning = 0;
      continue;
    }
    warning = 0;
    cpt_sat = 0;
    n = 1;

    if (Q[i]>=254)
    {
      cpt_sat++;
      while (Q[i+n]>=254)
      {
        cpt_sat++;
        n++;
      }
    }
    satInCluster.push_back(cpt_sat);
    
    if (Q[i]>=254 && Q[i+1]>=254) warning++;
  }

  cpt_sat = *max_element(satInCluster.begin(), satInCluster.end());

  return cpt_sat;
}

void MuonPU_checkCorrection::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  int index;
  double DeltaR=0, DeltaR_val=0;
  int cpt_FLFR_bad=0, cpt_FLFR_bad_id=0, cpt_LR_bad=0, cpt_center_bad=0;

  TFile* SaveData = new TFile("ROOT_SVG/MuonPU_checkCorrection_Reco_MuonPU.root", "RECREATE");


  TH2F *IhTRUE_VS_IhCORR = new TH2F("IhTRUE_VS_IhCORR","IhTRUE_VS_IhCORR",300,0,15,200,0.9,1.1);
  TH2F *IhTRUE_VS_IhOLDCORR = new TH2F("IhTRUE_VS_IhOLDCORR","IhTRUE_VS_IhOLDCORR",300,0,15,200,0.9,1.1);
  TH2F *IhTRUE_VS_IhSAT = new TH2F("IhTRUE_VS_IhSAT","IhTRUE_VS_IhSAT",300,0,15,200,0.5,1.5);
  IhTRUE_VS_IhCORR->GetXaxis()->SetTitle("I_{h}^{sat} [MeV/cm]");
  IhTRUE_VS_IhCORR->GetYaxis()->SetTitle("I_{h}^{corr new} / I_{h}^{sat}");
  IhTRUE_VS_IhOLDCORR->GetXaxis()->SetTitle("I_{h}^{sat} [MeV/cm]");
  IhTRUE_VS_IhOLDCORR->GetYaxis()->SetTitle("I_{h}^{corr old} / I_{h}^{sat}");
  IhTRUE_VS_IhSAT->GetXaxis()->SetTitle("I_{h}^{true} [MeV/cm]");
  IhTRUE_VS_IhSAT->GetYaxis()->SetTitle("I_{h}^{sat} / I_{h}^{true}");

  TH2F *BarycenterCheck_newcorr = new TH2F("BarycenterCheck_newcorr","BarycenterCheck_newcorr",200,0,15,200,-5,5);
  TH2F *BarycenterCheck_sat = new TH2F("BarycenterCheck_sat","BarycenterCheck_sat",200,0,15,200,-5,5);
  BarycenterCheck_newcorr->GetXaxis()->SetTitle("bar(Q_{sat})");
  BarycenterCheck_newcorr->GetYaxis()->SetTitle("bar(Q_{corr new}) - bar(Q_{sat})");
  BarycenterCheck_sat->GetXaxis()->SetTitle("bar(Q_{true})");
  BarycenterCheck_sat->GetYaxis()->SetTitle("bar(Q_{sat}) - bar(Q_{true})");
  TH1F *BarycenterCheck_overcorrected = new TH1F("BarycenterCheck_overcorrected","BarycenterCheck_overcorrected",200,-5,5);
  BarycenterCheck_overcorrected->GetXaxis()->SetTitle("bar(Q_{corr new}) - bar(Q_{true})");
  BarycenterCheck_overcorrected->SetTitle("Cluster overcorrected 2.5 times");

  TH1F* MaxCorr_VS_MaxSat = new TH1F("MaxCorr_VS_MaxSat","MaxCorr_VS_MaxSat",250,-0.2,2);
  MaxCorr_VS_MaxSat->GetXaxis()->SetTitle("Max^{corr new} / Max^{sat}");
  TH1F* MaxOldCorr_VS_MaxSat = new TH1F("MaxOldCorr_VS_MaxSat","MaxOldCorr_VS_MaxSat",250,-0.2,2);
  MaxOldCorr_VS_MaxSat->GetXaxis()->SetTitle("Q^{corr old} / Q^{sat}");

  TH1F* ChargeDistributionReco = new TH1F("ChargeDistributionReco","ChargeDistributionReco",250,0,10);
  ChargeDistributionReco->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  TH1F* ChargeDistributionSigdigi = new TH1F("ChargeDistributionSigdigi","ChargeDistributionSigdigi",250,0,10);
  ChargeDistributionSigdigi->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  TH1F* ChargeDistributionSigdigi_zerosupress = new TH1F("ChargeDistributionSigdigi_zerosupress","ChargeDistributionSigdigi_zerosupress",250,0,10);
  ChargeDistributionSigdigi_zerosupress->GetXaxis()->SetTitle("dE/dx [MeV/cm]");

  TH2F* dEdx_sigdigi_VS_recoPU = new TH2F("dEdx_sigdigi_VS_recoPU","dEdx_sigdigi_VS_recoPU",250,0,10,250,0,10);
  dEdx_sigdigi_VS_recoPU->GetXaxis()->SetTitle("dE/dx_{reco} [MeV/cm]");
  dEdx_sigdigi_VS_recoPU->GetYaxis()->SetTitle("dE/dx_{sigdigi} [MeV/cm]");
  
  int cpt_cluster_strip = 0;
  int cpt_sat_1=0, cpt_sat_2=0, cpt_sat_3=0, cpt_sat_4=0, cpt_sat_5=0, cpt_sat_more5=0;
  

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
    if (track_nvalidhits[index]>8 && track_qual[index])
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
          std::vector <int> ClusterSigDigi_zerosupress;
          int cpt_sat = 0;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
          int Nleft = 0, Nright = 0;
          int Qcorr = -1;
          int Qoldcorr = -1;

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              ClusterCharge.push_back(sigdigi_strip_adc[isigdigi]);
          
              if (sigdigi_strip_adc[isigdigi]>=8) ClusterSigDigi_zerosupress.push_back(sigdigi_strip_adc[isigdigi]);
              /*if (sigdigi_strip_adc[isigdigi]<=253) ClusterChargeSAT.push_back(sigdigi_strip_adc[isigdigi]);
              else if (sigdigi_strip_adc[isigdigi]>253 && sigdigi_strip_adc[isigdigi]<=1023) ClusterChargeSAT.push_back(254);
              else if (sigdigi_strip_adc[isigdigi]>1023) ClusterChargeSAT.push_back(255);

              if (sigdigi_strip_adc[isigdigi]>=254) cpt_sat++;*/
            }
          }
          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            ClusterChargeSAT.push_back(strip_ampl[istrip]);
            if (strip_ampl[istrip]>=254) cpt_sat++;
          }

          
          if (!ClusterCharge.empty() && !ClusterChargeSAT.empty())
          {
            int value_max = *max_element(ClusterChargeSAT.begin(), ClusterChargeSAT.end());
            int i_max = getIndex(ClusterChargeSAT,value_max);
            int sum_true = accumulate(ClusterCharge.begin(), ClusterCharge.end(), 0);
            int sum_sat = accumulate(ClusterChargeSAT.begin(), ClusterChargeSAT.end(), 0);
            int sum_sigdigisupress = accumulate(ClusterSigDigi_zerosupress.begin(), ClusterSigDigi_zerosupress.end(), 0);
            if (ClusterChargeSAT.size()>=3) {Nleft = ClusterChargeSAT[i_max-1]; Nright = ClusterChargeSAT[i_max+1];}

            if (cpt_sat == 0)
            {
              ChargeDistributionReco->Fill(sum_sat*1./dedx_pathlength[idedx] *3.61*Power(10,-6)*265);
              ChargeDistributionSigdigi->Fill(sum_true*1./dedx_pathlength[idedx] *3.61*Power(10,-6)*265);
              ChargeDistributionSigdigi_zerosupress->Fill(sum_sigdigisupress*1./dedx_pathlength[idedx] *3.61*Power(10,-6)*265);
              dEdx_sigdigi_VS_recoPU->Fill(sum_sat*1./dedx_pathlength[idedx] *3.61*Power(10,-6)*265 , sum_true*1./dedx_pathlength[idedx] *3.61*Power(10,-6)*265);
            }

                // SHAPE
            if (ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft>1.1*Nright) left=true;
            if (ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft<0.9*Nright) right=true;
            if (ClusterChargeSAT.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterChargeSAT.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
            if (ClusterChargeSAT.size()>=2 && cpt_sat==1 && i_max==0) FullLeft=true;
            if (ClusterChargeSAT.size()>=2 && cpt_sat==1 && i_max==ClusterChargeSAT.size()-1) FullRight=true;

               // Saturated clusters stats
            int cpt_sat_insidecluster = ConsecutiveSaturatedStrips(ClusterChargeSAT);
            cpt_cluster_strip++;
            if (cpt_sat_insidecluster==1) cpt_sat_1++;
            if (cpt_sat_insidecluster==2) cpt_sat_2++;
            if (cpt_sat_insidecluster==3) cpt_sat_3++;
            if (cpt_sat_insidecluster==4) cpt_sat_4++;
            if (cpt_sat_insidecluster==5) cpt_sat_5++;
            if (cpt_sat_insidecluster>5) cpt_sat_more5++;

                // CORRECTION
            Qoldcorr = Correction_OLD(ClusterChargeSAT);
            if (Qoldcorr > 0) MaxOldCorr_VS_MaxSat->Fill(1. * Qoldcorr/sum_sat);

            if (Qoldcorr > 0) Ih_oldcorr += 1./Power(Qoldcorr*1./dedx_pathlength[idedx],2);
            else Ih_oldcorr += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);
            
            if (center) Qcorr = Correction_wNeighbours(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
            if (left || right) Qcorr = Correction_wNeighbours(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
            if (FullLeft || FullRight) Qcorr = Correction_FL_FR_xtalk(ClusterChargeSAT, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");
            
            if (Qcorr >= 254)
            {
              for (int i=0; i<ClusterChargeSAT.size(); i++)
              {
                if (i != i_max) ClusterCorr.push_back(ClusterChargeSAT[i]);
                else ClusterCorr.push_back(Qcorr);
              }
            }
                // BARYCENTER
            if (Qcorr >= 254) BarycenterCheck_newcorr->Fill(ClusterBarycenter(ClusterCorr), ClusterBarycenter(ClusterCorr) - ClusterBarycenter(ClusterChargeSAT));
            BarycenterCheck_sat->Fill(ClusterBarycenter(ClusterCharge), ClusterBarycenter(ClusterChargeSAT) - ClusterBarycenter(ClusterCharge));
            
            if (Qcorr > 0) MaxCorr_VS_MaxSat->Fill(1. * Qcorr/value_max);
            if (Qcorr >= 254)
            {
              for (int i=0; i<ClusterChargeSAT.size(); i++)
              {
                if (i != i_max) Qcorr += ClusterChargeSAT[i];
              }
              Ih_corr += 1./Power(Qcorr*1./dedx_pathlength[idedx],2);
						}
           

            if (Qcorr < 254 || (!left && !right && !center && !FullLeft && !FullRight)) Ih_corr += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);
            Ih_true += 1./Power(sum_true*1./dedx_pathlength[idedx],2);
            Ih_sat += 1./Power(sum_sat*1./dedx_pathlength[idedx],2);


            if (Qcorr >= 2.5*sum_true) BarycenterCheck_overcorrected->Fill(ClusterBarycenter(ClusterCorr) - ClusterBarycenter(ClusterCharge));
            

          } // end check ClusterCharge.empty()
        } // end check dedx_isstrip
      } // end dedx loop

      //IhTRUE_VS_IhCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_corr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true)));
      //IhTRUE_VS_IhOLDCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_oldcorr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true)));
      IhTRUE_VS_IhCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_corr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)));
      IhTRUE_VS_IhOLDCORR->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_oldcorr)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)));
      IhTRUE_VS_IhSAT->Fill(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true), 3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_sat)/(3.61*Power(10,-6)*265*Sqrt(cpt_isstrip/Ih_true)));

    } // end track loop

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  }

  SaveData->cd();

	MaxCorr_VS_MaxSat->Write();
  MaxOldCorr_VS_MaxSat->Write();
  IhTRUE_VS_IhCORR->Write();
  IhTRUE_VS_IhOLDCORR->Write();
  //IhTRUE_VS_IhSAT->Write();

  BarycenterCheck_newcorr->Write();
  ChargeDistributionReco->Write();
  ChargeDistributionSigdigi->Write();
  ChargeDistributionSigdigi_zerosupress->Write();
  dEdx_sigdigi_VS_recoPU->Write();
  //BarycenterCheck_sat->Write();
  //BarycenterCheck_overcorrected->Write();

  SaveData->Close();

  cout << " --------- Saturated clusters --------- " << endl;
  int cpt = cpt_sat_1 + cpt_sat_2 + cpt_sat_3 + cpt_sat_4 + cpt_sat_5 + cpt_sat_more5;
  cout << "Total cluster in strip: " << cpt_cluster_strip << endl;
  cout << "Total saturated clusters: " << cpt << endl;
  cout << "1 saturated strip: " << cpt_sat_1 << " / " << cpt << " (" << 100*(float)cpt_sat_1/(float)cpt << " %)" << endl;
  cout << "2 consecutive saturated strips: " << cpt_sat_2 << " / " << cpt << " (" << 100*(float)cpt_sat_2/(float)cpt << " %)" << endl;
  cout << "3 consecutive saturated strips: " << cpt_sat_3 << " / " << cpt << " (" << 100*(float)cpt_sat_3/(float)cpt << " %)" << endl;
  cout << "4 consecutive saturated strips: " << cpt_sat_4 << " / " << cpt << " (" << 100*(float)cpt_sat_4/(float)cpt << " %)" << endl;
  cout << "5 consecutive saturated strips: " << cpt_sat_5 << " / " << cpt << " (" << 100*(float)cpt_sat_5/(float)cpt << " %)" << endl;
  cout << "More than 5 consecutive saturated strips: " << cpt_sat_more5 << " / " << cpt << " (" << 100*(float)cpt_sat_more5/(float)cpt << " %)" << endl;
}