/*
    -----------------------------------------------

    Code de avec la fausse saturation ne prenant
    en compte que les clusters tq Max in [240,250]

    LES FONCTIONS SONT MODIFIEES POUR, ATTENTION

    -----------------------------------------------
*/


#define CorrectionOnData_cxx
#include "CorrectionOnData.h"
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
  //if (Q.size()<2 || Q.size()>8 || *(mQ-1) <= 15 || *(mQ+1) <= 15
  //    || *(mQ-1)>=240 || *(mQ+1)>=240 || abs(*(mQ-1) - *(mQ+1))>=40) return -1;
  if (Q.size()<2 || Q.size()>8 || *(mQ-1) <= 15 || *(mQ+1) <= 15 || *(mQ-1)>=240 || *(mQ+1)>=240) return -1;
  
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
  //if (Q.size()<2 || Q.size()>8 || *(max_Q-1) <= 15 || *(max_Q+1) <= 15
  //    || *(max_Q-1)>=240 || *(max_Q+1)>=240 || abs(*(max_Q-1) - *(max_Q+1))>=40) return -1;
  if (Q.size()<2 || Q.size()>8 || *(max_Q-1) <= 15 || *(max_Q+1) <= 15 || *(max_Q-1)>=240 || *(max_Q+1)>=240) return -1;

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
    if (*(mQ-1)>15 && *(mQ+1)>15 && *(mQ-1)<240 && *(mQ+1)<240)// && abs(*(mQ-1) - *(mQ+1))<40)
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
    if (*(mQ-1)>15 && *(mQ+1)>15 && *(mQ-1)<240 && *(mQ+1)<240)// && abs(*(mQ-1) - *(mQ+1))<40)
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

void CorrectionOnData::Loop()
{
  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  int cpt_OldCorr_under = 0, cpt_OldCorrNC = 0, cpt_CorrNew_under = 0;
  TFile* SaveData = new TFile("ROOT_SVG/OnData_Vpm_applied.root", "RECREATE");

  TH1F* QNoCorr_VS_Qtrue = new TH1F("QNoCorr_VS_Qtrue","QNoCorr_VS_Qtrue",250,-0.2,2);
  QNoCorr_VS_Qtrue->GetXaxis()->SetTitle("Q^{no corr} / Q^{true} with Max#in[240,250]");
  TH1F* QOldCorr_VS_Qtrue = new TH1F("QOldCorr_VS_Qtrue","QOldCorr_VS_Qtrue",250,-0.2,2);
  QOldCorr_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr old} / Q^{true} with Max#in[240,250]");
  TH1F* QOldCorrNC_VS_Qtrue = new TH1F("QOldCorrNC_VS_Qtrue","QOldCorrNC_VS_Qtrue",250,-0.2,2);
  QOldCorrNC_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr old NC} / Q^{true} with Max#in[240,250]");
  TH1F* QCorrNew_VS_Qtrue = new TH1F("QCorrNew_VS_Qtrue","QCorrNew_VS_Qtrue",250,-0.2,2);
  QCorrNew_VS_Qtrue->GetXaxis()->SetTitle("Q^{corr new} / Q^{true} with Max#in[240,250]");
  
  TH2F *BarycenterCheck_newcorr = new TH2F("BarycenterCheck_newcorr","BarycenterCheck_newcorr",200,-5,5,200,-5,5);
  BarycenterCheck_newcorr->GetXaxis()->SetTitle("bar(Q_{true})");
  BarycenterCheck_newcorr->GetYaxis()->SetTitle("bar(Q_{corr new})");


  cout << "Histograms saved in: " << SaveData->GetName() << endl;
  //nentries=100000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

        // TRACK LOOP
    for (int itr=0; itr<ntracks; itr++)
    {
      bool selection = true;

      // all except HLT 
      if (track_pt[itr] < 30 || track_pt[itr] > 55) selection=false;  // NEW
      if (abs(track_eta[itr]) >= 2.4) selection=false;                // NEW
      if (track_npixhits[itr] < 2) selection=false;
      if (track_validfraction[itr] <= 0.8) selection=false;
      if (track_nvalidhits[itr] < 10) selection=false;
      if (!track_qual[itr]) selection=false;
      if (track_chi2[itr] >= 5) selection=false;
      if (abs(track_dz[itr]) >= 0.1) selection=false;
      if (abs(track_dxy[itr]) >= 0.02) selection=false;
      //if (track_miniRelIso[itr] >= 0.02) selection=false;           // not in the data tree
      //if (track_TkRelIso[itr] >= 15) selection=false;               // not in the data tree
      //if (track_EoP[itr] >= 0.3) selection=false;                   // not in the data tree
      if (track_pterr[itr]/(track_pt[itr]*track_pt[itr]) >= 0.0008) selection=false;
      if (1 - track_probQNoL1[itr] <= 0.3) selection=false;

      if (selection)
      {
        float Ih_true = 0, Ih_sat = 0, Ih_corr = 0, Ih_oldcorr = 0;
        int cpt_isstrip = 0;

        for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
        {
          if (dedx_isstrip[idedx])
          {
            cpt_isstrip++;

                // LAYER
            int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

                // SETUP
            std::vector <int> ClusterChargeSAT;
            std::vector <int> ClusterCorr;
            int cpt_sat = 0;
            bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
            int Nleft = 0, Nright = 0;
            int Qcorr = -1, Qoldcorr = -1, Qoldcorr_NeighboursChanged = -1;

            for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
            {
              ClusterChargeSAT.push_back(strip_ampl[istrip]);
              if (strip_ampl[istrip]>=240 && strip_ampl[istrip]<=250) cpt_sat++;
            }

            
            if (!ClusterChargeSAT.empty())
            {
              int value_max = *max_element(ClusterChargeSAT.begin(), ClusterChargeSAT.end());
              int i_max = getIndex(ClusterChargeSAT,value_max);
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

                  // NO CORRECTION
              if (value_max_InRange) 
              {
                int QNoCorr = 0;
                for (int i=0; i<ClusterChargeSAT.size(); i++)
                {
                  if (i != i_max) QNoCorr += ClusterChargeSAT[i];
                  else QNoCorr += 240;
                }
                QNoCorr_VS_Qtrue->Fill(1.*QNoCorr/sum_sat);
              }

                  // OLD CORRECTION
              if (center || left || right || FullLeft || FullRight)
              {
                Qoldcorr = Correction_OLD(ClusterChargeSAT);
                Qoldcorr_NeighboursChanged = Correction_OLD_NeighboursChanged(ClusterChargeSAT);
              }

              if (Qoldcorr >= 0) QOldCorr_VS_Qtrue->Fill(1.*Qoldcorr/sum_sat);
              if (Qoldcorr_NeighboursChanged >= 0) QOldCorrNC_VS_Qtrue->Fill(1.*Qoldcorr_NeighboursChanged/sum_sat);
              
              if (Qoldcorr >= 0 && Qoldcorr<240) cpt_OldCorr_under++;
              if (Qoldcorr_NeighboursChanged >= 0 && Qoldcorr_NeighboursChanged<240) cpt_OldCorrNC++;

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
                QCorrNew_VS_Qtrue->Fill(1.* Qcorr/sum_sat);
              }

              if (Qcorr >= 0 && Qcorr<240) cpt_CorrNew_under++;          
              
            } // end check ClusterChargeSAT.empty()
          } // end check dedx_isstrip
        } // end dedx loop
      } // end selection on track
    } // end track loop

    if (jentry%100000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  }

  SaveData->cd();

  QNoCorr_VS_Qtrue->Write();
  QOldCorr_VS_Qtrue->Write();
  QOldCorrNC_VS_Qtrue->Write();
  QCorrNew_VS_Qtrue->Write();

  BarycenterCheck_newcorr->Write();

  SaveData->Close();


  cout << endl;
  cout << "       BAD CLUSTER CORRECTION" << endl;
  cout << "Old corr, Qcorr < 240: " << cpt_OldCorr_under << endl;
  cout << "Old corr NC, Qcorr < 240: " << cpt_OldCorrNC << endl;
  cout << "New corr, Qcorr < 240: " << cpt_CorrNew_under << endl;
  cout << endl;

}