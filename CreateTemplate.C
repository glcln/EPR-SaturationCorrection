#define CreateTemplate_cxx
#include "CreateTemplate.h"
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
  TH1F *hist=new TH1F(Form("Cluster_%d_%d_%d__%4.3f-%4.3f",entry,itr,layer,VmaxOverVmin,MaxOverVmax),Form("Cluster_%d_%d_%d__%4.3f-%4.3f",entry,itr,layer,VmaxOverVmin,MaxOverVmax),C.size()+2,0,C.size()+2);
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

double getMaxHisto_X(TH1D* h)
{
  return  h->GetMaximumBin() * (h->GetXaxis()->GetBinCenter(h->GetNbinsX()) - h->GetXaxis()->GetBinCenter(1)) * 1./h->GetNbinsX();
}

void CreateTemplate::Loop()
{
  gROOT->SetBatch(kTRUE);

  TFile* SaveData = new TFile("ROOT_SVG/Check_Template.root", "RECREATE");

  ofstream Template_FLFR("ROOT_SVG/Template_correction/Template_FLFR.txt", std::ofstream::out);
  ofstream Template_CENTER("ROOT_SVG/Template_correction/Template_CENTER.txt", std::ofstream::out);
  ofstream Template_FirstNeighbour("ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", std::ofstream::out);
  
  
  // Template First Neighbour at a coarse cut
  ifstream Template_FirstNeighbour_cut("ROOT_SVG/Template_correction/Template_LEFTRIGHT_firstcut.txt");
  std::string line;
  std::vector<double> template_a1;
  std::vector<double> template_a2;
  while (std::getline(Template_FirstNeighbour_cut, line))
  {
    std::istringstream iss(line);
    int i_layer=-1;
    double a1=-1, a2=-1;
    if (!(iss >> i_layer >> a1 >> a2 )) break; // error
      
    template_a1.push_back(a1);
    template_a2.push_back(a2);
  }
  Template_FirstNeighbour_cut.close();


  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  std::vector <TH1D*> FL_TH2;
  std::vector <TH1D*> FR_TH2;
  std::vector <TH2D*> LR_TH2;
  std::vector <TH1D*> Center_TH2;
  std::vector <TH2D*> LR_254_TH2;
  std::vector <TH2D*> LR_255_TH2;
  double DeltaR=0;
  double DeltaR_val;
  int index;

  for (int i=0; i<=20; i++)
  {
    TString histName = Form("FL_a_%d", i);
    TH1D *histo = new TH1D(histName,histName,200,0,20);
    histo->GetXaxis()->SetTitle("Max/Vmax");
    FL_TH2.push_back(histo);

    TString histName2 = Form("FR_a_%d", i);
    TH1D *histo2 = new TH1D(histName2,histName2,200,0,20);
    histo2->GetXaxis()->SetTitle("Max/Vmax");
    FR_TH2.push_back(histo2);

    TString histName3 = Form("LR_%d", i);
    TH2D *histo3 = new TH2D(histName3,histName3,100,0,10,200,0,20);
    histo3->GetXaxis()->SetTitle("Vmax/Vmin");
    histo3->GetYaxis()->SetTitle("Max/Vmax");
    LR_TH2.push_back(histo3);

    TString histName4 = Form("Center_%d", i);
    TH1D *histo4 = new TH1D(histName4,histName4,200,0,20);
    histo4->GetXaxis()->SetTitle("Max/Vmax");
    Center_TH2.push_back(histo4);

    TString histName5 = Form("LR_254_%d", i);
    TH2D *histo5 = new TH2D(histName5,histName5,100,0,10,200,0,20);
    histo5->GetXaxis()->SetTitle("Vmax/Vmin");
    histo5->GetYaxis()->SetTitle("Max/Vmax");
    LR_254_TH2.push_back(histo5);

    TString histName6 = Form("LR_255_%d", i);
    TH2D *histo6 = new TH2D(histName6,histName6,100,0,10,200,0,20);
    histo6->GetXaxis()->SetTitle("Vmax/Vmin");
    histo6->GetYaxis()->SetTitle("Max/Vmax");
    LR_255_TH2.push_back(histo6);
  }

  TH2D *MaxVMax_VS_size_LR = new TH2D("MaxVMax_VS_size_LR","MaxVMax_VS_size_LR",20,0,20,200,0,20);


  //nentries = 20000;
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
          vector<int> index_sig;
          vector<int> value_sig;
          double value_max=-1;
          int i_max=-1;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
          
          int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              index_sig.push_back(isigdigi);
              value_sig.push_back(sigdigi_strip_adc[isigdigi]);
            }
          }
          if (!index_sig.empty() && !value_sig.empty())
          {
            value_max = *max_element(value_sig.begin(), value_sig.end());
            i_max = index_sig.front()+getIndex(value_sig,value_max);

            int check_OneMaxSat=0;
            for (int i=0; i<value_sig.size(); i++)
            {
              if (value_sig[i]>253) check_OneMaxSat++;
            }

            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>253 && i_max>index_sig.front()
            && i_max<index_sig.back() && sigdigi_strip_adc[i_max-1]>1.1*sigdigi_strip_adc[i_max+1]) left=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>253 && i_max>index_sig.front()
            && i_max<index_sig.back() && sigdigi_strip_adc[i_max-1]<0.9*sigdigi_strip_adc[i_max+1]) right=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>253 && i_max>index_sig.front()
            && i_max<index_sig.back() && sigdigi_strip_adc[i_max-1]<=1.1*sigdigi_strip_adc[i_max+1]
            && sigdigi_strip_adc[i_max-1]>=0.9*sigdigi_strip_adc[i_max+1]) center=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && value_max>253 && i_max==index_sig.front()) FullLeft=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && value_max>253 && i_max==index_sig.back()) FullRight=true;
            
                // FL AND FR TEMPLATES
            if (FullLeft && ( ( sigdigi_strip_adc[i_max+1]>=16 && (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==4 || (layer>=14 && layer<=17)) )
            || ( sigdigi_strip_adc[i_max+1]>=25 && (dedx_subdetid[idedx]==5 || layer>=18) ) ) )
            {
              FL_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max] / sigdigi_strip_adc[i_max+1]);  
            }

            if (FullRight && ( ( sigdigi_strip_adc[i_max-1]>=16 && (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==4 || (layer>=14 && layer<=17)) )
            || ( sigdigi_strip_adc[i_max-1]>=25 && (dedx_subdetid[idedx]==5 || layer>=18) ) ) )
            {
              FR_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max] / sigdigi_strip_adc[i_max-1]);
            }

                // LR TEMPLATES
            if (left || right)
            {

              if (layer == 2 && (sigdigi_strip_adc[i_max+1]>50 || sigdigi_strip_adc[i_max-1]>50)) MaxVMax_VS_size_LR->Fill(index_sig.size(), sigdigi_strip_adc[i_max]*1./(sigdigi_strip_adc[i_max+1]+sigdigi_strip_adc[i_max-1]));
                
              if (sigdigi_strip_adc[i_max-1]>=sigdigi_strip_adc[i_max+1])
              {
                if ((float)sigdigi_strip_adc[i_max] >= 0.85*(template_a1[layer-1]*(float)sigdigi_strip_adc[i_max-1] + template_a2[layer-1]*(float)sigdigi_strip_adc[i_max+1]) )
                LR_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);

                //if (layer == 2) MaxVMax_VS_size_LR->Fill(index_sig.size(), (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (value_max>253 && value_max<=1023) LR_254_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (value_max>1023) LR_255_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
              }
              else
              {
                if ((float)sigdigi_strip_adc[i_max] >= 0.85*(template_a1[layer-1]*(float)sigdigi_strip_adc[i_max+1] + template_a2[layer-1]*(float)sigdigi_strip_adc[i_max-1]) ) 
                LR_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
              
                //if (layer == 2) MaxVMax_VS_size_LR->Fill(index_sig.size(), (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (value_max>253 && value_max<=1023) LR_254_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (value_max>1023) LR_255_TH2[layer]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
              }
            }

                // CENTER TEMPLATES
            if (center) Center_TH2[layer]->Fill( 2.*sigdigi_strip_adc[i_max] / (sigdigi_strip_adc[i_max+1] + sigdigi_strip_adc[i_max-1]));

          } //end if index_sig
        } //end if strip
      } //end hit loop
    } //end track loop

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  } // end entry loop

      // FL AND FR TEMPLATES
  for (int i=1;i<FL_TH2.size();i++) Template_FLFR << i << " " << getMaxHisto_X(FL_TH2[i]) << " " << getMaxHisto_X(FR_TH2[i]) << "\n";
  Template_FLFR.close();
  

      // Max/Vmax TEMPLATES
  std::vector <float> p0;
  std::vector <float> p1;
  for (int i=0; i<LR_TH2.size(); i++)
  {
    TF1 *fit = new TF1("fit","[0] + [1]/x",1,10);
    LR_TH2[i]->Fit("fit","REMQ");

    p0.push_back(fit->GetParameter(0));
    p1.push_back(fit->GetParameter(1));
  }
  for (int i=1; i<p0.size(); i++) Template_FirstNeighbour << i << " " << p0[i] << " " << p1[i] << "\n";
  Template_FirstNeighbour.close();

  for (int i=1;i<Center_TH2.size();i++) Template_CENTER << i << " " << getMaxHisto_X(Center_TH2[i]) << "\n";
  Template_CENTER.close();


  SaveData->cd();
  MaxVMax_VS_size_LR->Write();
  for (int i=1; i<LR_TH2.size(); i++) LR_TH2[i]->Write();
  for (int i=1; i<Center_TH2.size(); i++) Center_TH2[i]->Write();
  for (int i=0; i<FL_TH2.size(); i++) {FR_TH2[i]->Write(); FL_TH2[i]->Write();}
  for (int i=0; i<LR_254_TH2.size(); i++) {LR_254_TH2[i]->Write(); LR_255_TH2[i]->Write();}
  for (int i=0; i<LR_TH2.size(); i++) {delete LR_TH2[i]; delete Center_TH2[i];delete FL_TH2[i]; delete FR_TH2[i];}
  for (int i=0; i<LR_254_TH2.size(); i++) {delete LR_254_TH2[i]; delete LR_255_TH2[i];}
}

