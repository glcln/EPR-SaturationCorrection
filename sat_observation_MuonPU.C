#define sat_observation_MuonPU_cxx
#include "sat_observation_MuonPU.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
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

double Correction_wNeighbours(int Vmin, int Vmax, int layer)
{
  std::vector<double> template_a = {-0.67307, -0.687187, -0.890144, -0.855589, -0.552413, -0.564462, -0.584204, -0.595558,
  -0.525664, -0.509628, -0.533439, -1.431693, -1.249519, -0.991215, -1.093429, -0.893097, -0.833976, -0.681379, -0.675934, -0.491908};
  std::vector<double> template_b = {12.68, 12.5859, 15.7077, 15.6552, 10.4165, 10.415, 10.4359, 10.4078, 8.17375, 8.146,
  14.775403, 19.473629, 19.286545, 18.035135, 19.712891, 15.955571, 17.185255, 11.332367, 11.431352, 9.636855};

  // templates index:
  //    0->3: TIB (layer)
  //    4->9: TOB (layer)
  //    10->12: TID (ring)
  //    13->19: TEC (ring)

  return template_a[layer-1]*Vmax + template_b[layer-1]*Vmin;
}

double Correction_wBarycenter(int Vmin, int i_max, double barycenter, int layer)
{
  std::vector<double> template_a = {-0.67307, -0.687187, -0.890144, -0.855589, -0.552413, -0.564462, -0.584204, -0.595558,
  -0.525664, -0.509628, -1.063325, -1.142731, -1.054355, -0.991215, -1.093429, -0.893097, -0.833976, -0.681379, -0.675934, -0.491908};
  std::vector<double> template_b = {12.68, 12.5859, 15.7077, 15.6552, 10.4165, 10.415, 10.4359, 10.4078, 8.17375, 8.146,
  17.678024, 17.680134, 17.776665, 18.035135, 19.712891, 15.955571, 17.185255, 11.332367, 11.431352, 9.636855};
  
  std::vector<double> template_c = {0.868764, 0.880221, 0.862678, 0.855892, 0.856870, 0.861826, 0.866404, 0.862338, 
  0.878096, 0.881639, 0.871144, 0.949353, 0.917536, 1.006906, 0.974363, 0.902437, 0.906889, 0.824679, 0.855054, 0.909103};
  std::vector<double> template_d = {14.641308, 14.333699, 17.128500, 17.519157, 11.296318, 11.216012, 11.038068, 11.175527, 
  8.644450, 8.602039, 16.058081, 17.335613, 18.797688, 17.049026, 17.574888, 15.430065, 18.386696, 11.838357, 12.115040, 10.681962};
  
  std::vector<double> template_cp = {0.802739, 0.792635, 0.800417, 0.856294, 0.602448, 0.616582, 0.621663, 0.615896, 
  0.738393, 0.736678, 0.946892, 1.093929, 1.349523, 1.547499, 1.267664, 0.844524, 1.140316, 0.565715, 0.650845, 0.784318};
  std::vector<double> template_dp = {-14.870478, -14.921204, -17.388657, -17.115129, -12.534080, -12.470596, -12.290698, -12.187500, 
  -9.132001, -9.221056, -17.031996, -15.913959, -15.21862, -10.304313, -16.601250, -15.731742, -17.790249, -12.963305, -12.760837, -11.228178};

  // templates index:
  //    0->3: TIB (layer)
  //    4->9: TOB (layer)
  //    10->12: TID (ring)
  //    13->19: TEC (ring)

  if (i_max-barycenter >= 0) return (template_a[layer-1]*template_c[layer-1] + template_b[layer-1])*Vmin + template_a[layer-1]*template_d[layer-1]*(i_max-barycenter)*Vmin;
  else return (template_a[layer-1]*template_cp[layer-1] + template_b[layer-1])*Vmin + template_a[layer-1]*template_dp[layer-1]*(i_max-barycenter)*Vmin;
}

double Correction_wNeighboursBarycenter(int Vmax, int i_max, double barycenter, int layer)
{
  // Exponential
  std::vector<double> template_a = { 13.277793, 13.059519, 16.484966, 16.312445, 11.343474, 11.244805, 11.235435, 11.185284,
  8.484281, 8.354471, 15.631947, 21.338182, 25.304216, 16.511023, 22.838003, 18.094494, 19.447824, 13.050611, 12.010463};
  std::vector<double> template_b = {-10.490966, -10.642560, -9.496873, -10.034931, -7.733625, -7.533386, -7.448637, -7.553052,
  -7.590789, -7.451965, -10.429406, -11.113056, -12.736379, -8.846722, -11.935238, -10.763234, -11.946908, -8.423126, -7.240906, -7.468538};
  std::vector<double> template_c = {13.277793, 13.099030, 16.774704, 16.583477, 11.326892, 11.220656, 11.235435, 11.214873,
  8.511932, 8.354471, 16.189613, 20.779060, 23.918007, 19.870031, 23.448797, 18.232422, 19.214231, 13.286536, 12.166878, 9.546692};
  std::vector<double> template_d = {10.490966, 10.657241, 9.608710, 10.117137, 7.691601, 7.575212, 7.448637, 7.595020,
  7.639862, 7.451965, 11.097799, 10.970179, 12.002810, 10.159702, 12.022460, 11.243575, 11.702321, 8.665884, 7.501427, 7.938479};

  // Absolute value
  std::vector<double> template_e = { 11.374518, 11.226698, 14.298002, 14.129365, 9.773895, 9.753185, 9.764444, 9.727222, 7.528874, 7.495228,
  13.373905, 17.337006, 17.716879, 16.434248, 17.719839, 14.487938, 15.585187, 10.576308, 10.683282, 8.943758};
  std::vector<double> template_f = {-35.670902, -35.017609, -44.983925, -44.641060, -27.419739, -27.210995, -26.723938, -26.782194, -19.670015,
  -19.439096, -40.951248, -54.639954, -56.932945, -54.235996, -56.534191, -43.377903, -49.286282, -29.325024, -29.476278, -24.483458};
  

  // templates index:
  //    0->3: TIB (layer)
  //    4->9: TOB (layer)
  //    10->12: TID (ring)
  //    13->19: TEC (ring)
  
  return (Vmax * (template_e[layer-1] + template_f[layer-1] * TMath::Abs(i_max - barycenter)) ); 
  /*if (i_max-barycenter > 0.03) return (Vmax * (template_a[layer-1] * TMath::Exp(template_b[layer-1] * (i_max-barycenter)) + 1) );
  else if (i_max-barycenter < -0.03) return (Vmax * (template_c[layer-1] * TMath::Exp(template_d[layer-1] * (i_max-barycenter)) + 1) );
  else if (i_max-barycenter <= 0.03 || i_max-barycenter >= -0.03) return (Vmax * (template_e[layer-1] + template_f[layer-1] * TMath::Abs(i_max - barycenter)) );
  else 
  {
    cout<<"ERROR correction with Correction_wNeighboursBarycenter: barycenter position "<<endl;
    return -1;
  }*/
}

std::vector<int> Correction_FL_FR_xtalk(const std::vector <int>&  Q, int layer, bool &NotCorrected, std::string file)
{
  std::vector<int> QII;
  int cpt_sat = 0;
  float thresholdSat = 25;
  
  // CORRECTION TEMPLATES
  std::string line;
  std::vector<double> template_FL;
  std::vector<double> template_FR;

  std::ifstream TemplateFILE(file);
  while (std::getline(TemplateFILE, line))
  {
    std::istringstream iss(line);
    int i_layer=-1;
    double FL_a=-1, FR_a=-1;
    if (!(iss >> i_layer >> FL_a >> FR_a)) break; // error
      
    template_FL.push_back(FL_a);
    template_FR.push_back(FR_a);
  }
  TemplateFILE.close();

  // NUMBER OF MAX
  for (unsigned int i=0;i<Q.size();i++)
  {
    if (Q[i]>253) cpt_sat++;
  }

  // NO CORRECTION IF TOO SMALL OR LARGE
  if (Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    NotCorrected = true;
    return QII;
  }

  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
  
  if (cpt_sat == 1)
  {
    // FULL LEFT
    if (*mQ==Q[0] && *(mQ+1)>thresholdSat && *(mQ+1)<254)
    {
      QII.push_back(template_FL[layer-1] * (*(mQ+1)));
      NotCorrected = false;
      return QII;
    }

    // FULL RIGHT
    if (*mQ==Q[Q.size()-1] && *(mQ-1)>thresholdSat && *(mQ-1)<254)
    {
      QII.push_back(template_FR[layer-1] * (*(mQ-1)));
      NotCorrected = false;
      return QII;
    }   
    
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    NotCorrected = true;
    return QII;
  }

  // NO SATURATION --> no x-talk inversion
  else
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    NotCorrected = true;
    return QII;
  }
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
    if(((detid>>5)&0x7)==1) return 1;
    else if(((detid>>5)&0x7)==2) return 2;
    else if(((detid>>5)&0x7)==3) return 3;
    else if(((detid>>5)&0x7)==4) return 4;
    else if(((detid>>5)&0x7)==5) return 5;
    else if(((detid>>5)&0x7)==6) return 6;
    else if(((detid>>5)&0x7)==7) return 7;
    else return -1;
  }
  return -1;
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

int getIndex(vector<int> v, int K) 
{ 
  auto it = find(v.begin(), v.end(), K); 
  int index = it - v.begin(); 

  return index; 
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

std::vector<int> CrossTalkInv(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) 
{
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  bool debugbool=false;
  TMatrix A(N,N);

//---  que pour 1 max bien net
 if(Q.size()<2 || Q.size()>8){
	for (unsigned int i=0;i<Q.size();i++){
		QII.push_back((int) Q[i]);
  	}
	return QII;
  }
 if(way){
	  std::vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())	;
	  if(*mQ>253){
	 	 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
	 	 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
		     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
	  }
  }
//---

  for(unsigned int i=0; i<N; i++) {
        A(i,i) =a;
        if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
        else continue;
        if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
  }

  if(N==1) A(0,0)=1/a;
  else  A.InvertFast();

  for(unsigned int i=0; i<N; i++) {
        for(unsigned int j=0; j<N; j++) {
        QI[i]+=A(i,j)*(float)Q[j];
        }
  }

 for (unsigned int i=0;i<QI.size();i++){
	if(QI[i]<threshold) QI[i]=0;
	QII.push_back((int) QI[i]);
  }

return QII;
}

double getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, double dropHigherDeDxValue, int & nv, int & ns)
{
  double result=-1;
//double dropLowerDeDxValue=0.15;
  size_t MaxStripNOM=99;
  bool usePixel=true;
  bool useStrip=true;

  std::vector<double> vect;

  bool debugprint=false;
  unsigned int SiStripNOM = 0;
  ns=0;

  for(unsigned int h=0;h<charge.size();h++){
    if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
    if(!usePixel && subdetId[h]<3)continue; // skip pixels
    if(!useStrip && subdetId[h]>=3)continue; // skip strips
    if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;

    if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
    if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
    if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

    int ClusterCharge = charge[h];
    if (subdetId[h]>=3 && charge[h]>=254) ns++;

    //To uncomment and delete two lines above
    double scaleFactor = 1;
    if (subdetId[h]<3) scaleFactor *= 1; // add pixel scaling
//         double scaleFactor = scaleFactors[0];
//         if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
    if (debugprint) std::cout << " after SF " << std::endl;

    if(templateHisto){  //save discriminator probability
        double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
        int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
        int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
        int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
        double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
        vect.push_back(Prob); //save probability
        if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
    }else{
        double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
        double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/pathlength[h];
        vect.push_back(ChargeOverPathlength); //save charge
        if (debugprint) std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
    }
  }

  if(dropLowerDeDxValue>0){
      std::vector <double> tmp (vect.size());
      std::copy (vect.begin(), vect.end(), tmp.begin());
      std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
      int nTrunc = tmp.size()*dropLowerDeDxValue;

      vect.clear();
      for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
  }
  if (debugprint) std::cout << " after dropLowerDeDxValue " << std::endl;

  if(dropHigherDeDxValue>0){
      std::vector <double> tmp (vect.size());
      std::copy (vect.begin(), vect.end(), tmp.begin());
      std::sort(tmp.begin(), tmp.end(), std::less<double>() );
      int nTrunc = tmp.size()*dropHigherDeDxValue;

      vect.clear();
      for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
  }
  if (debugprint) std::cout << " after dropHigherDeDxValue " << std::endl;
  int size = vect.size();
  nv = size;

  if(size>0){
    if(templateHisto){  //dEdx discriminator
      //Ias discriminator
      result = 1.0/(12*size);
        std::sort(vect.begin(), vect.end(), std::less<double>() );
        for(int i=1;i<=size;i++){
          result += vect[i-1] * pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
        }
        result *= (3.0/size);
        if (debugprint) std::cout << " Ias discriminator " << result << std::endl;
    }else{  //dEdx estimator
        //harmonic2 estimator
        result=0;
//           double expo = -2;
        double expo = -1* n_estim;
        for(int i = 0; i< size; i ++){
          if (vect[i]!=0) result+=pow(vect[i],expo);
      }
        //cout<<"done, result: "<<result<<endl;
        result = pow(result/size,1./expo);
        if (debugprint) std::cout << " harmonic discriminator " << result << " with expo " << expo << std::endl;
    }
  }else{
    result = -1;
  }
  if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;

  return result;
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

std::vector<int> New_SaturationCorrection(const std::vector <int>&  Q, bool way, float thresholdSat, std::vector <int>& cpt_cout)
{
  std::vector<int> QII;
  unsigned int cpt_sat=0;
  bool mleft=false, mright=false;

  // NUMBER OF MAX
  for (unsigned int i=0;i<Q.size();i++)
  {
    if (Q[i]>253) cpt_sat++;
  }

  if(Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    cpt_cout[0]++;
    return QII;
  }

  if(way)
  {
    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
    
    // ONE SATURATED STRIP
    if(cpt_sat==1)
    {
      if(*mQ==Q[0] && *(mQ+1)>thresholdSat && *(mQ+1)<254)  //FullLeft one sat
      {
        QII.push_back(9.26*(*(mQ+1)));
        cpt_cout[1]++;
        return QII;
        // max = 0.879 (sat)
        // 1st right neighbor = 0.108 -> inverse = 1/0.108 = 9.26
      }

      if(*mQ==Q[Q.size()-1] && *(mQ-1)>thresholdSat && *(mQ-1)<254)  //FullRight one sat
      {
        QII.push_back(9.26*(*(mQ-1)));
        cpt_cout[2]++;  
        return QII;
        // max = 0.872 (sat)
        // 1st left neighbor = 0.108 -> inverse = 9.26
      }

      if(Q.size()>2 && *(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 && (*(mQ-1)>1.1*(*(mQ+1)) || *(mQ+1)>1.1*(*(mQ-1)) || abs(*(mQ-1)/(*(mQ+1)) -1)<=0.1))
      {
        if (*(mQ-1)>1.1*(*(mQ+1)))
        {
          QII.push_back((4.78*(*(mQ-1))+13.51*(*(mQ+1)))/2); //Left one sat
          cpt_cout[3]++;
          return QII;
        }
        else if (*(mQ+1)>1.1*(*(mQ-1)))
        {
          QII.push_back((13.51*(*(mQ-1))+4.78*(*(mQ+1)))/2); //Right one sat
          cpt_cout[4]++;
          return QII;
        }
        else if (abs(*(mQ-1)/(*(mQ+1)) -1)<=0.1)
        {
          QII.push_back((11.90*(*(mQ-1))+11.90*(*(mQ+1)))/2);  //Center one sat
          cpt_cout[5]++;
          return QII;
        }
        // Left/Right: max = 0.617 (sat)
        // Left/Right: 1st left/right neighbor = 0.209 (no sat) -> inverse = 4.78
        // Left/Right: 1st right/left neighbor = 0.074  (no sat) -> inverse = 13.51
        
        // Center: max = 0.811 (sat)
        // Center: 1st neighbors = 0.084 (no sat) -> inverse = 11.90
      }
    }
    
    // TWO CONSECUTIVE SATURATED STRIP
    if (Q.size()>2 && cpt_sat==2)
    {
      //position set
      if (*mQ>253 && *(mQ+1)>253) mleft=true;
      if (*mQ>253 && *(mQ-1)>253) mright=true;

      if((mleft && *mQ==Q[0] && *(mQ+2)>thresholdSat) || (mright && *(mQ-1)==Q[0] && *(mQ+1)>thresholdSat))  //FullLeft two sat
      {
        if (mleft)
        {
          QII.push_back(0.54*10.64*(*(mQ+2)));
          QII.push_back(0.46*10.64*(*(mQ+2)));
        }
        if (mright)
        {
          QII.push_back(0.54*10.64*(*(mQ+1)));
          QII.push_back(0.46*10.64*(*(mQ+1)));
        }
        cpt_cout[6]++;
        return QII;
        // max = 0.466 (sat)
        // 1st right neighbor = 0.403 (sat)
        // 2nd right neighbor = 0.094  (no sat) -> inverse = 10.64
      }

      if((mleft && *(mQ+1)==Q[Q.size()-1] && *(mQ-1)>thresholdSat) || (mright && *mQ==Q[Q.size()-1] && *(mQ-2)>thresholdSat))  //FullRight two sat
      {
        if (mleft)
        {
          QII.push_back(0.42*13.9*(*(mQ-1)));
          QII.push_back(0.58*13.9*(*(mQ-1)));
        }
        if (mright)
        {
          QII.push_back(0.42*13.9*(*(mQ-2)));
          QII.push_back(0.58*13.9*(*(mQ-2)));
        }
        cpt_cout[7]++;
        return QII;
        // max = 0.517 (sat)
        // 1st left neighbor = 0.386 (sat)
        // 2nd left neighbor = 0.072  (no sat) -> inverse = 13.9
      }

      if(Q.size()>3 && ((mleft && *(mQ-1)>thresholdSat && *(mQ+2)>thresholdSat && *(mQ-1)<254 && *(mQ+2)<254) || (mright && *(mQ-2)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-2)<254 && *(mQ+1)<254))
         && (((mleft && *(mQ-1)>1.1*(*(mQ+2))) || (mright && *(mQ-2)>1.1*(*(mQ+1)))) 
          || ((mleft && *(mQ+2)>1.1*(*(mQ-1))) || (mright && *(mQ+1)>1.1*(*(mQ-2))))
          || ((mleft && abs(*(mQ-1)/(*(mQ+2)) -1)<=0.1) || (mright && abs(*(mQ-2)/(*(mQ+1)) -1)<=0.1))
            )) 
      {
        if (mleft && *(mQ-1)>1.1*(*(mQ+2))) //Left mleft 2 sat
        {
          QII.push_back((4.4*(*(mQ-1))+10*(*(mQ+2)))/2);
          QII.push_back((4.4*(*(mQ-1))+10*(*(mQ+2)))/2);
          cpt_cout[8]++;
        }
        else if (mright && *(mQ-2)>1.1*(*(mQ+1))) //Left mright 2 sat
        {
          QII.push_back((4.4*(*(mQ-2))+10*(*(mQ+1)))/2);
          QII.push_back((4.4*(*(mQ-2))+10*(*(mQ+1)))/2);
          cpt_cout[8]++;
        }
        else if (mleft && *(mQ+2)>1.1*(*(mQ-1))) //Right mleft 2 sat
        {
          QII.push_back((10.2*(*(mQ-1))+4.5*(*(mQ+2)))/2);
          QII.push_back((10.2*(*(mQ-1))+4.5*(*(mQ+2)))/2);
          cpt_cout[9]++;
        }
        else if (mright && *(mQ+1)>1.1*(*(mQ-2))) //Right mright 2 sat
        {
          QII.push_back((10.2*(*(mQ-2))+4.5*(*(mQ+1)))/2);
          QII.push_back((10.2*(*(mQ-2))+4.5*(*(mQ+1)))/2);
          cpt_cout[9]++;
        } 
        else if (mleft && abs(*(mQ-1)/(*(mQ+2)) -1)<=0.1) //Center mleft 2 sat
        {
          QII.push_back((8.2*(*(mQ-1))+8.2*(*(mQ+2)))/2);
          QII.push_back((8.2*(*(mQ-1))+8.2*(*(mQ+2)))/2);
          cpt_cout[10]++;
        } 
        else if (mright && abs(*(mQ-2)/(*(mQ+1)) -1)<=0.1) //Center mright 2 sat
        {
          QII.push_back((8.2*(*(mQ-2))+8.2*(*(mQ+1)))/2);
          QII.push_back((8.2*(*(mQ-2))+8.2*(*(mQ+1)))/2);
          cpt_cout[10]++;
        } 
        return QII;
        // Left mright: max = 0.367 ; 1st left neighbor = 0.292 (mean 0.330)
        // Left mright: 2nd left neighbor = 0.141 ; 1st right neighbor = 0.054 -> inverse = 7.1 and 18.5
        // Left mleft: max = 0.445 ; 1st right neighbor = 0.329 (mean 0.387)
        // Left mleft: 1st left neighbor = 0.095 ; 2nd right neighbor = 0.046 -> inverse = 10.5 and 21.7

        // Right mright: max = 0.451 ; 1st left neighbor = 0.329 (mean 0.390)
        // Right mright: 2nd left neighbor = 0.045 ; 1st right neighbor = 0.094 -> inverse = 22.2 and 10.6
        // Right mleft: max = 0.365 ; 1st right neighbor = 0.300 (mean 0.333)
        // Right mleft: 1st left neighbor = 0.054 ; 2nd right neighbor = 0.136 -> inverse = 18.5 and 7.35
        
        // Because of lost information above 255, assume the same reconstruction for Left mleft and mright, same for Right mleft and mright:
        // Left mleft/mright inverse : 8.8/2=4.4 and 20.1/2=10, mean max = 0.359
        // Right mleft/mright inverse : 20.4/2=10.2 and 8.98/2=4.5, mean max = 0.361
        // Finally : adapt. 
        

        // Center mright: max = 0.423 ; 1st left neighbor = 0.384 (mean 0.404)
        // Center mright: 2nd left and 1st right neighbor = 0.060 -> inverse = 16.67
        // Center mleft: max = 0.417 ; 1st right neighbor = 0.383 (mean 0.400)
        // Center mleft: 1st left and 2nd right neighbor = 0.062 -> inverse = 16.13
        
        // Because of lost information above 255, assume the same reconstruction for the two saturated strips:
        // inverse fixed at = 16.4/2 (two strip) = 8.2 
      }
    }

    else return Q; // no saturation --> no x-talk inversion
  }
  
  return Q;
}

// ORIGINAL
bool clusterCleaning_UNCHANGED(std::vector<int> ampls,  int crosstalkInv, uint8_t * exitCode, Float_t coeff1, Float_t coeff2, Float_t coeff3)//, Float_t coeffnn, Float_t noise)
{
/*
   if(!cluster) return true;
   std::vector<int>  ampls = convert(cluster->amplitudes());
   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true);
*/


  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // Début avec max
         if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255)
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }



        // Maximum entouré
         if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
         }
        // Fin avec un max
         if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] )
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touchée
        if(ampls.size()==1){    NofMax=1;}


  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
//
//               ____
//              |    |____
//          ____|    |    |
//         |    |    |    |____
//     ____|    |    |    |    |
//    |    |    |    |    |    |____
//  __|____|____|____|____|____|____|__
//    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
//
//   bool shapetest=true;
   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
//--------------debug:
//        if(shapecdtn==0) cout<<"NofMax: "<<NofMax<<endl;
        return shapecdtn;
      }

//      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        //Float_t coeff1=1.7;     Float_t coeff2=2.0;
        //Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255)
                                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3 ){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;

                              // debug
                                 //if (shapecdtn==0) cout<<"MaxOnStart"<<endl;

                              }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255)
                                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                 // debug
                                 //if (shapecdtn==0) cout<<"MaxOnEnd"<<endl;
                                }
                }

                if(MaxInMiddle==true){

                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*C_M+coeff2*C_D+2*coeff3 || C_M==255){
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;
                                        }
                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255)
                                           && C_Mnn<=coeff1*C_Mn+coeff2*C_M+2*coeff3
                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3) {
                                                shapecdtn=true;
                                        }
                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

//----debug:
//   if(shapecdtn==0) cout<<"Return 2, NofMax: "<<NofMax<<endl;
   return shapecdtn;
}

// WITH SATURATION CORRECTED
//bool clusterCleaning(const SiStripCluster*   cluster,  int crosstalkInv, uint8_t * exitCode)
bool clusterCleaning(std::vector<int> ampls,  int crosstalkInv, uint8_t * exitCode, Float_t coeff1, Float_t coeff2, Float_t coeff3)//, Float_t coeffnn, Float_t noise)
{
/*
   if(!cluster) return true;
   std::vector<int>  ampls = convert(cluster->amplitudes());
   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true);
*/


  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // Début avec max
         if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255)
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }



        // Maximum entouré
         if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
         }
        // Fin avec un max
         if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] )
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touchée
        if(ampls.size()==1){    NofMax=1;}


  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
//
//               ____
//              |    |____
//          ____|    |    |
//         |    |    |    |____
//     ____|    |    |    |    |
//    |    |    |    |    |    |____
//  __|____|____|____|____|____|____|__
//    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
//
//   bool shapetest=true;
   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
//--------------debug:
//        if(shapecdtn==0) cout<<"NofMax: "<<NofMax<<endl;
        return shapecdtn;
      }

//      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        //Float_t coeff1=1.7;     Float_t coeff2=2.0;
        //Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255)
                                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3 ){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;

                              // debug
                                 //if (shapecdtn==0) cout<<"MaxOnStart"<<endl;

                              }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D==255)
                                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                 // debug
                                 //if (shapecdtn==0) cout<<"MaxOnEnd"<<endl;
                                }
                }

                if(MaxInMiddle==true){

                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*C_M+coeff2*C_D+2*coeff3 || C_M>=254){
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;
                                        }
                                        if((C_Dn<=coeff1*C_D+coeff2*C_M+2*coeff3 || C_D>=254)
                                           && C_Mnn<=coeff1*C_Mn+coeff2*C_M+2*coeff3
                                           && C_Dnn<=coeff1*C_Dn+coeff2*C_D+2*coeff3) {
                                                shapecdtn=true;
                                        }
                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

//----debug:
//   if(shapecdtn==0) cout<<"Return 2, NofMax: "<<NofMax<<endl;
   return shapecdtn;
}

// NEW VERSION
bool myNewCleaning(std::vector<int> ampls, int crosstalkInv = 0, uint8_t* exitCode = nullptr)
{

  // --------------------- COUNT THE NUMBER OF MAXIMAS  ------------------------
  //----------------------------------------------------------------------------
  Int_t NofMax = 0;
  Int_t recur255 = 1;
  Int_t recur254 = 1;
  bool MaxOnStart = false, MaxInMiddle = false, MaxOnEnd = false;
  Int_t MaxPos = 0;

  // If only one strip is hit
  if (ampls.size() == 1) NofMax = 1;

  if ( crosstalkInv==0 ) {

    // Max at the start of the cluster
    if (ampls.size() > 1 &&
        ((ampls[0] > ampls[1]) ||
        (ampls.size() > 2 && ampls[0] == ampls[1] && ampls[1] > ampls[2] && ampls[0] != 254 && ampls[0] != 255) ||
        (ampls.size() == 2 && ampls[0] == ampls[1] && ampls[0] != 254 && ampls[0] != 255))) {
      NofMax = NofMax + 1;
      MaxOnStart = true;
    }

    // Max in the middle (strip between other ones)
    if (ampls.size() > 2) {
      for (unsigned int i = 1; i < ampls.size() - 1; i++) {
        if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
            (ampls.size() > 3 && i < ampls.size() - 2 && ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
            ampls[i] > ampls[i + 2] && ampls[i] != 254 && ampls[i] != 255)) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
          MaxPos = i;
        }
        if (ampls[i] == 255 && ampls[i] == ampls[i - 1]) {
          recur255 = recur255 + 1;
          MaxPos = i - (recur255 / 2);
          if (ampls[i] > ampls[i + 1]) {
            NofMax = NofMax + 1;
            MaxInMiddle = true;
          }
        }
        if (ampls[i] == 254 && ampls[i] == ampls[i - 1]) {
          recur254 = recur254 + 1;
          MaxPos = i - (recur254 / 2);
          if (ampls[i] > ampls[i + 1]) {
            NofMax = NofMax + 1;
            MaxInMiddle = true;
          }
        }
      }
    }

    // Max at the end of the cluster
    if ( (ampls.size() > 1 && (ampls[ampls.size() - 1] > ampls[ampls.size() - 2])) ||
        (ampls.size() > 2 && ampls[ampls.size() - 1] == ampls[ampls.size() - 2] &&
        ampls[ampls.size() - 2] > ampls[ampls.size() - 3] && ampls[ampls.size() - 1] != 254 && ampls[ampls.size() - 1] != 255)) {
      NofMax = NofMax + 1;
      MaxOnEnd = true;
    }

  }
  else {

    //check if double peak or not
    if (ampls.size()==2) NofMax = NofMax + 1;

    if (ampls.size() > 2) {
      for (unsigned int i = 1; i < ampls.size() - 1; i++) {
        if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
            (ampls.size() > 3 && i < ampls.size() - 2 && ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
            ampls[i] > ampls[i + 2])) {
          NofMax = NofMax + 1;
        }
      }
    }

  }



  //------------------------  SHAPE SELECTION  -----------------------------
  //------------------------------------------------------------------------
  //                              ____
  //                             |    |____
  //                         ____|    |    |
  //                        |    |    |    |____
  //                    ____|    |    |    |    |
  //                   |    |    |    |    |    |____
  //                 __|____|____|____|____|____|____|__
  //                    Vmm   Vm   M1   M2   Vp   Vpp
  //

  bool shapecdtn = false;
  if (exitCode) *exitCode = 255;

  if (crosstalkInv == 1) {
    if (NofMax == 1) {
      shapecdtn = true;
      if (exitCode) *exitCode = 0;
    }
    return shapecdtn;
  }

  Float_t M1 = 0.0;
  Float_t M2 = 0.0;
  Float_t Vm = 10000;
  Float_t Vp = 10000;
  Float_t Vmm = 10000;
  Float_t Vpp = 10000;
  Int_t CDPos;

  Float_t coeff1 = 0.17, coeff2 = 0.02, noise = 8.0;


  if (NofMax == 1) {
    if (MaxOnStart == true) {
      M1 = (Float_t)ampls[0];
      M2 = (Float_t)ampls[1];

      if (ampls.size() < 3) shapecdtn = true;
      else if (ampls.size() == 3) {
        Vp = (Float_t)ampls[2];
        if (Vp <= coeff1 * M2 + coeff2 * M1 + noise) shapecdtn = true;
        else if (exitCode) *exitCode = 2;
      }
      else if (ampls.size() > 3) {
        Vp = (Float_t)ampls[2];
        Vpp = (Float_t)ampls[3];
        if ((Vp <= coeff1 * M2 + coeff2 * M1 + noise) && Vpp <= coeff1 * Vp + coeff2 * M2 + noise) shapecdtn = true;
        else if (exitCode) *exitCode = 3;
      }
    }

    if (MaxOnEnd == true) {
      M1 = (Float_t)ampls[ampls.size() - 1];
      M2 = (Float_t)ampls[ampls.size() - 2];

      if (ampls.size() < 3) shapecdtn = true;
      else if (ampls.size() == 3) {
        Vp = (Float_t)ampls[0];
        if (Vp <= coeff1 * M2 + coeff2 * M1 + noise) shapecdtn = true;
        else if (exitCode) *exitCode = 4;
      }
      else if (ampls.size() > 3) {
        Vp = (Float_t)ampls[ampls.size() - 3];
        Vpp = (Float_t)ampls[ampls.size() - 4];
        if ((Vp <= coeff1 * M2 + coeff2 * M1 + noise) && Vpp <= coeff1 * Vp + coeff2 * M2 + noise) shapecdtn = true;
        else if (exitCode) *exitCode = 5;
      }
    }

    if (MaxInMiddle == true) {
      M1 = (Float_t)ampls[MaxPos];

      int LeftOfMaxPos = MaxPos - 1;
      if (LeftOfMaxPos <= 0) LeftOfMaxPos = 0;

      int RightOfMaxPos = MaxPos + 1;
      if (RightOfMaxPos >= (int)ampls.size()) RightOfMaxPos = ampls.size() - 1;
      //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;

      if (ampls[LeftOfMaxPos] < ampls[RightOfMaxPos]) {
        M2 = (Float_t)ampls[RightOfMaxPos];
        Vm = (Float_t)ampls[LeftOfMaxPos];
        CDPos = RightOfMaxPos;
      }
      else {
        M2 = (Float_t)ampls[LeftOfMaxPos];
        Vm = (Float_t)ampls[RightOfMaxPos];
        CDPos = LeftOfMaxPos;
      }

      if (Vm < coeff1 * M1 + coeff2 * M2 + noise || M1 >= 254) {
        if (ampls.size() == 3) shapecdtn = true;
        else if (ampls.size() > 3) {
          if (CDPos > MaxPos) {
            if (ampls.size() - CDPos - 1 == 0) {
              Vp = 0;
              Vpp = 0;
            }
            if (ampls.size() - CDPos - 1 == 1) {
              Vp = (Float_t)ampls[CDPos + 1];
              Vpp = 0;
            }
            if (ampls.size() - CDPos - 1 > 1) {
              Vp = (Float_t)ampls[CDPos + 1];
              Vpp = (Float_t)ampls[CDPos + 2];
            }
            if (MaxPos >= 2) {
              Vmm = (Float_t)ampls[MaxPos - 2];
            }
            else if (MaxPos < 2) Vmm = 0;
          }
          if (CDPos < MaxPos) {
            if (CDPos == 0) {
              Vp = 0;
              Vpp = 0;
            }
            if (CDPos == 1) {
              Vp = (Float_t)ampls[0];
              Vpp = 0;
            }
            if (CDPos > 1) {
              Vp = (Float_t)ampls[CDPos - 1];
              Vpp = (Float_t)ampls[CDPos - 2];
            }
            if (ampls.size() - LeftOfMaxPos > 1 && MaxPos + 2 < (int)(ampls.size()) - 1) {
              Vmm = (Float_t)ampls[MaxPos + 2];
            }
            else Vmm = 0;
          }
          if ((Vp <= coeff1 * M2 + coeff2 * M1 + noise || M2 >= 254) &&
              Vmm <= coeff1 * Vm + coeff2 * M1 + noise &&
              Vpp <= coeff1 * Vp + coeff2 * M2 + noise) {
            shapecdtn = true;
          }
        }
      }
      else if (exitCode) *exitCode = 6;
    }
  }
  else if (NofMax > 1 && exitCode) *exitCode = 1;  // more than one maximum

  if (ampls.size() == 1) shapecdtn = true;
  if (shapecdtn && exitCode) *exitCode = 0;

  return shapecdtn;
}


void sat_observation_MuonPU::Loop()
{
  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  //ofstream CheckDeadStrip("ROOT_SVG/CheckDeadStrip.txt", std::ofstream::out);
  ifstream CheckDeadStrip("ROOT_SVG/CheckDeadStrip.txt");
  ofstream NotDeadStrip("ROOT_SVG/NotDeadStrip.txt", std::ofstream::out);
  
  //ofstream Template_FLFR("ROOT_SVG/Template_FLFR.txt", std::ofstream::out);
  

  TFile* SaveData = new TFile("ROOT_SVG/Sat_MuonPU_alltracks.root", "RECREATE");
  //TFile* SaveData = new TFile("ROOT_SVG/PourRaph.root", "RECREATE"); 

  TH1F *Compare_eloss_sigdigi=new TH1F("Compare_eloss_over_sigdigi","Compare_eloss_over_sigdigi",1000,0,10);
  TH1F *Compare_eloss_sigdigi_nosat=new TH1F("Compare_eloss_over_sigdigi_nosat","Compare_eloss_over_sigdigi_nosat",1000,0,10);
  TH2F *Ih_VS_p_1=new TH2F("Ih_VS_p_1","Ih_VS_p_1",500,0,5,800,2,10);
  TH2F *Ih_VS_p_2=new TH2F("Ih_VS_p_2","Ih_VS_p_2",1000,0,100,800,2,10);
  TH2F *sat_cluster_1_VS_p=new TH2F("sat_cluster_1_VS_p","sat_cluster_1_VS_p",200,0,2000,10,0,10);
  TH2F *sat_cluster_1_VS_pt=new TH2F("sat_cluster_1_VS_pt","sat_cluster_1_VS_pt",300,0,300,10,0,10);
  TH2F *sat_cluster_2_VS_p=new TH2F("sat_cluster_2_VS_p","sat_cluster_2_VS_p",200,0,2000,10,0,10);
  TH2F *sat_cluster_2_VS_pt=new TH2F("sat_cluster_2_VS_pt","sat_cluster_2_VS_pt",300,0,300,10,0,10);

  TH1F *Sigdigi_left=new TH1F("Sigdigi_left","Sigdigi_left",50,0,50);
  TH1F *Sigdigi_right=new TH1F("Sigdigi_right","Sigdigi_right",50,0,50);
  TH1F *Sigdigi_center=new TH1F("Sigdigi_center","Sigdigi_center",50,0,50);
  TH1F *Sigdigi_FullLeft=new TH1F("Sigdigi_FullLeft","Sigdigi_FullLeft",50,0,50);
  TH1F *Sigdigi_FullRight=new TH1F("Sigdigi_FullRight","Sigdigi_FullRight",50,0,50);
  TH1F *Sigdigi_FullLeft_noBorder=new TH1F("Sigdigi_FullLeft_noBorder","Sigdigi_FullLeft_noBorder",50,0,50);
  TH1F *Sigdigi_FullRight_noBorder=new TH1F("Sigdigi_FullRight_noBorder","Sigdigi_FullRight_noBorder",50,0,50);
  TH1F *Sigdigi_FullLeft_Border=new TH1F("Sigdigi_FullLeft_Border","Sigdigi_FullLeft_Border",50,0,50);
  TH1F *Sigdigi_FullRight_Border=new TH1F("Sigdigi_FullRight_Border","Sigdigi_FullRight_Border",50,0,50);
  TH1F *Sigdigi_FullLeft_channel=new TH1F("Sigdigi_FullLeft_channel","Sigdigi_FullLeft_channel",768,0,768);
  TH1F *Sigdigi_FullRight_channel=new TH1F("Sigdigi_FullRight_channel","Sigdigi_FullRight_channel",768,0,768);
  TH1F *Sigdigi_FullLeft_layer=new TH1F("Sigdigi_FullLeft_layer","Sigdigi_FullLeft_layer",25,0,25);
  TH1F *Sigdigi_FullRight_layer=new TH1F("Sigdigi_FullRight_layer","Sigdigi_FullRight_layer",25,0,25);
  TH1F *Sigdigi_layer_maxsat=new TH1F("Sigdigi_layer_maxsat","Sigdigi_layer_maxsat",25,0,25);
  TH1F *Sigdigi_FullLeft_recons=new TH1F("Sigdigi_FullLeft_recons","Sigdigi_FullLeft_recons",50,0,50);
  TH1F *Sigdigi_FullRight_recons=new TH1F("Sigdigi_FullRight_recons","Sigdigi_FullRight_recons",50,0,50);
  
  TH1F *cpt_sat=new TH1F("cpt_sat","cpt_sat",10,0,10);

  TH1F *Compare_QII_Over_SigDigi=new TH1F("Compare_QII_Over_SigDigi","Compare_QII_Over_SigDigi",100,0.5,1.5);
  TH1F *Compare_QII_Over_SigDigi_L=new TH1F("Compare_QII_Over_SigDigi_L","Compare_QII_Over_SigDigi_L",100,0.5,1.5);
  TH1F *Compare_QII_Over_SigDigi_R=new TH1F("Compare_QII_Over_SigDigi_R","Compare_QII_Over_SigDigi_R",100,0.5,1.5);
  TH1F *Compare_QII_Over_SigDigi_C=new TH1F("Compare_QII_Over_SigDigi_C","Compare_QII_Over_SigDigi_C",100,0.5,1.5);
  TH1F *Compare_QII_Over_SigDigi_FL=new TH1F("Compare_QII_Over_SigDigi_FL","Compare_QII_Over_SigDigi_FL",100,0.5,1.5);
  TH1F *Compare_QII_Over_SigDigi_FR=new TH1F("Compare_QII_Over_SigDigi_FR","Compare_QII_Over_SigDigi_FR",100,0.5,1.5);

  TH1F *Compare_QII_Minus_SigDigi=new TH1F("Compare_QII_Minus_SigDigi","Compare_QII_Minus_SigDigi",80,-400,400);
  TH1F *Compare_QII_Minus_SigDigi_L=new TH1F("Compare_QII_Minus_SigDigi_L","Compare_QII_Minus_SigDigi_L",80,-400,400);
  TH1F *Compare_QII_Minus_SigDigi_R=new TH1F("Compare_QII_Minus_SigDigi_R","Compare_QII_Minus_SigDigi_R",80,-400,400);
  TH1F *Compare_QII_Minus_SigDigi_C=new TH1F("Compare_QII_Minus_SigDigi_C","Compare_QII_Minus_SigDigi_C",80,-400,400);
  TH1F *Compare_QII_Minus_SigDigi_FL=new TH1F("Compare_QII_Minus_SigDigi_FL","Compare_QII_Minus_SigDigi_FL",80,-400,400);
  TH1F *Compare_QII_Minus_SigDigi_FR=new TH1F("Compare_QII_Minus_SigDigi_FR","Compare_QII_Minus_SigDigi_FR",80,-400,400);
  
  TH1F *check_L_Original=new TH1F("check_L_Original","check_L_Original",100,0,100);
  TH1F *check_R_Original=new TH1F("check_R_Original","check_R_Original",100,0,100);
  TH1F *check_C_Original=new TH1F("check_C_Original","check_C_Original",100,0,100);
  TH1F *check_FL_Original=new TH1F("check_FL_Original","check_FL_Original",100,0,100);
  TH1F *check_FR_Original=new TH1F("check_FR_Original","check_FR_Original",100,0,100);


  TH1F *Compare_new_Correction=new TH1F("Compare_new_Correction","Compare_new_Correction",200,0,2);
  TH1F *Compare_old_Correction=new TH1F("Compare_old_Correction","Compare_old_Correction",200,0,2);
  TH1F *Compare_new_Correction_TIB=new TH1F("Compare_new_Correction_TIB","Compare_new_Correction_TIB",200,0,2);
  TH1F *Compare_new_Correction_TOB=new TH1F("Compare_new_Correction_TOB","Compare_new_Correction_TOB",200,0,2);
  TH1F *Compare_new_Correction_TID=new TH1F("Compare_new_Correction_TID","Compare_new_Correction_TID",200,0,2);
  TH1F *Compare_new_Correction_TID_R1=new TH1F("Compare_new_Correction_TID_R1","Compare_new_Correction_TID_R1",200,0,2);
  TH1F *Compare_new_Correction_TID_R2=new TH1F("Compare_new_Correction_TID_R2","Compare_new_Correction_TID_R2",200,0,2);
  TH1F *Compare_new_Correction_TID_R3=new TH1F("Compare_new_Correction_TID_R3","Compare_new_Correction_TID_R3",200,0,2);
  TH1F *Compare_new_Correction_TEC=new TH1F("Compare_new_Correction_TEC","Compare_new_Correction_TEC",200,0,2);
  TH1F *Compare_old_Correction_TIB=new TH1F("Compare_old_Correction_TIB","Compare_old_Correction_TIB",200,0,2);
  TH1F *Compare_old_Correction_TOB=new TH1F("Compare_old_Correction_TOB","Compare_old_Correction_TOB",200,0,2);
  TH1F *Compare_old_Correction_TID=new TH1F("Compare_old_Correction_TID","Compare_old_Correction_TID",200,0,2);
  TH1F *Compare_old_Correction_TID_R1=new TH1F("Compare_old_Correction_TID_R1","Compare_old_Correction_TID_R1",200,0,2);
  TH1F *Compare_old_Correction_TID_R2=new TH1F("Compare_old_Correction_TID_R2","Compare_old_Correction_TID_R2",200,0,2);
  TH1F *Compare_old_Correction_TID_R3=new TH1F("Compare_old_Correction_TID_R3","Compare_old_Correction_TID_R3",200,0,2);
  TH1F *Compare_old_Correction_TEC=new TH1F("Compare_old_Correction_TEC","Compare_old_Correction_TEC",200,0,2);

  TH1F *Compare_new_Correction_B=new TH1F("Compare_new_Correction_B","Compare_new_Correction_B",200,0,2);
  TH1F *Compare_old_Correction_B=new TH1F("Compare_old_Correction_B","Compare_old_Correction_B",200,0,2);
  TH1F *Compare_new_Correction_TIB_B=new TH1F("Compare_new_Correction_TIB_B","Compare_new_Correction_TIB_B",200,0,2);
  TH1F *Compare_new_Correction_TOB_B=new TH1F("Compare_new_Correction_TOB_B","Compare_new_Correction_TOB_B",200,0,2);
  TH1F *Compare_new_Correction_TID_B=new TH1F("Compare_new_Correction_TID_B","Compare_new_Correction_TID_B",200,0,2);
  TH1F *Compare_new_Correction_TID_R1_B=new TH1F("Compare_new_Correction_TID_R1_B","Compare_new_Correction_TID_R1_B",200,0,2);
  TH1F *Compare_new_Correction_TID_R2_B=new TH1F("Compare_new_Correction_TID_R2_B","Compare_new_Correction_TID_R2_B",200,0,2);
  TH1F *Compare_new_Correction_TID_R3_B=new TH1F("Compare_new_Correction_TID_R3_B","Compare_new_Correction_TID_R3_B",200,0,2);
  TH1F *Compare_new_Correction_TEC_B=new TH1F("Compare_new_Correction_TEC_B","Compare_new_Correction_TEC_B",200,0,2);
  TH1F *Compare_old_Correction_TIB_B=new TH1F("Compare_old_Correction_TIB_B","Compare_old_Correction_TIB_B",200,0,2);
  TH1F *Compare_old_Correction_TOB_B=new TH1F("Compare_old_Correction_TOB_B","Compare_old_Correction_TOB_B",200,0,2);
  TH1F *Compare_old_Correction_TID_B=new TH1F("Compare_old_Correction_TID_B","Compare_old_Correction_TID_B",200,0,2);
  TH1F *Compare_old_Correction_TID_R1_B=new TH1F("Compare_old_Correction_TID_R1_B","Compare_old_Correction_TID_R1_B",200,0,2);
  TH1F *Compare_old_Correction_TID_R2_B=new TH1F("Compare_old_Correction_TID_R2_B","Compare_old_Correction_TID_R2_B",200,0,2);
  TH1F *Compare_old_Correction_TID_R3_B=new TH1F("Compare_old_Correction_TID_R3_B","Compare_old_Correction_TID_R3_B",200,0,2);
  TH1F *Compare_old_Correction_TEC_B=new TH1F("Compare_old_Correction_TEC_B","Compare_old_Correction_TEC_B",200,0,2);



  TH2F *FirstNeighbour_VS_Barycenter=new TH2F("FirstNeighbour_VS_Barycenter","FirstNeighbour_VS_Barycenter",300,-3,3,100,0,10);
  TH2F *left_VS_Barycenter=new TH2F("left_VS_Barycenter","left_VS_Barycenter",300,-3,3,100,0,10);
  TH2F *right_VS_Barycenter=new TH2F("right_VS_Barycenter","right_VS_Barycenter",300,-3,3,100,0,10);
  TH2F *center_VS_Barycenter=new TH2F("center_VS_Barycenter","center_VS_Barycenter",300,-3,3,100,0,10);
  TH2F *MaxNeighbourRatio_VS_Barycenter=new TH2F("MaxNeighbourRatio_VS_Barycenter","MaxNeighbourRatio_VS_Barycenter",300,-3,3,100,0,10);
  TH2F *MaxOverMaxNeighbour_VS_MaxNeihbourRatio=new TH2F("MaxOverMaxNeighbour_VS_MaxNeihbourRatio","MaxOverMaxNeighbour_VS_MaxNeihbourRatio",100,0,10,200,0,20);
  TH2F *TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio=new TH2F("TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio","TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio",100,0,10,200,0,20);
  TH2F *TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio=new TH2F("TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio","TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio",100,0,10,200,0,20);
  TH2F *TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio=new TH2F("TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio","TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio",100,0,10,200,0,20);
  TH2F *TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio=new TH2F("TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio","TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio",100,0,10,200,0,20);

  TH2F *MaxOverMaxNeighbour_VS_Barycenter=new TH2F("MaxOverMaxNeighbour_VS_Barycenter","MaxOverMaxNeighbour_VS_Barycenter",300,-3,3,200,0,20);
  TH1F *StripSat_Over5_layer=new TH1F("StripSat_Over5_layer","StripSat_Over5_layer",21,0,21);
  TH2F *StripSat_Over5_layer_eta=new TH2F("StripSat_Over5_layer_eta","StripSat_Over5_layer_eta",21,0,21,30,-3,3);
  

  TH1F *Compare_new_Correction_BN=new TH1F("Compare_new_Correction_BN","Compare_new_Correction_BN",200,0,2);
  TH1F *Compare_old_Correction_BN=new TH1F("Compare_old_Correction_BN","Compare_old_Correction_BN",200,0,2);
  TH1F *Compare_new_Correction_TIB_BN=new TH1F("Compare_new_Correction_TIB_BN","Compare_new_Correction_TIB_BN",200,0,2);
  TH1F *Compare_new_Correction_TOB_BN=new TH1F("Compare_new_Correction_TOB_BN","Compare_new_Correction_TOB_BN",200,0,2);
  TH1F *Compare_new_Correction_TID_BN=new TH1F("Compare_new_Correction_TID_BN","Compare_new_Correction_TID_BN",200,0,2);
  TH1F *Compare_new_Correction_TEC_BN=new TH1F("Compare_new_Correction_TEC_BN","Compare_new_Correction_TEC_BN",200,0,2);
  TH1F *Compare_old_Correction_TIB_BN=new TH1F("Compare_old_Correction_TIB_BN","Compare_old_Correction_TIB_BN",200,0,2);
  TH1F *Compare_old_Correction_TOB_BN=new TH1F("Compare_old_Correction_TOB_BN","Compare_old_Correction_TOB_BN",200,0,2);
  TH1F *Compare_old_Correction_TID_BN=new TH1F("Compare_old_Correction_TID_BN","Compare_old_Correction_TID_BN",200,0,2);
  TH1F *Compare_new_Correction_TID_R1_BN=new TH1F("Compare_new_Correction_TID_R1_BN","Compare_new_Correction_TID_R1_BN",200,0,2);
  TH1F *Compare_new_Correction_TID_R2_BN=new TH1F("Compare_new_Correction_TID_R2_BN","Compare_new_Correction_TID_R2_BN",200,0,2);
  TH1F *Compare_new_Correction_TID_R3_BN=new TH1F("Compare_new_Correction_TID_R3_BN","Compare_new_Correction_TID_R3_BN",200,0,2);
  TH1F *Compare_old_Correction_TID_R1_BN=new TH1F("Compare_old_Correction_TID_R1_BN","Compare_old_Correction_TID_R1_BN",200,0,2);
  TH1F *Compare_old_Correction_TID_R2_BN=new TH1F("Compare_old_Correction_TID_R2_BN","Compare_old_Correction_TID_R2_BN",200,0,2);
  TH1F *Compare_old_Correction_TID_R3_BN=new TH1F("Compare_old_Correction_TID_R3_BN","Compare_old_Correction_TID_R3_BN",200,0,2);
  TH1F *Compare_old_Correction_TEC_BN=new TH1F("Compare_old_Correction_TEC_BN","Compare_old_Correction_TEC_BN",200,0,2);

  TH1F *ClusterUnderCorrection_layer_B=new TH1F("ClusterUnderCorrection_layer_B","ClusterUnderCorrection_layer_B",21,0,21);
  TH1F *ClusterUnderCorrection_layer_BandN=new TH1F("ClusterUnderCorrection_layer_BandN","ClusterUnderCorrection_layer_BandN",21,0,21);

  TH1F *ClusterUnderCorrection_layer=new TH1F("ClusterUnderCorrection_layer","ClusterUnderCorrection_layer",21,0,21);
  TH1F *ClusterUnderCorrection_layer_all=new TH1F("ClusterUnderCorrection_layer_all","ClusterUnderCorrection_layer_all",21,0,21);
  TH1F *ClusterUnderCorrection_shape=new TH1F("ClusterUnderCorrection_shape","ClusterUnderCorrection_shape",5,0,5);
  TH1F *ClusterUnderCorrection_shape_good=new TH1F("ClusterUnderCorrection_shape_good","ClusterUnderCorrection_shape_good",5,0,5);
  TH2F *ClusterUnderCorrection_sizeVSeta=new TH2F("ClusterUnderCorrection_sizeVSeta","ClusterUnderCorrection_sizeVSeta",30,0,30,30,-3,3);
  TH2F *ClusterGoodCorrection_sizeVSeta=new TH2F("ClusterGoodCorrection_sizeVSeta","ClusterGoodCorrection_sizeVSeta",30,0,30,30,-3,3);
  

  TH1F *Compare_Qcorr_OVER_Qini_FL=new TH1F("Compare_Qcorr_OVER_Qini_FL","Compare_Qcorr_OVER_Qini_FL",100,0.5,1.5);
  TH1F *Compare_Qcorr_MINUS_Qini_FL=new TH1F("Compare_Qcorr_MINUS_Qini_FL","Compare_Qcorr_MINUS_Qini_FL",80,-400,400);
  TH1F *Compare_Qcorr_OVER_Qini_FR=new TH1F("Compare_Qcorr_OVER_Qini_FR","Compare_Qcorr_OVER_Qini_FR",100,0.5,1.5);
  TH1F *Compare_Qcorr_MINUS_Qini_FR=new TH1F("Compare_Qcorr_MINUS_Qini_FR","Compare_Qcorr_MINUS_Qini_FR",80,-400,400);

  TH2F *MaxcorrOverMax_VS_shape_TIBL2=new TH2F("MaxcorrOverMax_VS_shape_TIBL2","MaxcorrOverMax_VS_shape_TIBL2",5,0,5,200,0.5,1.5);
  TH2F *MaxcorrOverMax_VS_size_TIBL2=new TH2F("MaxcorrOverMax_VS_size_TIBL2","MaxcorrOverMax_VS_size_TIBL2",30,0,30,200,0.5,1.5);
  TH2F *MaxcorrOverMax_VS_eta_TIBL2=new TH2F("MaxcorrOverMax_VS_eta_TIBL2","MaxcorrOverMax_VS_eta_TIBL2",30,-3,3,200,0.5,1.5);

  TH2F *MaxcorrOverMax_VS_shape_TIDR2=new TH2F("MaxcorrOverMax_VS_shape_TIDR2","MaxcorrOverMax_VS_shape_TIDR2",5,0,5,200,0.5,1.5);
  TH2F *MaxcorrOverMax_VS_size_TIDR2=new TH2F("MaxcorrOverMax_VS_size_TIDR2","MaxcorrOverMax_VS_size_TIDR2",30,0,30,200,0.5,1.5);
  TH2F *MaxcorrOverMax_VS_eta_TIDR2=new TH2F("MaxcorrOverMax_VS_eta_TIDR2","MaxcorrOverMax_VS_eta_TIDR2",30,-3,3,200,0.5,1.5);
  
  TH1F *size_FLFR=new TH1F("size_FLFR","size_FLFR",15,0,15);
  TH1F *FRFL_layer=new TH1F("FRFL_layer","FRFL_layer",25,0,25);

  std::vector <TH2D*> histoVector_IonizationNeighbour;
  std::vector <TH2D*> histoVector_MaxNeighbourBarycenter;
  std::vector <TH2D*> histoVector_MaxVmaxBarycenter;
  std::vector <TH1D*> histoCluster_FL;
  std::vector <TH1D*> histoCluster_FR;
  std::vector <int> Sum_Qini_FL;
  std::vector <int> Sum_Qini_FR;
  std::vector <TH1D*> FL_a;
  std::vector <TH1D*> FR_a;
  std::vector <TH1D*> Compare_Qcorr_OVER_Qini_FL_layer;
  for (int i=0; i<=29; i++)
  {
    TString histName = Form("MaxOverMaxNeighbour_VS_MaxNeihbourRatio_layer_%d", i);
    TH2D *histo = new TH2D(histName,histName,100,0,10,200,0,20);
    histoVector_IonizationNeighbour.push_back(histo);
  }
  for (int i=0; i<=29; i++)
  {
    TString histName = Form("MaxNeighbourRatio_VS_Barycenter_layer_%d", i);
    TH2D *histo = new TH2D(histName,histName,300,-3,3,100,0,10);
    histoVector_MaxNeighbourBarycenter.push_back(histo);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("MaxOverMaxNeighbour_VS_Barycenter_layer_%d", i);
    TH2D *histo = new TH2D(histName,histName,300,-3,3,200,0,20);
    histoVector_MaxVmaxBarycenter.push_back(histo);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("histoCluster_FL_layer_%d", i);
    TH1D *histo = new TH1D(histName,histName,50,0,50);
    histoCluster_FL.push_back(histo);
    Sum_Qini_FL.push_back(0);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("histoCluster_FR_layer_%d", i);
    TH1D *histo = new TH1D(histName,histName,50,0,50);
    histoCluster_FR.push_back(histo);
    Sum_Qini_FR.push_back(0);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("FL_a_%d", i);
    TH1D *histo = new TH1D(histName,histName,50,0,50);
    FL_a.push_back(histo);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("FR_a_%d", i);
    TH1D *histo = new TH1D(histName,histName,50,0,50);
    FR_a.push_back(histo);
  }
  for (int i=0; i<=20; i++)
  {
    TString histName = Form("Compare_Qcorr_OVER_Qini_FL_layer_%d", i);
    TH1D *histo = new TH1D(histName,histName,100,0.5,1.5);
    Compare_Qcorr_OVER_Qini_FL_layer.push_back(histo);
  }
  

  int Sum_left=0, Sum_right=0, Sum_center=0, Sum_FullLeft=0, Sum_FullRight=0, Sum_FullLeft_noBorder=0, Sum_FullRight_noBorder=0, Sum_FullLeft_Border=0, Sum_FullRight_Border=0;
  int cpt_left=0, cpt_right=0, cpt_center=0, cpt_FullLeft=0, cpt_FullRight=0, cpt_FullLeft_noBorder=0, cpt_FullRight_noBorder=0, cpt_FullLeft_Border=0, cpt_FullRight_Border=0;
  int cpt_1sat=0, cpt_sat_tot=0, cpt_only1sat=0, cpt_cluster_strip=0;
  int cpt_channel128_FullLeft=0, cpt_channel128_FullRight=0;
  int cpt_left_idem=0, cpt_right_idem=0, cpt_center_idem=0, cpt_FullLeft_idem=0, cpt_FullRight_idem=0;

  int cpt_FR_idem=0, cpt_FL_idem=0, cpt_FL_all=0, cpt_FR_all=0;
  int Sum_fit_FR=0, Sum_fit_FL=0;

  vector <int> cpt_cout_left;
  vector <int> cpt_cout_right;
  vector <int> cpt_cout_center;
  vector <int> cpt_cout_FullLeft;
  vector <int> cpt_cout_FullRight;
  vector <int> cpt_cout_all;
  vector<int> osef;
  for (int i=0; i<=10; i++)
  {
    cpt_cout_left.push_back(0);
    cpt_cout_right.push_back(0);
    cpt_cout_center.push_back(0);
    cpt_cout_FullLeft.push_back(0);
    cpt_cout_FullRight.push_back(0);
    osef.push_back(0);
    cpt_cout_all.push_back(0);
  }

  vector<int> NotDeadDetid;
  vector<int> FRFL_Detid;
  double DeltaR=0;
  double DeltaR_val;
  int index;

  int cpt_ring1=0, cpt_ring2=0, cpt_ring3=0, cpt_ring4=0, cpt_ring5=0, cpt_ring6=0, cpt_ring7=0;
  int cpt_1strip=0, cpt_2strip=0, cpt_3strip=0, cpt_4strip=0, cpt_5strip=0, cpt_6strip=0, cpt_7strip=0, cpt_more7strip=0;
  


  nentries=20000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    bool sat_cluster_1=false, sat_cluster_2=false;
    int cpt_sat=0;

    // Muon matching
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

    // Ih
    for (int itr=0; itr<ntracks; itr++)
    {
      vector <float> charge_corr;
      vector <float> pathlength;
      vector <int> subdetId;
      vector <int> moduleGeometry;
      vector <bool> bool_cleaning;
      vector <bool> mustBeInside;  

      if (track_nvalidhits[itr]>8 && track_qual[itr])
      {
        for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
        {
          if (dedx_isstrip[idedx])
          {
            float ch1=dedx_charge[idedx];
            bool clean1=true;
            bool no_in_L1_pixel=true;

            if (dedx_ispixel[idedx]) //"clean1" pixel
            {
              if (dedx_probQ[idedx]<0.8) clean1=true;
              else clean1=false;
            }
            else  //"clean1" and "ch1" strip
            {
              float check_charge=0;
              vector<int> Quncor;
              for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
              {
                check_charge+=strip_ampl[istrip];
                Quncor.push_back(strip_ampl[istrip]);
              }
              
              vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25,osef);
              float newcharge =0;
              for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
              ch1=newcharge;
              clean1=true; //true;
            }

            if (clean1 && dedx_insideTkMod[idedx] && dedx_isstrip[idedx])
            {
              charge_corr.push_back(ch1);
              pathlength.push_back(dedx_pathlength[idedx]);
              subdetId.push_back(dedx_subdetid[idedx]);
              moduleGeometry.push_back(dedx_modulgeom[idedx]);
              mustBeInside.push_back(dedx_insideTkMod[idedx]);
              bool_cleaning.push_back(clean1);
            }
          
          } //end isstrip
        } //end idedx

        int nval=0;
        int nsat=0;
        double Ih = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, NULL,2, 0., 0, nval, nsat);

        if (Ih >= -5.5*track_p[itr]+10 && track_p[itr]<=2)
        {
          Ih_VS_p_1->Fill(track_p[itr], Ih, track_prescale[itr]);
          sat_cluster_1=true;
        }
        if (track_p[itr]>2)
        {
          Ih_VS_p_2->Fill(track_p[itr], Ih, track_prescale[itr]);
          sat_cluster_2=true;
        }
      }

    }//end tracks



    // Check if a strip is dead or not:
    string line;
    while (std::getline(CheckDeadStrip, line))
    {
      std::istringstream iss(line);
      int w_detid=-1, idx_dead_strip=-1;
      if (!(iss >> w_detid >> idx_dead_strip)) break; // error
      
      for (int itr=0; itr<ntracks; itr++)
      { 
        if (track_nvalidhits[itr]>8 && track_qual[itr])
        {
          for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
          {
            if (dedx_isstrip[idedx])
            {
              if (w_detid==dedx_detid[idedx] && idx_dead_strip>=0 && idx_dead_strip<768 && sigdigi_strip_adc[idx_dead_strip]>0) NotDeadDetid.push_back(dedx_detid[idedx]);
            }
          }
        }
      }
    }


    // Full Cluster analysis
    if (track_nvalidhits[index]>8 && track_qual[index])
    {
      for (int idedx=track_index_hit[index]; idedx<track_index_hit[index]+track_nhits[index]; idedx++)
      {
        int cpt_sat_muon=0;

        if (dedx_isstrip[idedx])
        {
          cpt_cluster_strip++;
          double sum_sigdigi=0;
          double Qloss=(double)sclus_eloss[idedx]/(3.61*pow(10,-9)*265);
          vector<int> index_sig;
          vector<int> value_sig;
          vector<int> value_sig_nosat;
          double barycenter_sig=0;
          double value_max=-1;
          int i_max=-1;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              index_sig.push_back(isigdigi);
              value_sig.push_back(sigdigi_strip_adc[isigdigi]);

              if (sigdigi_strip_adc[isigdigi]<=253) value_sig_nosat.push_back(sigdigi_strip_adc[isigdigi]);
              else if (sigdigi_strip_adc[isigdigi]>253 && sigdigi_strip_adc[isigdigi]<=1023) value_sig_nosat.push_back(254);
              else if (sigdigi_strip_adc[isigdigi]>1023) value_sig_nosat.push_back(255);

              if (sigdigi_strip_adc[isigdigi]>255) cpt_sat_muon++;
            }
          }

          if (cpt_sat_muon==1) cpt_1strip++;
          else if (cpt_sat_muon==2) cpt_2strip++;
          else if (cpt_sat_muon==3) cpt_3strip++;
          else if (cpt_sat_muon==4) cpt_4strip++;
          else if (cpt_sat_muon==5) cpt_5strip++;
          else if (cpt_sat_muon==6) cpt_6strip++;
          else if (cpt_sat_muon==7) cpt_7strip++;
          else if (cpt_sat_muon>7) cpt_more7strip++;

          if (cpt_sat_muon>=5)
          {
            if (dedx_subdetid[idedx] == 3 || dedx_subdetid[idedx] == 5)
            {
              StripSat_Over5_layer->Fill(dedx_layer[idedx],1);
              StripSat_Over5_layer_eta->Fill(dedx_layer[idedx], track_eta[index], 1);
            }
            if (dedx_subdetid[idedx] == 4)
            { 
              StripSat_Over5_layer->Fill(FindLayer(4,dedx_detid[idedx]), 1);
              StripSat_Over5_layer_eta->Fill(FindLayer(4,dedx_detid[idedx]), track_eta[index], 1);
            }
            if (dedx_subdetid[idedx] == 6)
            {
              StripSat_Over5_layer->Fill(13 + FindLayer(6,dedx_detid[idedx]), 1);
              StripSat_Over5_layer_eta->Fill(13 + FindLayer(6,dedx_detid[idedx]), track_eta[index], 1);
            }
          }

          if (!index_sig.empty() && !value_sig.empty())
          {
            value_max = *max_element(value_sig.begin(), value_sig.end());
            i_max = index_sig.front()+getIndex(value_sig,value_max);
            sum_sigdigi = accumulate(value_sig.begin(), value_sig.end(), 0);

            Compare_eloss_sigdigi->Fill(Qloss/sum_sigdigi);

            for (int i=0; i<value_sig.size(); i++) barycenter_sig+=value_sig[i]*index_sig[i];
            barycenter_sig=barycenter_sig/sum_sigdigi;

            int check_OneMaxSat=0;
            for (int i=0; i<value_sig.size(); i++)
            {
              if (value_sig_nosat[i]>=254) check_OneMaxSat++;
            }

            if (check_OneMaxSat>=1)
            {
              cpt_sat_tot++;
              cpt_sat++;
            }

            if (check_OneMaxSat==0) Compare_eloss_sigdigi_nosat->Fill(Qloss/sum_sigdigi);

            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>255 && i_max>index_sig.front() && i_max<index_sig.back() && sigdigi_strip_adc[i_max-1]>1.1*sigdigi_strip_adc[i_max+1]) left=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>255 && i_max>index_sig.front() && i_max<index_sig.back() && sigdigi_strip_adc[i_max+1]>1.1*sigdigi_strip_adc[i_max-1]) right=true;
            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>255 && i_max>index_sig.front() && i_max<index_sig.back() && sigdigi_strip_adc[i_max-1]<=1.1*sigdigi_strip_adc[i_max+1] && sigdigi_strip_adc[i_max-1]>=0.9*sigdigi_strip_adc[i_max+1]) center=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && value_max>255 && i_max==index_sig.front()) FullLeft=true;
            if (index_sig.size()>=2 && check_OneMaxSat==1 && value_max>255 && i_max==index_sig.back()) FullRight=true;
            

            // 1st Max Neighbour / 1st Min neighbour VS imax-barycenter
            if (index_sig.size()>=3 && check_OneMaxSat==1 && value_max>255 && i_max>index_sig.front() && i_max<index_sig.back())
            {
              if (sigdigi_strip_adc[i_max-1]>=sigdigi_strip_adc[i_max+1])
              {
                MaxNeighbourRatio_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);
                MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (dedx_subdetid[idedx]==3) TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (dedx_subdetid[idedx]==5) TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (dedx_subdetid[idedx]==4) TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                if (dedx_subdetid[idedx]==6) MaxNeighbourRatio_VS_Barycenter->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);     
              
                if (dedx_subdetid[idedx] != 4)
                {
                  histoVector_IonizationNeighbour[dedx_layer[idedx]]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  histoVector_MaxNeighbourBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);
                }

                if (dedx_subdetid[idedx] == 4)
                {
                  int ring = FindLayer(4,dedx_detid[idedx]);
                  histoVector_IonizationNeighbour[ring]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  histoVector_MaxNeighbourBarycenter[ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);
                }
                if (dedx_subdetid[idedx] == 6)
                {
                  int ring = FindLayer(6,dedx_detid[idedx]);
                  histoVector_IonizationNeighbour[22+ring]->Fill((float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  histoVector_MaxNeighbourBarycenter[22+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);
                }
              }
              if (sigdigi_strip_adc[i_max+1]>sigdigi_strip_adc[i_max-1])
              {
                MaxNeighbourRatio_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1]);
                MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (dedx_subdetid[idedx]==3) TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (dedx_subdetid[idedx]==5) TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (dedx_subdetid[idedx]==4) TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                if (dedx_subdetid[idedx]==6) TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
              
                if (dedx_subdetid[idedx] != 4)
                {
                  histoVector_IonizationNeighbour[dedx_layer[idedx]]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  histoVector_MaxNeighbourBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1]);
                }

                if (dedx_subdetid[idedx] == 4)
                {
                  int ring = FindLayer(4,dedx_detid[idedx]);
                  histoVector_IonizationNeighbour[dedx_layer[idedx]]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  histoVector_MaxNeighbourBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1]);
                }
                if (dedx_subdetid[idedx]==6)
                {
                  int ring = FindLayer(6,dedx_detid[idedx]);
                  histoVector_IonizationNeighbour[22+ring]->Fill((float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1], (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  histoVector_MaxNeighbourBarycenter[22+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max+1]/(float)sigdigi_strip_adc[i_max-1]);
                }
              }
            }


            // Max/Vmax = f(i_max - barycenter)
            if (index_sig.size()>=2 && check_OneMaxSat==1 && value_max>255)
            {
              if (index_sig.size() == 2)
              {
                if (i_max == index_sig.front())
                {
                  MaxOverMaxNeighbour_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  
                  if (dedx_subdetid[idedx] == 3 || dedx_subdetid[idedx] == 5) histoVector_MaxVmaxBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  if (dedx_subdetid[idedx] == 4)
                  {
                    int ring = FindLayer(4,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  }
                  if (dedx_subdetid[idedx] == 6)
                  {
                    int ring = FindLayer(6,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[13+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  }
                }
                else if (i_max == index_sig.back())
                {
                  MaxOverMaxNeighbour_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                
                  if (dedx_subdetid[idedx] == 3 || dedx_subdetid[idedx] == 5) histoVector_MaxVmaxBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  if (dedx_subdetid[idedx] == 4)
                  {
                    int ring = FindLayer(4,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  }
                  if (dedx_subdetid[idedx] == 6)
                  {
                    int ring = FindLayer(6,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[13+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  }
                }
                else cout<<"Max/Vmax template: problem on clusters with size == 2"<<endl;
              }
              else
              {
                if (sigdigi_strip_adc[i_max-1]>=sigdigi_strip_adc[i_max+1])
                {
                  MaxOverMaxNeighbour_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                
                  if (dedx_subdetid[idedx] == 3 || dedx_subdetid[idedx] == 5) histoVector_MaxVmaxBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  if (dedx_subdetid[idedx] == 4)
                  {
                    int ring = FindLayer(4,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  }
                  if (dedx_subdetid[idedx] == 6)
                  {
                    int ring = FindLayer(6,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[13+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max-1]);
                  }
                }
                if (sigdigi_strip_adc[i_max+1]>sigdigi_strip_adc[i_max-1])
                {
                  MaxOverMaxNeighbour_VS_Barycenter->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                
                  if (dedx_subdetid[idedx] == 3 || dedx_subdetid[idedx] == 5) histoVector_MaxVmaxBarycenter[dedx_layer[idedx]]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  if (dedx_subdetid[idedx] == 4)
                  {
                    int ring = FindLayer(4,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  }
                  if (dedx_subdetid[idedx] == 6)
                  {
                    int ring = FindLayer(6,dedx_detid[idedx]);
                    histoVector_MaxVmaxBarycenter[13+ring]->Fill(i_max-barycenter_sig, (float)sigdigi_strip_adc[i_max]/(float)sigdigi_strip_adc[i_max+1]);
                  }
                }
              }
            }

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

            // FL and FR fit and correction

            int layer = -1;
            if (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==5) layer = dedx_layer[idedx];
            else if (dedx_subdetid[idedx]==4) layer = FindLayer(4, dedx_detid[idedx]);
            else if (dedx_subdetid[idedx]==6) layer = 13 + FindLayer(6, dedx_detid[idedx]);

            if (FullLeft)// && value_sig_nosat[i_max-index_sig.front()+1]>25)
            {
              size_FLFR->Fill(index_sig.size());
              FRFL_layer->Fill(layer);
              Sum_Qini_FL[layer]+=sum_sigdigi;
              cpt_FL_all++;

              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i < i_max) 
                {
                  histoCluster_FL[layer]->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+histoCluster_FL[layer]->GetBinContent(25+(i-i_max)));
                }
                if (i == i_max)
                {
                  histoCluster_FL[layer]->SetBinContent(25,sigdigi_strip_adc[i_max]+histoCluster_FL[layer]->GetBinContent(25));
                }
                if (i > i_max)
                {
                  histoCluster_FL[layer]->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+histoCluster_FL[layer]->GetBinContent(25+(i-i_max)));
                }
              }

              bool NotCorrected=false;
              vector<int> MaxCorr_vec = Correction_FL_FR_xtalk(value_sig_nosat, layer, NotCorrected, "ROOT_SVG/Template_FLFR.txt");
              int MaxCorr = accumulate(MaxCorr_vec.begin(), MaxCorr_vec.end(), 0);
              int Max = value_sig[i_max-index_sig.front()];

              if (NotCorrected) cpt_FL_idem++;
              else
              {
                FL_a[layer]->Fill((double) value_sig[i_max-index_sig.front()] / value_sig[i_max-index_sig.front()+1]);
                Compare_Qcorr_OVER_Qini_FL->Fill((double) MaxCorr / Max);
                Compare_Qcorr_OVER_Qini_FL_layer[layer]->Fill((double) MaxCorr / Max);
                Compare_Qcorr_MINUS_Qini_FL->Fill(MaxCorr - Max);
              }
            }

            if (FullRight) // && value_sig_nosat[i_max-index_sig.front()-1]>25)
            {
              size_FLFR->Fill(index_sig.size());
              FRFL_layer->Fill(layer);
              Sum_Qini_FR[layer]+=sum_sigdigi;
              cpt_FR_all++;

              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i < i_max) 
                {
                  histoCluster_FR[layer]->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+histoCluster_FR[layer]->GetBinContent(25+(i-i_max)));
                }
                if (i == i_max)
                {
                  histoCluster_FR[layer]->SetBinContent(25,sigdigi_strip_adc[i_max]+histoCluster_FR[layer]->GetBinContent(25));
                }
                if (i > i_max)
                {
                  histoCluster_FR[layer]->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+histoCluster_FR[layer]->GetBinContent(25+(i-i_max)));
                }
              }

              bool NotCorrected=true;
              vector<int> MaxCorr_vec = Correction_FL_FR_xtalk(value_sig_nosat, layer, NotCorrected, "ROOT_SVG/Template_FLFR.txt");
              int MaxCorr = accumulate(MaxCorr_vec.begin(), MaxCorr_vec.end(), 0);
              int Max = value_sig[i_max-index_sig.front()];

              if (NotCorrected) cpt_FR_idem++;
              else
              {
                FR_a[layer]->Fill((double) value_sig[i_max-index_sig.front()] / value_sig[i_max-index_sig.front()-1]);
                Compare_Qcorr_OVER_Qini_FR->Fill((double) MaxCorr / Max);
                Compare_Qcorr_MINUS_Qini_FR->Fill(MaxCorr - Max);
              }
            }
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


            if ( (jentry==0 || jentry==5 || jentry==16 || jentry==18) && value_sig.size()>=3 && value_sig[i_max-index_sig.front()]<=250)
            {
              //DisplayCluster(value_sig, SaveData, jentry, index, 0, 0, 0, true);
            }


// --------------------------------------------------------------------------

            // Test new correction with Vmax and Vmin
            if (check_OneMaxSat==1 && value_max>255 && i_max<index_sig.back() && i_max>index_sig.front() && index_sig.size()>2)
            {
              int Vmax = value_sig_nosat[i_max-index_sig.front()-1], Vmin = value_sig_nosat[i_max-index_sig.front()+1];
              double NewValueCorrected = -1;
              double new_TIB=-1, new_TOB=-1, new_TID=-1, new_TEC=-1;
              double old_TIB=-1, old_TOB=-1, old_TID=-1, old_TEC=-1;
              int ring = -13; 
              if (dedx_subdetid[idedx]==6) ring = FindLayer(6,dedx_detid[idedx]);

              // Vmin
              if (value_sig_nosat[i_max-index_sig.front()-1]<value_sig_nosat[i_max-index_sig.front()+1])
              {
                Vmax = value_sig_nosat[i_max-index_sig.front()+1];
                Vmin = value_sig_nosat[i_max-index_sig.front()-1];
              }
              double VmaxOverVmin = (double)Vmax/Vmin;
              double MaxOverVmax = (double)value_sig[i_max-index_sig.front()]/Vmax;

              // correction
              if (dedx_subdetid[idedx]==6) NewValueCorrected = Correction_wNeighbours(Vmin, Vmax, ring+13); // TEC (ring)
              else if (dedx_subdetid[idedx]==4) NewValueCorrected = Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx])); // TID (ring)
              else NewValueCorrected = Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx]); // TIB TOB  

              vector<int> OldCorrection = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
              int OldValueCorrected = accumulate(OldCorrection.begin(), OldCorrection.end(), 0);

              for (int i=0; i<value_sig_nosat.size(); i++)
              {
                if (i != i_max-index_sig.front()) NewValueCorrected += value_sig_nosat[i];
              }

              Compare_new_Correction->Fill(NewValueCorrected/sum_sigdigi);
              Compare_old_Correction->Fill(OldValueCorrected/sum_sigdigi);

              


              // comparison if under correction

              if (dedx_subdetid[idedx]==6)
              {
                if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, ring+13)
                 || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, ring+13))
                {
                  ClusterUnderCorrection_layer->Fill(ring+13);
                  ClusterUnderCorrection_layer_all->Fill(ring+13);
                  ClusterUnderCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape->Fill(1);
                  if (right) ClusterUnderCorrection_shape->Fill(2);
                  if (center) ClusterUnderCorrection_shape->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape->Fill(5);

                  if (jentry < 10000 && value_sig[i_max-index_sig.front()] <= 0.2*Correction_wNeighbours(Vmin, Vmax, ring+13))
                  {
                    DisplayCluster(value_sig, SaveData, jentry, index, ring+13, VmaxOverVmin, MaxOverVmax, false);
                  }
                }
                else 
                {
                  ClusterUnderCorrection_layer_all->Fill(ring+13);
                  ClusterGoodCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape_good->Fill(1);
                  if (right) ClusterUnderCorrection_shape_good->Fill(2);
                  if (center) ClusterUnderCorrection_shape_good->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape_good->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape_good->Fill(5);
                }
              }
              if (dedx_subdetid[idedx]==4)
              {
                if (dedx_layer[idedx]==11)
                {
                  MaxcorrOverMax_VS_eta_TIDR2->Fill(track_eta[index], Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))/value_sig[i_max-index_sig.front()]);
                  MaxcorrOverMax_VS_size_TIDR2->Fill(value_sig.size(), Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))/value_sig[i_max-index_sig.front()]);
                  if (left) MaxcorrOverMax_VS_shape_TIDR2->Fill(1, Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))/value_sig[i_max-index_sig.front()]);
                  if (right) MaxcorrOverMax_VS_shape_TIDR2->Fill(2, Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))/value_sig[i_max-index_sig.front()]);
                  if (center) MaxcorrOverMax_VS_shape_TIDR2->Fill(3, Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))/value_sig[i_max-index_sig.front()]);
                }

                if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))
                || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx])))
                {
                  ClusterUnderCorrection_layer->Fill(FindLayer(4, dedx_detid[idedx]));
                  ClusterUnderCorrection_layer_all->Fill(FindLayer(4, dedx_detid[idedx]));
                  ClusterUnderCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape->Fill(1);
                  if (right) ClusterUnderCorrection_shape->Fill(2);
                  if (center) ClusterUnderCorrection_shape->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape->Fill(5);

                  if (jentry < 10000 && value_sig[i_max-index_sig.front()] <= 0.2*Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx])))
                  {
                    DisplayCluster(value_sig, SaveData, jentry, index, FindLayer(4, dedx_detid[idedx]), VmaxOverVmin, MaxOverVmax, false);
                  }
                }
                else 
                {
                  ClusterUnderCorrection_layer_all->Fill(FindLayer(4, dedx_detid[idedx]));
                  ClusterGoodCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape_good->Fill(1);
                  if (right) ClusterUnderCorrection_shape_good->Fill(2);
                  if (center) ClusterUnderCorrection_shape_good->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape_good->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape_good->Fill(5);
                }
              }
              if (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==5)
              {
                if (dedx_layer[idedx]==2)
                {
                  MaxcorrOverMax_VS_eta_TIBL2->Fill(track_eta[index], Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])/value_sig[i_max-index_sig.front()]);
                  MaxcorrOverMax_VS_size_TIBL2->Fill(value_sig.size(), Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])/value_sig[i_max-index_sig.front()]);
                  if (left) MaxcorrOverMax_VS_shape_TIBL2->Fill(1, Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])/value_sig[i_max-index_sig.front()]);
                  if (right) MaxcorrOverMax_VS_shape_TIBL2->Fill(2, Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])/value_sig[i_max-index_sig.front()]);
                  if (center) MaxcorrOverMax_VS_shape_TIBL2->Fill(3, Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])/value_sig[i_max-index_sig.front()]);
                }

                if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])
                || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx]))
                {
                  ClusterUnderCorrection_layer->Fill(dedx_layer[idedx]);
                  ClusterUnderCorrection_layer_all->Fill(dedx_layer[idedx]);
                  ClusterUnderCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape->Fill(1);
                  if (right) ClusterUnderCorrection_shape->Fill(2);
                  if (center) ClusterUnderCorrection_shape->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape->Fill(5);

                  if (jentry < 10000 && value_sig[i_max-index_sig.front()] <= 0.2*Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx]))
                  {
                    DisplayCluster(value_sig, SaveData, jentry, index, dedx_layer[idedx], VmaxOverVmin, MaxOverVmax, false);
                  }
                }
                else 
                {
                  ClusterUnderCorrection_layer_all->Fill(dedx_layer[idedx]);
                  ClusterGoodCorrection_sizeVSeta->Fill(value_sig.size(), track_eta[index]);
                  if (left) ClusterUnderCorrection_shape_good->Fill(1);
                  if (right) ClusterUnderCorrection_shape_good->Fill(2);
                  if (center) ClusterUnderCorrection_shape_good->Fill(3);
                  if (FullLeft) ClusterUnderCorrection_shape_good->Fill(4);
                  if (FullRight) ClusterUnderCorrection_shape_good->Fill(5);
                }
              }


              // correction layer by layer
              if (dedx_subdetid[idedx]==3)
              {
                new_TIB=Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TIB += value_sig_nosat[i];
                }

                vector<int> Vold_TIB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TIB = accumulate(Vold_TIB.begin(), Vold_TIB.end(), 0);

                Compare_new_Correction_TIB->Fill(new_TIB/sum_sigdigi);
                Compare_old_Correction_TIB->Fill(old_TIB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==5)
              {
                new_TOB=Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TOB += value_sig_nosat[i];
                }

                vector<int> Vold_TOB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TOB = accumulate(Vold_TOB.begin(), Vold_TOB.end(), 0);

                Compare_new_Correction_TOB->Fill(new_TOB/sum_sigdigi);
                Compare_old_Correction_TOB->Fill(old_TOB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==4)
              {
                int ring_TID = FindLayer(4, dedx_detid[idedx]);
                new_TID=Correction_wNeighbours(Vmin, Vmax, ring_TID);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TID += value_sig_nosat[i];
                }

                vector<int> Vold_TID = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TID = accumulate(Vold_TID.begin(), Vold_TID.end(), 0);

                Compare_new_Correction_TID->Fill(new_TID/sum_sigdigi);
                Compare_old_Correction_TID->Fill(old_TID/sum_sigdigi);

                if (ring_TID==11) Compare_new_Correction_TID_R1->Fill(new_TID/sum_sigdigi);
                if (ring_TID==11) Compare_old_Correction_TID_R1->Fill(old_TID/sum_sigdigi);
                if (ring_TID==12) Compare_new_Correction_TID_R2->Fill(new_TID/sum_sigdigi);
                if (ring_TID==12) Compare_old_Correction_TID_R2->Fill(old_TID/sum_sigdigi);
                if (ring_TID==13) Compare_new_Correction_TID_R3->Fill(new_TID/sum_sigdigi);
                if (ring_TID==13) Compare_old_Correction_TID_R3->Fill(old_TID/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==6)
              {
                new_TEC=Correction_wNeighbours(Vmin, Vmax, ring+13);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TEC += value_sig_nosat[i];
                }

                vector<int> Vold_TEC = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TEC = accumulate(Vold_TEC.begin(), Vold_TEC.end(), 0);
                
                Compare_new_Correction_TEC->Fill(new_TEC/sum_sigdigi);
                Compare_old_Correction_TEC->Fill(old_TEC/sum_sigdigi);
              }
              else cout<<"Wrong subdetID in Correction_wNeighbours";

            }
// --------------------------------------------------------------------------


// ==========================================================================

            // Test new correction with Vmin and barycenter (old fit, correlation)
            if (check_OneMaxSat==1 && value_max>255 && index_sig.size()>1)
            {
              int Vmin=-1;
              int barycenter=0, sum_charge_nosat=0;
              double NewValueCorrected = -1;
              double new_TIB=-1, new_TOB=-1, new_TID=-1, new_TEC=-1;
              double old_TIB=-1, old_TOB=-1, old_TID=-1, old_TEC=-1;
              int ring = -13;
              if (dedx_subdetid[idedx]==6) ring = FindLayer(6,dedx_detid[idedx]);

              // Vmin
              if (i_max==index_sig.front()) Vmin = value_sig_nosat[1];  // FR and FL clusters
              else if (i_max==index_sig.back()) Vmin = value_sig_nosat[i_max-index_sig.front()-1];  // FR and FL clusters
              else
              {
                if (value_sig_nosat[i_max-index_sig.front()-1]<value_sig_nosat[i_max-index_sig.front()+1]) Vmin = value_sig_nosat[i_max-index_sig.front()-1];
                else Vmin = value_sig_nosat[i_max-index_sig.front()+1];
                if (index_sig.size()==2) cout<<"Problem size clusters"<<endl;
              }
              
              // barycenter
              for (int i=0; i<value_sig_nosat.size(); i++)
              {
                barycenter+=value_sig_nosat[i]*index_sig[i];
                sum_charge_nosat+=value_sig_nosat[i];
              }
              barycenter=barycenter/sum_charge_nosat;

              // correction
              if (dedx_subdetid[idedx]==6) NewValueCorrected = Correction_wBarycenter(Vmin, i_max, barycenter, ring+13); //TEC (ring)
              else if (dedx_subdetid[idedx]==4) NewValueCorrected = Correction_wBarycenter(Vmin, i_max, barycenter, FindLayer(4, dedx_detid[idedx])); //TID (ring)
              else NewValueCorrected = Correction_wBarycenter(Vmin, i_max, barycenter, dedx_layer[idedx]); // TIB TOB  

              vector<int> OldCorrection = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
              int OldValueCorrected = accumulate(OldCorrection.begin(), OldCorrection.end(), 0);

              for (int i=0; i<value_sig_nosat.size(); i++)
              {
                if (i != i_max-index_sig.front()) NewValueCorrected += value_sig_nosat[i];
              }

              Compare_new_Correction_B->Fill(NewValueCorrected/sum_sigdigi);
              Compare_old_Correction_B->Fill(OldValueCorrected/sum_sigdigi);

              // correction layer by layer
              if (dedx_subdetid[idedx]==3)
              {
                new_TIB = Correction_wBarycenter(Vmin, i_max, barycenter, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TIB += value_sig_nosat[i];
                }

                vector<int> Vold_TIB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TIB = accumulate(Vold_TIB.begin(), Vold_TIB.end(), 0);
                
                Compare_new_Correction_TIB_B->Fill(new_TIB/sum_sigdigi);
                Compare_old_Correction_TIB_B->Fill(old_TIB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==5)
              {
                new_TOB = Correction_wBarycenter(Vmin, i_max, barycenter, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TOB += value_sig_nosat[i];
                }

                vector<int> Vold_TOB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TOB = accumulate(Vold_TOB.begin(), Vold_TOB.end(), 0);
                
                Compare_new_Correction_TOB_B->Fill(new_TOB/sum_sigdigi);
                Compare_old_Correction_TOB_B->Fill(old_TOB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==4)
              {
                new_TID = Correction_wBarycenter(Vmin, i_max, barycenter, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TID += value_sig_nosat[i];
                }

                vector<int> Vold_TID = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TID = accumulate(Vold_TID.begin(), Vold_TID.end(), 0);
                
                Compare_new_Correction_TID_B->Fill(new_TID/sum_sigdigi);
                Compare_old_Correction_TID_B->Fill(old_TID/sum_sigdigi);

                if (FindLayer(4, dedx_detid[idedx])==11) Compare_new_Correction_TID_R1_B->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==11) Compare_old_Correction_TID_R1_B->Fill(old_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==12) Compare_new_Correction_TID_R2_B->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==12) Compare_old_Correction_TID_R2_B->Fill(old_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==13) Compare_new_Correction_TID_R3_B->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==13) Compare_old_Correction_TID_R3_B->Fill(old_TID/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==6)
              {
                new_TEC = Correction_wBarycenter(Vmin, i_max, barycenter, ring+13);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TEC += value_sig_nosat[i];
                }

                vector<int> Vold_TEC = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TEC = accumulate(Vold_TEC.begin(), Vold_TEC.end(), 0);
                
                Compare_new_Correction_TEC_B->Fill(new_TEC/sum_sigdigi);
                Compare_old_Correction_TEC_B->Fill(old_TEC/sum_sigdigi);
              }
              else cout<<"Wrong subdetID in Correction_wBarycenter";

            }
// ==========================================================================


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            // Test new correction with Vmax and barycenter (new fit, no correlation)
            if (check_OneMaxSat==1 && value_max>255 && index_sig.size()>1)
            {
              int Vmax=-1;
              int barycenter=0, sum_charge_nosat=0;
              double NewValueCorrected = -1;
              double new_TIB=-1, new_TOB=-1, new_TID=-1, new_TEC=-1;
              double old_TIB=-1, old_TOB=-1, old_TID=-1, old_TEC=-1;
              int ring = -13;
              if (dedx_subdetid[idedx]==6) ring = FindLayer(6,dedx_detid[idedx]);

              // Vmax
              if (i_max==index_sig.front()) Vmax = value_sig_nosat[1];  // FR and FL clusters
              else if (i_max==index_sig.back()) Vmax = value_sig_nosat[i_max-index_sig.front()-1];  // FR and FL clusters
              else
              {
                if (value_sig_nosat[i_max-index_sig.front()-1]<value_sig_nosat[i_max-index_sig.front()+1]) Vmax = value_sig_nosat[i_max-index_sig.front()+1];
                else Vmax = value_sig_nosat[i_max-index_sig.front()-1];
                if (index_sig.size()==2) cout<<"Problem size clusters"<<endl;
              }
              
              // barycenter
              for (int i=0; i<value_sig_nosat.size(); i++)
              {
                barycenter+=value_sig_nosat[i]*index_sig[i];
                sum_charge_nosat+=value_sig_nosat[i];
              }
              barycenter=barycenter/sum_charge_nosat;

              // correction
              if (dedx_subdetid[idedx]==6) NewValueCorrected = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, ring+13); // TEC (ring)
              else if (dedx_subdetid[idedx]==4) NewValueCorrected = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, FindLayer(4, dedx_detid[idedx])); // TID (ring)
              else NewValueCorrected = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, dedx_layer[idedx]); // TIB TOB  

              vector<int> OldCorrection = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
              int OldValueCorrected = accumulate(OldCorrection.begin(), OldCorrection.end(), 0);

              for (int i=0; i<value_sig_nosat.size(); i++)
              {
                if (i != i_max-index_sig.front()) NewValueCorrected += value_sig_nosat[i];
              }

              Compare_new_Correction_BN->Fill(NewValueCorrected/sum_sigdigi);
              Compare_old_Correction_BN->Fill(OldValueCorrected/sum_sigdigi);

              //if (NewValueCorrected < 0) cout<<"Problem : Correction negative in Correction_wNeighboursBarycenter"<<endl;

              // comparison if under correction
              if (i_max<index_sig.back() && i_max>index_sig.front() && index_sig.size()>2)
              {
                int Vmin = value_sig_nosat[i_max-index_sig.front()+1];
                if (value_sig_nosat[i_max-index_sig.front()-1]<value_sig_nosat[i_max-index_sig.front()+1])
                {
                  Vmax = value_sig_nosat[i_max-index_sig.front()+1];
                  Vmin = value_sig_nosat[i_max-index_sig.front()-1];
                }

                if (dedx_subdetid[idedx]==6)
                {
                  if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, ring+13)
                  || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, ring+13) )
                  {
                    ClusterUnderCorrection_layer_B->Fill(ring+13);
                    if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, ring+13)
                      || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, ring+13)) ClusterUnderCorrection_layer_BandN->Fill(ring+13);
                  }
                }
                if (dedx_subdetid[idedx]==4)
                {
                  if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, FindLayer(4, dedx_detid[idedx]))
                  || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, FindLayer(4, dedx_detid[idedx])) )
                  {
                    ClusterUnderCorrection_layer_B->Fill(FindLayer(4, dedx_detid[idedx]));
                    if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))
                      || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, FindLayer(4, dedx_detid[idedx]))) ClusterUnderCorrection_layer_BandN->Fill(FindLayer(4, dedx_detid[idedx]));
                  }
                }
                if (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==5)
                {
                  if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, dedx_layer[idedx])
                  || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, dedx_layer[idedx]) )
                  {
                    ClusterUnderCorrection_layer_B->Fill(dedx_layer[idedx]);
                    if (value_sig[i_max-index_sig.front()] <= 0.9*Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])
                      || value_sig[i_max-index_sig.front()] >= 1.1*Correction_wNeighbours(Vmin, Vmax, dedx_layer[idedx])) ClusterUnderCorrection_layer_BandN->Fill(dedx_layer[idedx]);
                  }
                }
              }

              // correction layer by layer
              if (dedx_subdetid[idedx]==3)
              {
                new_TIB = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TIB += value_sig_nosat[i];
                }

                vector<int> Vold_TIB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TIB = accumulate(Vold_TIB.begin(), Vold_TIB.end(), 0);
                
                Compare_new_Correction_TIB_BN->Fill(new_TIB/sum_sigdigi);
                Compare_old_Correction_TIB_BN->Fill(old_TIB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==5)
              {
                new_TOB = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, dedx_layer[idedx]);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TOB += value_sig_nosat[i];
                }

                vector<int> Vold_TOB = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TOB = accumulate(Vold_TOB.begin(), Vold_TOB.end(), 0);
                
                Compare_new_Correction_TOB_BN->Fill(new_TOB/sum_sigdigi);
                Compare_old_Correction_TOB_BN->Fill(old_TOB/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==4)
              {
                new_TID = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, FindLayer(4, dedx_detid[idedx]));
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TID += value_sig_nosat[i];
                }

                vector<int> Vold_TID = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TID = accumulate(Vold_TID.begin(), Vold_TID.end(), 0);
                
                Compare_new_Correction_TID_BN->Fill(new_TID/sum_sigdigi);
                Compare_old_Correction_TID_BN->Fill(old_TID/sum_sigdigi);

                if (FindLayer(4, dedx_detid[idedx])==11) Compare_new_Correction_TID_R1_BN->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==11) Compare_old_Correction_TID_R1_BN->Fill(old_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==12) Compare_new_Correction_TID_R2_BN->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==12) Compare_old_Correction_TID_R2_BN->Fill(old_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==13) Compare_new_Correction_TID_R3_BN->Fill(new_TID/sum_sigdigi);
                if (FindLayer(4, dedx_detid[idedx])==13) Compare_old_Correction_TID_R3_BN->Fill(old_TID/sum_sigdigi);
              }
              else if (dedx_subdetid[idedx]==6)
              {
                new_TEC = Correction_wNeighboursBarycenter(Vmax, i_max, barycenter, ring+13);
                for (int i=0; i<value_sig_nosat.size(); i++)
                {
                  if (i != i_max-index_sig.front()) new_TEC += value_sig_nosat[i];
                }

                vector<int> Vold_TEC = New_SaturationCorrection(value_sig_nosat, true, 25, osef);
                old_TEC = accumulate(Vold_TEC.begin(), Vold_TEC.end(), 0);
                
                Compare_new_Correction_TEC_BN->Fill(new_TEC/sum_sigdigi);
                Compare_old_Correction_TEC_BN->Fill(old_TEC/sum_sigdigi);
              }
              else cout<<"Wrong subdetID in Correction_wBarycenter";
            }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


            // 1st Left Neighbour / 1st Right neighbour VS imax-barycenter
            if (left || right || center || FullLeft || FullRight)
            {
              cpt_1sat++;
              Sigdigi_layer_maxsat->Fill(sigdigi_strip_layer[i_max]);

              if (left || right || center) FirstNeighbour_VS_Barycenter->Fill(i_max-barycenter_sig,(float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);
            }
            if (value_max>255 && index_sig.size()==1) cpt_only1sat++;

            if (left)
            {
              left_VS_Barycenter->Fill(i_max-barycenter_sig,(float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);

              Sum_left+=sum_sigdigi;
              cpt_left++;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_left->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_left->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_left->GetBinContent(25));
                if (i>i_max) Sigdigi_left->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left->GetBinContent(25+(i-i_max)));
              }
              //vector<int> QLeft = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_left);
              vector<int> QLeft = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QLeft = accumulate(QLeft.begin(), QLeft.end(), 0);
              if (compareVectors(QLeft,value_sig_nosat)) cpt_left_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QLeft/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QLeft-sum_sigdigi);
                Compare_QII_Minus_SigDigi_L->Fill(Sum_QLeft-sum_sigdigi);
                Compare_QII_Over_SigDigi_L->Fill((float)Sum_QLeft/sum_sigdigi);
                vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                check_L_Original->Fill(100*(float)(*(mQ-1))/sum_sigdigi);

                /*if ((jentry==4428 && idedx==13) || (jentry==3420 && idedx==12))
                {
                  cout<<"Original: "<<endl;
                  for (int i=0; i<value_sig.size();i++) cout<<value_sig[i]<<" ";
                  cout<<endl;
                  vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                  cout<<*(mQ-1)<<"   ("<<100*(float)(*(mQ-1))/sum_sigdigi<<" %)"<<endl;
                  for (int i=0; i<value_sig_nosat.size();i++) cout<<value_sig_nosat[i]<<" ";
                  cout<<endl;
                  cout<<"Corrected: "<<endl;
                  for (int i=0; i<QLeft.size();i++) cout<<QLeft[i]<<" ";
                  cout<<endl;
                }*/
              }
            }

            if (right)
            {
              right_VS_Barycenter->Fill(i_max-barycenter_sig,(float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);

              Sum_right+=sum_sigdigi;
              cpt_right++;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_right->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_right->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_right->GetBinContent(25));
                if (i>i_max) Sigdigi_right->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right->GetBinContent(25+(i-i_max)));
              }
              //vector<int> QRight = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_right);
              vector<int> QRight = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QRight = accumulate(QRight.begin(), QRight.end(), 0);
              if (compareVectors(QRight,value_sig_nosat)) cpt_right_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QRight/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QRight-sum_sigdigi);
                Compare_QII_Minus_SigDigi_R->Fill(Sum_QRight-sum_sigdigi);
                Compare_QII_Over_SigDigi_R->Fill((float)Sum_QRight/sum_sigdigi);
                vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                check_R_Original->Fill(100*(float)(*(mQ+1))/sum_sigdigi);
              }
            }

            if (center)
            {
              center_VS_Barycenter->Fill(i_max-barycenter_sig,(float)sigdigi_strip_adc[i_max-1]/(float)sigdigi_strip_adc[i_max+1]);

              Sum_center+=sum_sigdigi;
              cpt_center++;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_center->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_center->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_center->GetBinContent(25));
                if (i>i_max) Sigdigi_center->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center->GetBinContent(25+(i-i_max)));
              }
              //vector<int> QCenter = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_center);
              vector<int> QCenter = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QCenter = accumulate(QCenter.begin(), QCenter.end(), 0);
              if (compareVectors(QCenter,value_sig_nosat)) cpt_center_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QCenter/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QCenter-sum_sigdigi);
                Compare_QII_Minus_SigDigi_C->Fill(Sum_QCenter-sum_sigdigi);
                Compare_QII_Over_SigDigi_C->Fill((float)Sum_QCenter/sum_sigdigi);
                vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                check_C_Original->Fill(100*(float)(*(mQ-1))/sum_sigdigi);
              }
            }

            if (FullLeft)
            {
              Sum_FullLeft+=sum_sigdigi;
              cpt_FullLeft++;
              int istrip=sclus_index_strip[idedx];
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (istrip>sclus_index_strip[idedx]+sclus_nstrip[idedx]) break;
                if (i<i_max) 
                {
                  Sigdigi_FullLeft->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft->GetBinContent(25+(i-i_max)));
                  Sigdigi_FullLeft_recons->SetBinContent(25+(i-i_max),strip_ampl[istrip]+Sigdigi_FullLeft_recons->GetBinContent(25+(i-i_max)));
                }
                if (i==i_max)
                {
                  Sigdigi_FullLeft->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullLeft->GetBinContent(25));
                  Sigdigi_FullLeft_recons->SetBinContent(25,strip_ampl[istrip]+Sigdigi_FullLeft_recons->GetBinContent(25));
                }
                if (i>i_max)
                {
                  Sigdigi_FullLeft->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft->GetBinContent(25+(i-i_max)));
                  Sigdigi_FullLeft_recons->SetBinContent(25+(i-i_max),strip_ampl[istrip]+Sigdigi_FullLeft_recons->GetBinContent(25+(i-i_max)));
                }
                istrip++;
              }

              Sigdigi_FullLeft_layer->Fill(sigdigi_strip_layer[i_max]);
              Sigdigi_FullLeft_channel->Fill(sigdigi_strip_channel[i_max]);
              if (sigdigi_strip_channel[i_max]%128!=0)
              {
                Sum_FullLeft_noBorder+=sum_sigdigi;
                cpt_FullLeft_noBorder++;
                for (int i=index_sig.front(); i<=index_sig.back(); i++)
                {
                  if (i<i_max) Sigdigi_FullLeft_noBorder->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft_noBorder->GetBinContent(25+(i-i_max)));
                  if (i==i_max) Sigdigi_FullLeft_noBorder->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullLeft_noBorder->GetBinContent(25));
                  if (i>i_max) Sigdigi_FullLeft_noBorder->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft_noBorder->GetBinContent(25+(i-i_max)));
                }
              }
              else
              {
                cpt_channel128_FullLeft++;
                
                Sum_FullLeft_Border+=sum_sigdigi;
                cpt_FullLeft_Border++;
                for (int i=index_sig.front(); i<=index_sig.back(); i++)
                {
                  if (i<i_max) Sigdigi_FullLeft_Border->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft_Border->GetBinContent(25+(i-i_max)));
                  if (i==i_max) Sigdigi_FullLeft_Border->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullLeft_Border->GetBinContent(25));
                  if (i>i_max) Sigdigi_FullLeft_Border->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft_Border->GetBinContent(25+(i-i_max)));
                }
              }

              //vector<int> QFullLeft = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_FullLeft);
              vector<int> QFullLeft = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QFullLeft = accumulate(QFullLeft.begin(), QFullLeft.end(), 0);
              if (compareVectors(QFullLeft,value_sig_nosat)) cpt_FullLeft_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QFullLeft/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QFullLeft-sum_sigdigi);
                Compare_QII_Minus_SigDigi_FL->Fill(Sum_QFullLeft-sum_sigdigi);
                Compare_QII_Over_SigDigi_FL->Fill((float)Sum_QFullLeft/sum_sigdigi);
                vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                check_FL_Original->Fill(100*(float)(*(mQ+1))/sum_sigdigi);
              }

              //CheckDeadStrip << dedx_detid[idedx] << " "<< i_max-1 << "\n";
              FRFL_Detid.push_back(dedx_detid[idedx]);
            }

            if (FullRight)
            {
              Sum_FullRight+=sum_sigdigi;
              cpt_FullRight++;
              int istrip=sclus_index_strip[idedx];
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (istrip>sclus_index_strip[idedx]+sclus_nstrip[idedx]) break;
                if (i<i_max) 
                {
                  Sigdigi_FullRight->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight->GetBinContent(25+(i-i_max)));
                  Sigdigi_FullRight_recons->SetBinContent(25+(i-i_max),strip_ampl[istrip]+Sigdigi_FullRight_recons->GetBinContent(25+(i-i_max)));
                }
                if (i==i_max)
                {
                  Sigdigi_FullRight->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullRight->GetBinContent(25));
                  Sigdigi_FullRight_recons->SetBinContent(25,strip_ampl[istrip]+Sigdigi_FullRight_recons->GetBinContent(25));
                }
                if (i>i_max)
                {
                  Sigdigi_FullRight->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight->GetBinContent(25+(i-i_max)));
                  Sigdigi_FullRight_recons->SetBinContent(25+(i-i_max),strip_ampl[istrip]+Sigdigi_FullRight_recons->GetBinContent(25+(i-i_max)));
                }
                istrip++;
              }

              Sigdigi_FullRight_layer->Fill(sigdigi_strip_layer[i_max]);
              Sigdigi_FullRight_channel->Fill(sigdigi_strip_channel[i_max]);
              if ((sigdigi_strip_channel[i_max]+1)%128!=0)
              {
                Sum_FullRight_noBorder+=sum_sigdigi;
                cpt_FullRight_noBorder++;
                for (int i=index_sig.front(); i<=index_sig.back(); i++)
                {
                  if (i<i_max) Sigdigi_FullRight_noBorder->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight_noBorder->GetBinContent(25+(i-i_max)));
                  if (i==i_max) Sigdigi_FullRight_noBorder->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullRight_noBorder->GetBinContent(25));
                  if (i>i_max) Sigdigi_FullRight_noBorder->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight_noBorder->GetBinContent(25+(i-i_max)));
                }
                //cout<<"r "<<i_max+1<<" "<<dedx_detid[idedx]<<endl;
                //check_dead_strip1->Fill(dedx_detid[idedx]);
                //check_dead_strip2->Fill(i_max+1);
                //check_dead_strip->Fill(dedx_detid[idedx],i_max+1);
                //cout<<"r end"<<endl;
              }
              else
              {
                cpt_channel128_FullRight++;

                Sum_FullRight_Border+=sum_sigdigi;
                cpt_FullRight_Border++;
                for (int i=index_sig.front(); i<=index_sig.back(); i++)
                {
                  if (i<i_max) Sigdigi_FullRight_Border->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight_Border->GetBinContent(25+(i-i_max)));
                  if (i==i_max) Sigdigi_FullRight_Border->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullRight_Border->GetBinContent(25));
                  if (i>i_max) Sigdigi_FullRight_Border->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight_Border->GetBinContent(25+(i-i_max)));
                }
              }

              //vector<int> QFullRight = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_FullRight);
              vector<int> QFullRight = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QFullRight = accumulate(QFullRight.begin(), QFullRight.end(), 0);
              if (compareVectors(QFullRight,value_sig_nosat)) cpt_FullRight_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QFullRight/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QFullRight-sum_sigdigi);
                Compare_QII_Minus_SigDigi_FR->Fill(Sum_QFullRight-sum_sigdigi);
                Compare_QII_Over_SigDigi_FR->Fill((float)Sum_QFullRight/sum_sigdigi);
                vector<int>::const_iterator mQ = max_element(value_sig.begin(), value_sig.end());
                check_FR_Original->Fill(100*(float)(*(mQ-1))/sum_sigdigi);
              }

              //CheckDeadStrip << dedx_detid[idedx] << " "<< i_max+1 << "\n";
              FRFL_Detid.push_back(dedx_detid[idedx]);
            }
            
          } //end !vector.empty()

          if (dedx_subdetid[idedx]==6)
          {
            int ring = FindLayer(6,dedx_detid[idedx]);
            if (ring==1) cpt_ring1++;
            else if (ring==2) cpt_ring2++;
            else if (ring==3) cpt_ring3++;
            else if (ring==4) cpt_ring4++;
            else if (ring==5) cpt_ring5++;
            else if (ring==6) cpt_ring6++;
            else if (ring==7) cpt_ring7++;
            else cout<<"ERROR: ring number greater than 7"<<endl;
          }

        } //end dedx_isstrip
      } //end idedx
    } //end track

    if (sat_cluster_1)
    {
      sat_cluster_1_VS_p->Fill(track_p[index], cpt_sat, track_prescale[index]);
      sat_cluster_1_VS_pt->Fill(track_pt[index], cpt_sat, track_prescale[index]);
    }
    if (sat_cluster_2)
    {
      sat_cluster_2_VS_p->Fill(track_p[index], cpt_sat, track_prescale[index]);
      sat_cluster_2_VS_pt->Fill(track_pt[index], cpt_sat, track_prescale[index]);
    }

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;

  } //end EVENTS

  std::sort(FRFL_Detid.begin(),FRFL_Detid.end());
  std::sort(NotDeadDetid.begin(),NotDeadDetid.end());
  // All FL and FR detid
  NotDeadStrip << FRFL_Detid[0] << "\n";
  for (int i=1;i<FRFL_Detid.size();i++)
  {
    if (FRFL_Detid[i]!=FRFL_Detid[i-1]) NotDeadStrip << FRFL_Detid[i] << "\n";
  }
  // FL and FR detid for which neighbouring strip is alive
  NotDeadStrip << "\n";
  NotDeadStrip << NotDeadDetid[0] << "\n";
  for (int i=1;i<NotDeadDetid.size();i++)
  {
    if (NotDeadDetid[i]!=NotDeadDetid[i-1]) NotDeadStrip << NotDeadDetid[i] << "\n";
  }

  Sigdigi_left->Scale(1./Sum_left);
  Sigdigi_right->Scale(1./Sum_right);
  Sigdigi_center->Scale(1./Sum_center);
  Sigdigi_FullLeft->Scale(1./Sum_FullLeft);
  Sigdigi_FullRight->Scale(1./Sum_FullRight);
  Sigdigi_FullLeft_recons->Scale(1./Sum_FullLeft);
  Sigdigi_FullRight_recons->Scale(1./Sum_FullRight);
  Sigdigi_FullLeft_noBorder->Scale(1./Sum_FullLeft_noBorder);
  Sigdigi_FullRight_noBorder->Scale(1./Sum_FullRight_noBorder);
  Sigdigi_FullLeft_Border->Scale(1./Sum_FullLeft_Border);
  Sigdigi_FullRight_Border->Scale(1./Sum_FullRight_Border);
  
  cpt_sat->SetBinContent(1,(float)cpt_left/(float)Sum_left);
  cpt_sat->SetBinContent(2,(float)cpt_right/(float)Sum_right);
  cpt_sat->SetBinContent(3,(float)cpt_center/(float)Sum_center);
  cpt_sat->SetBinContent(4,(float)cpt_FullLeft/(float)Sum_FullLeft);
  cpt_sat->SetBinContent(5,(float)cpt_FullRight/(float)Sum_FullRight);
  cpt_sat->SetBinContent(6,(float)cpt_FullLeft_noBorder/(float)Sum_FullLeft_noBorder);
  cpt_sat->SetBinContent(7,(float)cpt_FullRight_noBorder/(float)Sum_FullRight_noBorder);
  cpt_sat->SetBinContent(8,(float)cpt_FullLeft_Border/(float)Sum_FullLeft_Border);
  cpt_sat->SetBinContent(9,(float)cpt_FullRight_Border/(float)Sum_FullRight_Border);
  

  SaveData->cd();
  Compare_eloss_sigdigi->Write();
  Compare_eloss_sigdigi_nosat->Write();
  Sigdigi_left->Write();
  Sigdigi_right->Write();
  Sigdigi_center->Write();
  Sigdigi_FullLeft->Write();
  Sigdigi_FullRight->Write();
  Sigdigi_FullLeft_recons->Write();
  Sigdigi_FullRight_recons->Write();
  Sigdigi_FullLeft_noBorder->Write();
  Sigdigi_FullRight_noBorder->Write();
  Sigdigi_FullLeft_Border->Write();
  Sigdigi_FullRight_Border->Write();
  Sigdigi_FullLeft_channel->Write();
  Sigdigi_FullRight_channel->Write();
  Sigdigi_FullLeft_layer->Write();
  Sigdigi_FullRight_layer->Write();
  Sigdigi_layer_maxsat->Write();
  cpt_sat->Write();
  Compare_QII_Over_SigDigi->Write();
  Compare_QII_Minus_SigDigi->Write();
  Compare_QII_Minus_SigDigi_L->Write();
  Compare_QII_Minus_SigDigi_R->Write();
  Compare_QII_Minus_SigDigi_C->Write();
  Compare_QII_Minus_SigDigi_FL->Write();
  Compare_QII_Minus_SigDigi_FR->Write();
  Compare_QII_Over_SigDigi_L->Write();
  Compare_QII_Over_SigDigi_R->Write();
  Compare_QII_Over_SigDigi_C->Write();
  Compare_QII_Over_SigDigi_FL->Write();
  Compare_QII_Over_SigDigi_FR->Write();

  check_L_Original->Write();
  check_R_Original->Write();
  check_C_Original->Write();
  check_FL_Original->Write();
  check_FR_Original->Write();

  Ih_VS_p_1->Write();
  Ih_VS_p_2->Write();
  sat_cluster_1_VS_p->Write();
  sat_cluster_1_VS_pt->Write();
  sat_cluster_2_VS_p->Write();
  sat_cluster_2_VS_pt->Write();

  FirstNeighbour_VS_Barycenter->Write();
  left_VS_Barycenter->Write();
  right_VS_Barycenter->Write();
  center_VS_Barycenter->Write();
  MaxNeighbourRatio_VS_Barycenter->Write();
  MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Write();
  TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Write();
  TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Write();
  TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Write();
  TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Write();

  for (int i=0; i<histoVector_IonizationNeighbour.size(); i++) histoVector_IonizationNeighbour[i]->Write();
  for (int i=0; i<histoVector_MaxNeighbourBarycenter.size(); i++) histoVector_MaxNeighbourBarycenter[i]->Write();

  MaxOverMaxNeighbour_VS_Barycenter->Write();
  for (int i=0; i<histoVector_MaxVmaxBarycenter.size(); i++) histoVector_MaxVmaxBarycenter[i]->Write();


  Compare_new_Correction->Write();
  Compare_old_Correction->Write();
  Compare_new_Correction_TIB->Write();
  Compare_old_Correction_TIB->Write();
  Compare_new_Correction_TOB->Write();
  Compare_old_Correction_TOB->Write();
  Compare_new_Correction_TID->Write();
  Compare_old_Correction_TID->Write();
  Compare_new_Correction_TID_R1->Write();
  Compare_old_Correction_TID_R1->Write();
  Compare_new_Correction_TID_R2->Write();
  Compare_old_Correction_TID_R2->Write();
  Compare_new_Correction_TID_R3->Write();
  Compare_old_Correction_TID_R3->Write();
  Compare_new_Correction_TEC->Write();
  Compare_old_Correction_TEC->Write();


  Compare_new_Correction_B->Write();
  Compare_old_Correction_B->Write();
  Compare_new_Correction_TIB_B->Write();
  Compare_old_Correction_TIB_B->Write();
  Compare_new_Correction_TOB_B->Write();
  Compare_old_Correction_TOB_B->Write();
  Compare_new_Correction_TID_B->Write();
  Compare_old_Correction_TID_B->Write();
  Compare_new_Correction_TID_R1_B->Write();
  Compare_old_Correction_TID_R1_B->Write();
  Compare_new_Correction_TID_R2_B->Write();
  Compare_old_Correction_TID_R2_B->Write();
  Compare_new_Correction_TID_R3_B->Write();
  Compare_old_Correction_TID_R3_B->Write();
  Compare_new_Correction_TEC_B->Write();
  Compare_old_Correction_TEC_B->Write();

  StripSat_Over5_layer->Write();
  StripSat_Over5_layer_eta->Write();

  Compare_new_Correction_BN->Write();
  Compare_old_Correction_BN->Write();
  Compare_new_Correction_TIB_BN->Write();
  Compare_old_Correction_TIB_BN->Write();
  Compare_new_Correction_TOB_BN->Write();
  Compare_old_Correction_TOB_BN->Write();
  Compare_new_Correction_TID_BN->Write();
  Compare_old_Correction_TID_BN->Write();
  Compare_new_Correction_TID_R1_BN->Write();
  Compare_old_Correction_TID_R1_BN->Write();
  Compare_new_Correction_TID_R2_BN->Write();
  Compare_old_Correction_TID_R2_BN->Write();
  Compare_new_Correction_TID_R3_BN->Write();
  Compare_old_Correction_TID_R3_BN->Write();
  Compare_new_Correction_TEC_BN->Write();
  Compare_old_Correction_TEC_BN->Write();

  ClusterUnderCorrection_layer_B->Write();
  ClusterUnderCorrection_layer_BandN->Write();

  ClusterUnderCorrection_layer->Write();
  ClusterUnderCorrection_layer_all->Write();
  ClusterUnderCorrection_shape->Write();
  ClusterUnderCorrection_shape_good->Write();
  ClusterUnderCorrection_sizeVSeta->Write();
  ClusterGoodCorrection_sizeVSeta->Write();

  Compare_Qcorr_OVER_Qini_FR->Write();
  Compare_Qcorr_MINUS_Qini_FR->Write();
  Compare_Qcorr_OVER_Qini_FL->Write();
  Compare_Qcorr_MINUS_Qini_FL->Write();
  for (int i=0; i<histoCluster_FL.size(); i++)
  {
    histoCluster_FL[i]->Scale(1./Sum_Qini_FL[i]);
    histoCluster_FL[i]->Write();
  }
  for (int i=0; i<histoCluster_FR.size(); i++)
  {
    histoCluster_FR[i]->Scale(1./Sum_Qini_FR[i]);
    histoCluster_FR[i]->Write();
  }

  MaxcorrOverMax_VS_eta_TIBL2->Write();
  MaxcorrOverMax_VS_size_TIBL2->Write();
  MaxcorrOverMax_VS_shape_TIBL2->Write();
  MaxcorrOverMax_VS_eta_TIDR2->Write();
  MaxcorrOverMax_VS_size_TIDR2->Write();
  MaxcorrOverMax_VS_shape_TIDR2->Write();
  size_FLFR->Write();
  for (int i=0; i<Compare_Qcorr_OVER_Qini_FL_layer.size(); i++) 
  {
    Compare_Qcorr_OVER_Qini_FL_layer[i]->Write();
    FR_a[i]->Write();
    FL_a[i]->Write();
  }
  FRFL_layer->Scale(1./FRFL_layer->Integral());
  FRFL_layer->Write();

  // Write TEMPLATE
  //for (int i=1;i<FL_a.size();i++) Template_FLFR << i << " " << FL_a[i]->GetMean() << " " << FR_a[i]->GetMean() << "\n";
  

  

  /*cout << "------------ FL and FR clusters NOT CORRECTED ------------" << endl;
  cout << "FL: " << cpt_FL_idem << " / " << cpt_FL_all << endl;
  cout << "FR: " << cpt_FR_idem << " / " << cpt_FR_all << endl;

  cout << "----------------------------------------" << endl;
  cout << "FL Max/Vmax coeff: " << endl;
  for (int i=1; i<FL_a.size(); i++) cout << FL_a[i]->GetMean() << endl;
  cout << "----------------------------------------" << endl;
  cout << "FR Max/Vmax coeff: " << endl;
  for (int i=1; i<FR_a.size(); i++) cout << FR_a[i]->GetMean() << endl;
  // */

/*
  cout<<"SATURATED STRIP COUNTER - MuonPU"<<endl;
  cout<<"cpt_1sat: "<<cpt_1sat<<endl;
  cout<<"1 strip: "<<cpt_1strip<<endl;
  cout<<"2 strips: "<<cpt_2strip<<endl;
  cout<<"3 strips: "<<cpt_3strip<<endl;
  cout<<"4 strips: "<<cpt_4strip<<endl;
  cout<<"5 strips: "<<cpt_5strip<<endl;
  cout<<"6 strips: "<<cpt_6strip<<endl;
  cout<<"7 strips: "<<cpt_7strip<<endl;
  cout<<"More than 7 strips: "<<cpt_more7strip<<endl;
  cout<<"TOTAL: "<<cpt_sat_tot<<endl;
  cout<<"TOTAL CLUSTERS IN STRIP: "<<cpt_cluster_strip<<endl;

  cout<<"TOT sat: "<<cpt_sat_tot<<endl;
  cout<<"1 strip sat: "<<cpt_1sat<<endl;
  cout<<"Only 1 strip sat: "<<cpt_only1sat<<endl;
  cout<<"left: "<<cpt_left<<endl;
  cout<<"right: "<<cpt_right<<endl;
  cout<<"center: "<<cpt_center<<endl;
  cout<<"FullLeft: "<<cpt_FullLeft<<endl;
  cout<<"FullRight: "<<cpt_FullRight<<endl;
  cout<<"FullLeft no Border: "<<cpt_FullLeft_noBorder<<endl;
  cout<<"FullRight no Border: "<<cpt_FullRight_noBorder<<endl;
  

  cout<<"Channel % 128 == 0 FullLeft: "<<cpt_channel128_FullLeft<<endl;
  cout<<"Channel+1 % 128 == 0 FullRight: "<<cpt_channel128_FullRight<<endl;


  cout<<"         LEFT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_left_idem<<" / "<<cpt_left<<"  ("<<100*(double)cpt_left_idem/cpt_left<<" %)"<<endl;
  for (int i=0; i<cpt_cout_left.size(); i++) cout<<"cpt_cout_left["<<i<<"]  "<<cpt_cout_left[i]<<"   "<<100*(double)cpt_cout_left[i]/cpt_left_idem<<" %"<<endl;
  cout<<"         RIGHT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_right_idem<<" / "<<cpt_right<<"  ("<<100*(double)cpt_right_idem/cpt_right<<" %)"<<endl;
  for (int i=0; i<cpt_cout_right.size(); i++) cout<<"cpt_cout_right["<<i<<"]  "<<cpt_cout_right[i]<<"   "<<100*(double)cpt_cout_right[i]/cpt_right_idem<<" %"<<endl;
  cout<<"         CENTER"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_center_idem<<" / "<<cpt_center<<"  ("<<100*(double)cpt_center_idem/cpt_center<<" %)"<<endl;
  for (int i=0; i<cpt_cout_center.size(); i++) cout<<"cpt_cout_center["<<i<<"]  "<<cpt_cout_center[i]<<"   "<<100*(double)cpt_cout_center[i]/cpt_center_idem<<" %"<<endl;
  cout<<"         FULL LEFT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_FullLeft_idem<<" / "<<cpt_FullLeft<<"  ("<<100*(double)cpt_FullLeft_idem/cpt_FullLeft<<" %)"<<endl;
  for (int i=0; i<cpt_cout_FullLeft.size(); i++) cout<<"cpt_cout_FullLeft["<<i<<"]  "<<cpt_cout_FullLeft[i]<<"   "<<100*(double)cpt_cout_FullLeft[i]/cpt_FullLeft_idem<<" %"<<endl;
  cout<<"         FULL RIGHT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_FullRight_idem<<" / "<<cpt_FullRight<<"  ("<<100*(double)cpt_FullRight_idem/cpt_FullRight<<" %)"<<endl;
  for (int i=0; i<cpt_cout_FullRight.size(); i++) cout<<"cpt_cout_FullRight["<<i<<"]  "<<cpt_cout_FullRight[i]<<"   "<<100*(double)cpt_cout_FullRight[i]/cpt_FullRight_idem<<" %"<<endl;  


  cout<<"         NBR OF CLUSTERS RECONSTRUCTED"<<endl;
  cout<<"nbr of clusters rejected  "<<(cpt_left_idem+cpt_right_idem+cpt_center_idem+cpt_FullLeft_idem+cpt_FullRight_idem)<<" / "<<(cpt_left+cpt_right+cpt_center+cpt_FullLeft+cpt_FullRight)<<"  ("<<100*(double)(cpt_left_idem+cpt_right_idem+cpt_center_idem+cpt_FullLeft_idem+cpt_FullRight_idem)/(cpt_left+cpt_right+cpt_center+cpt_FullLeft+cpt_FullRight)<<" %)"<<endl;
  for (int i=0; i<cpt_cout_all.size(); i++) cout<<"cpt_cout_all["<<i<<"]  "<<cpt_cout_all[i]<<endl;

  cout<<endl;
  cout<<"            RING COUNTER TEST"<<endl;
  cout<<"ring 1: "<<cpt_ring1<<endl;
  cout<<"ring 2: "<<cpt_ring2<<endl;
  cout<<"ring 3: "<<cpt_ring3<<endl;
  cout<<"ring 4: "<<cpt_ring4<<endl;
  cout<<"ring 5: "<<cpt_ring5<<endl;
  cout<<"ring 6: "<<cpt_ring6<<endl;
  cout<<"ring 7: "<<cpt_ring7<<endl;
*/
  
  CheckDeadStrip.close();
}
