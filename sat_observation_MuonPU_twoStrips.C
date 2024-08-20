#define sat_observation_MuonPU_twoStrips_cxx
#include "sat_observation_MuonPU_twoStrips.h"
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
        // 1st left neighbor = 0.108 -> inverse = 1/0.108 = 9.26
      }

      if(Q.size()>2 && *(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 && (*(mQ-1)>1.1*(*(mQ+1)) || *(mQ+1)>1.1*(*(mQ-1)) || abs(*(mQ-1)/(*(mQ+1)) -1)<=0.1))
      {
        if (*(mQ-1)>1.1*(*(mQ+1)))
        {
          QII.push_back((4.78*(*(mQ-1))+13.51*(*(mQ+1)))/2); //Left one sat
          cpt_cout[3]++;
        }
        else if (*(mQ+1)>1.1*(*(mQ-1)))
        {
          QII.push_back((13.51*(*(mQ-1))+4.78*(*(mQ+1)))/2); //Right one sat
          cpt_cout[4]++;
        }
        else if (abs(*(mQ-1)/(*(mQ+1)) -1)<=0.1)
        {
          QII.push_back((11.90*(*(mQ-1))+11.90*(*(mQ+1)))/2);  //Center one sat
          cpt_cout[5]++;
        }
        return QII;
        // Left/Right: max = 0.617 (sat)
        // Left/Right: 1st right neighbor = 0.209 (no sat) -> inverse = 4.78
        // Left/Right: 1st left neighbor = 0.074  (no sat) -> inverse = 13.51
        
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


void sat_observation_MuonPU_twoStrips::Loop()
{
  if (fChain == 0) return -1;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  TFile* SaveData = new TFile("ROOT_SVG/Sat_MuonPU_twoStrips_alltracks.root", "RECREATE");

  TH1F *Compare_eloss_sigdigi=new TH1F("Compare_eloss-sigdigi","Compare_eloss-sigdigi",1000,0,10);
  TH2F *Ih_VS_p=new TH2F("Ih_VS_p","Ih_VS_p",50,0,5,100,2,10);

  TH1F *Sigdigi_left_1=new TH1F("Sigdigi_left_1","Sigdigi_left_1",50,0,50);
  TH1F *Sigdigi_left_2=new TH1F("Sigdigi_left_2","Sigdigi_left_2",50,0,50);
  TH1F *Sigdigi_right_1=new TH1F("Sigdigi_right_1","Sigdigi_right_1",50,0,50);
  TH1F *Sigdigi_right_2=new TH1F("Sigdigi_right_2","Sigdigi_right_2",50,0,50);
  TH1F *Sigdigi_center_1=new TH1F("Sigdigi_center_1","Sigdigi_center_1",50,0,50);
  TH1F *Sigdigi_center_2=new TH1F("Sigdigi_center_2","Sigdigi_center_2",50,0,50);
  TH1F *Sigdigi_FullLeft=new TH1F("Sigdigi_FullLeft","Sigdigi_FullLeft",50,0,50);
  TH1F *Sigdigi_FullRight=new TH1F("Sigdigi_FullRight","Sigdigi_FullRight",50,0,50);
  TH1F *Sigdigi_FullLeft_channel=new TH1F("Sigdigi_FullLeft_channel","Sigdigi_FullLeft_channel",768,0,768);
  TH1F *Sigdigi_FullRight_channel=new TH1F("Sigdigi_FullRight_channel","Sigdigi_FullRight_channel",768,0,768);
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

  int Sum_left_1=0, Sum_left_2=0;
  int Sum_right_1=0, Sum_right_2=0;
  int Sum_center_1=0, Sum_center_2=0;
  int Sum_FullLeft=0;
  int Sum_FullRight=0;
  int cpt_left_1=0, cpt_left_2=0;
  int cpt_right_1=0, cpt_right_2=0;
  int cpt_center_1=0, cpt_center_2=0;
  int cpt_FullLeft=0;
  int cpt_FullRight=0;
  int cpt_2sat=0;
  int cpt_channel128_FullLeft=0, cpt_channel128_FullRight=0;
  int cpt_left_1_idem=0, cpt_left_2_idem=0;
  int cpt_right_1_idem=0, cpt_right_2_idem=0;
  int cpt_center_1_idem=0, cpt_center_2_idem=0;
  int cpt_FullLeft_idem=0, cpt_FullRight_idem=0;

  vector <int> cpt_cout_left_1;
  vector <int> cpt_cout_left_2;
  vector <int> cpt_cout_right_1;
  vector <int> cpt_cout_right_2;
  vector <int> cpt_cout_center_1;
  vector <int> cpt_cout_center_2;
  vector <int> cpt_cout_FullLeft;
  vector <int> cpt_cout_FullRight;
  vector <int> osef;
  vector <int> cpt_cout_all;
  for (int i=0; i<=10; i++)
  {
    cpt_cout_left_1.push_back(0);
    cpt_cout_left_2.push_back(0);
    cpt_cout_right_1.push_back(0);
    cpt_cout_right_2.push_back(0);
    cpt_cout_center_1.push_back(0);
    cpt_cout_center_2.push_back(0);
    cpt_cout_FullLeft.push_back(0);
    cpt_cout_FullRight.push_back(0);
    osef.push_back(0);
    cpt_cout_all.push_back(0);
  }
  

  double DeltaR=0;
  double DeltaR_val;
  int index;


  //nentries=1000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);   nbytes += nb;

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

    if (track_nvalidhits[index]>8 && track_qual[index])
    {
      for (int idedx=track_index_hit[index]; idedx<track_index_hit[index]+track_nhits[index]; idedx++)
      {
        if (dedx_isstrip[idedx])
        {
          double sum_sigdigi=0;
          double Qloss=(double)sclus_eloss[idedx]/(3.61*pow(10,-9)*265);
          vector<int> index_sig;
          vector<int> value_sig;
          vector<int> value_sig_nosat;
          double value_max=-1;
          int i_max=-1;
          bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
          bool Scnd_max_left=false, Scnd_max_right=false;

          for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          {
            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
            {
              index_sig.push_back(isigdigi);
              value_sig.push_back(sigdigi_strip_adc[isigdigi]);

              if (sigdigi_strip_adc[isigdigi]<=255) value_sig_nosat.push_back(sigdigi_strip_adc[isigdigi]);
              else if (sigdigi_strip_adc[isigdigi]>255 && sigdigi_strip_adc[isigdigi]<=512) value_sig_nosat.push_back(254);
              else if (sigdigi_strip_adc[isigdigi]>512) value_sig_nosat.push_back(255);
            }
          }

          if (!index_sig.empty() && !value_sig.empty())
          {
            value_max = *max_element(value_sig.begin(), value_sig.end());
            i_max = index_sig.front()+getIndex(value_sig,value_max);
            sum_sigdigi = accumulate(value_sig.begin(), value_sig.end(), 0);

            Compare_eloss_sigdigi->Fill(Qloss/sum_sigdigi);

            int check_TwoMaxSat=0;
            for (int i=1; i<value_sig.size(); i++)
            {
              if (value_sig[i-1]>255 && value_sig[i]>255) check_TwoMaxSat++;
            }


            if (check_TwoMaxSat==1 && value_sig[i_max-index_sig.front()-1]>255) Scnd_max_left=true;
            else if (check_TwoMaxSat==1 && value_sig[i_max-index_sig.front()+1]>255) Scnd_max_right=true;

            if (index_sig.size()>=4 && check_TwoMaxSat==1 
            && ((Scnd_max_left && i_max>index_sig[1] && i_max<index_sig[index_sig.size()-1] && sigdigi_strip_adc[i_max-2]>1.1*sigdigi_strip_adc[i_max+1])
            || (Scnd_max_right && i_max>index_sig[0] && i_max<index_sig[index_sig.size()-2] && sigdigi_strip_adc[i_max-1]>1.1*sigdigi_strip_adc[i_max+2])) 
            ) left=true;
            if (index_sig.size()>=4 && check_TwoMaxSat==1 
            && ((Scnd_max_left && i_max>index_sig[1] && i_max<index_sig[index_sig.size()-1] && sigdigi_strip_adc[i_max+1]>1.1*sigdigi_strip_adc[i_max-2])
            || (Scnd_max_right && i_max>index_sig[0] && i_max<index_sig[index_sig.size()-2] && sigdigi_strip_adc[i_max+2]>1.1*sigdigi_strip_adc[i_max-1])) 
            ) right=true;
            if (index_sig.size()>=4 && check_TwoMaxSat==1 
            && ((Scnd_max_left && i_max>index_sig[1] && i_max<index_sig[index_sig.size()-1] && sigdigi_strip_adc[i_max-2]<=1.1*sigdigi_strip_adc[i_max+1] && sigdigi_strip_adc[i_max-2]>=0.9*sigdigi_strip_adc[i_max+1])
            || (Scnd_max_right && i_max>index_sig[0] && i_max<index_sig[index_sig.size()-2] && sigdigi_strip_adc[i_max-1]<=1.1*sigdigi_strip_adc[i_max+2] && sigdigi_strip_adc[i_max-1]>=0.9*sigdigi_strip_adc[i_max+2])) 
            ) center=true;
            if (index_sig.size()>=3 && check_TwoMaxSat==1 && i_max==index_sig.front() && value_sig[i_max-index_sig.front()+1]>255) FullLeft=true;
            if (index_sig.size()>=3 && check_TwoMaxSat==1 && i_max==index_sig.back() && value_sig[i_max-index_sig.front()-1]>255) FullRight=true;
            
            if (left || right || center || FullLeft || FullRight) cpt_2sat++;

            if (left && Scnd_max_left)
            {
              Sum_left_1+=sum_sigdigi;
              cpt_left_1++;
              if (value_sig[i_max-index_sig.front()+1]>255) cout<<"probleme L left "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_left_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left_1->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_left_1->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_left_1->GetBinContent(25));
                if (i>i_max) Sigdigi_left_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left_1->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QLeft_L = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_left_1);
              vector <int> QLeft_L = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QLeft_L = accumulate(QLeft_L.begin(), QLeft_L.end(), 0);
              if (compareVectors(QLeft_L,value_sig_nosat)) cpt_left_1_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QLeft_L/sum_sigdigi);
                Compare_QII_Over_SigDigi_L->Fill((float)Sum_QLeft_L/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QLeft_L-sum_sigdigi);
                Compare_QII_Minus_SigDigi_L->Fill(Sum_QLeft_L-sum_sigdigi);
              }
            }
            if (left && Scnd_max_right)
            {
              Sum_left_2+=sum_sigdigi;
              cpt_left_2++;
              if (value_sig[i_max-index_sig.front()-1]>255) cout<<"probleme L right "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_left_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left_2->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_left_2->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_left_2->GetBinContent(25));
                if (i>i_max) Sigdigi_left_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_left_2->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QLeft_R = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_left_2);
              vector <int> QLeft_R = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QLeft_R = accumulate(QLeft_R.begin(), QLeft_R.end(), 0);
              if (compareVectors(QLeft_R,value_sig_nosat)) cpt_left_2_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QLeft_R/sum_sigdigi);
                Compare_QII_Over_SigDigi_L->Fill((float)Sum_QLeft_R/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QLeft_R-sum_sigdigi);
                Compare_QII_Minus_SigDigi_L->Fill(Sum_QLeft_R-sum_sigdigi);
              }
            }

            if (right && Scnd_max_left)
            {
              Sum_right_1+=sum_sigdigi;
              cpt_right_1++;
              if (value_sig[i_max-index_sig.front()+1]>255) cout<<"probleme R left "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_right_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right_1->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_right_1->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_right_1->GetBinContent(25));
                if (i>i_max) Sigdigi_right_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right_1->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QRight_L = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_right_1);
              vector <int> QRight_L = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QRight_L = accumulate(QRight_L.begin(), QRight_L.end(), 0);
              if (compareVectors(QRight_L,value_sig_nosat)) cpt_right_1_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QRight_L/sum_sigdigi);
                Compare_QII_Over_SigDigi_R->Fill((float)Sum_QRight_L/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QRight_L-sum_sigdigi);
                Compare_QII_Minus_SigDigi_R->Fill(Sum_QRight_L-sum_sigdigi);
              }
            }

            if (right && Scnd_max_right)
            {
              Sum_right_2+=sum_sigdigi;
              cpt_right_2++;
              if (value_sig[i_max-index_sig.front()-1]>255) cout<<"probleme R right "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_right_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right_2->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_right_2->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_right_2->GetBinContent(25));
                if (i>i_max) Sigdigi_right_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_right_2->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QRight_R = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_right_2);
              vector <int> QRight_R = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QRight_R = accumulate(QRight_R.begin(), QRight_R.end(), 0);
              if (compareVectors(QRight_R,value_sig_nosat)) cpt_right_2_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QRight_R/sum_sigdigi);
                Compare_QII_Over_SigDigi_R->Fill((float)Sum_QRight_R/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QRight_R-sum_sigdigi);
                Compare_QII_Minus_SigDigi_R->Fill(Sum_QRight_R-sum_sigdigi);
              }
            }

            if (center && Scnd_max_left)
            {
              Sum_center_1+=sum_sigdigi;
              cpt_center_1++;
              if (value_sig[i_max-index_sig.front()+1]>255) cout<<"probleme C left "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_center_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center_1->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_center_1->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_center_1->GetBinContent(25));
                if (i>i_max) Sigdigi_center_1->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center_1->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QCenter_L = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_center_1);
              vector <int> QCenter_L = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QCenter_L = accumulate(QCenter_L.begin(), QCenter_L.end(), 0);
              if (compareVectors(QCenter_L,value_sig_nosat)) cpt_center_1_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QCenter_L/sum_sigdigi);
                Compare_QII_Over_SigDigi_C->Fill((float)Sum_QCenter_L/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QCenter_L-sum_sigdigi);
                Compare_QII_Minus_SigDigi_C->Fill(Sum_QCenter_L-sum_sigdigi);
              }
            }

            if (center && Scnd_max_right)
            {
              Sum_center_2+=sum_sigdigi;
              cpt_center_2++;
              if (value_sig[i_max-index_sig.front()-1]>255) cout<<"probleme C right "<<jentry<<" "<<idedx<<" "<<"  index "<<" "<<index_sig.front()<<" "<<i_max<<" "<<index_sig.back()<<"  value "<<value_sig[i_max-index_sig.front()-1]<<" "<<value_sig[i_max-index_sig.front()]<<" "<<value_sig[i_max-index_sig.front()+1]<<endl;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_center_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center_2->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_center_2->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_center_2->GetBinContent(25));
                if (i>i_max) Sigdigi_center_2->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_center_2->GetBinContent(25+(i-i_max)));
              }

              //vector<int> QCenter_R = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_center_2);
              vector <int> QCenter_R = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QCenter_R = accumulate(QCenter_R.begin(), QCenter_R.end(), 0);
              if (compareVectors(QCenter_R,value_sig_nosat)) cpt_center_2_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QCenter_R/sum_sigdigi);
                Compare_QII_Over_SigDigi_C->Fill((float)Sum_QCenter_R/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QCenter_R-sum_sigdigi);
                Compare_QII_Minus_SigDigi_C->Fill(Sum_QCenter_R-sum_sigdigi);
              }
            }

            if (FullLeft)
            {
              Sum_FullLeft+=sum_sigdigi;
              cpt_FullLeft++;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_FullLeft->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_FullLeft->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullLeft->GetBinContent(25));
                if (i>i_max) Sigdigi_FullLeft->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullLeft->GetBinContent(25+(i-i_max)));
              }

              Sigdigi_FullLeft_channel->Fill(sigdigi_strip_channel[i_max]);
              if (sigdigi_strip_channel[i_max]%128==0) cpt_channel128_FullLeft++;

              //vector<int> QFullLeft = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_FullLeft);
              vector <int> QFullLeft = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QFullLeft = accumulate(QFullLeft.begin(), QFullLeft.end(), 0);
              if (compareVectors(QFullLeft,value_sig_nosat)) cpt_FullLeft_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QFullLeft/sum_sigdigi);
                Compare_QII_Over_SigDigi_FL->Fill((float)Sum_QFullLeft/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QFullLeft-sum_sigdigi);
                Compare_QII_Minus_SigDigi_FL->Fill(Sum_QFullLeft-sum_sigdigi);
              }
            }

            if (FullRight)
            {
              Sum_FullRight+=sum_sigdigi;
              cpt_FullRight++;
              for (int i=index_sig.front(); i<=index_sig.back(); i++)
              {
                if (i<i_max) Sigdigi_FullRight->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight->GetBinContent(25+(i-i_max)));
                if (i==i_max) Sigdigi_FullRight->SetBinContent(25,sigdigi_strip_adc[i_max]+Sigdigi_FullRight->GetBinContent(25));
                if (i>i_max) Sigdigi_FullRight->SetBinContent(25+(i-i_max),sigdigi_strip_adc[i]+Sigdigi_FullRight->GetBinContent(25+(i-i_max)));
              }

              Sigdigi_FullRight_channel->Fill(sigdigi_strip_channel[i_max]);
              if ((sigdigi_strip_channel[i_max]+1)%128==0) cpt_channel128_FullRight++;

              //vector<int> QFullRight = SaturationCorrection(value_sig_nosat,0.10,0.04, true,20,25, cpt_cout_FullRight);
              vector <int> QFullRight = New_SaturationCorrection(value_sig_nosat, true, 25, cpt_cout_all);
              int Sum_QFullRight = accumulate(QFullRight.begin(), QFullRight.end(), 0);
              if (compareVectors(QFullRight,value_sig_nosat)) cpt_FullRight_idem++;
              else
              {
                Compare_QII_Over_SigDigi->Fill((float)Sum_QFullRight/sum_sigdigi);
                Compare_QII_Over_SigDigi_FR->Fill((float)Sum_QFullRight/sum_sigdigi);
                Compare_QII_Minus_SigDigi->Fill(Sum_QFullRight-sum_sigdigi);
                Compare_QII_Minus_SigDigi_FR->Fill(Sum_QFullRight-sum_sigdigi);
              }
            }
            
          } //end !vector.empty()
        } //end dedx_isstrip
      } //end idedx
    }

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
        double Ih= getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, NULL,2, 0., 0, nval, nsat);

        Ih_VS_p->Fill(track_p[itr],Ih);
      }

    }//end tracks

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<"%)"<<endl;

  } //end EVENTS

  Sigdigi_left_1->Scale(1./Sum_left_1);
  Sigdigi_left_2->Scale(1./Sum_left_2);
  Sigdigi_right_1->Scale(1./Sum_right_1);
  Sigdigi_right_2->Scale(1./Sum_right_2);
  Sigdigi_center_1->Scale(1./Sum_center_1);
  Sigdigi_center_2->Scale(1./Sum_center_2);
  Sigdigi_FullLeft->Scale(1./Sum_FullLeft);
  Sigdigi_FullRight->Scale(1./Sum_FullRight);
  
  cpt_sat->SetBinContent(1,(float)cpt_left_1/(float)Sum_left_1);
  cpt_sat->SetBinContent(2,(float)cpt_left_2/(float)Sum_left_2);
  cpt_sat->SetBinContent(3,(float)cpt_right_1/(float)Sum_right_1);
  cpt_sat->SetBinContent(4,(float)cpt_right_2/(float)Sum_right_2);
  cpt_sat->SetBinContent(5,(float)cpt_center_1/(float)Sum_center_1);
  cpt_sat->SetBinContent(6,(float)cpt_center_2/(float)Sum_center_2);
  cpt_sat->SetBinContent(7,(float)cpt_FullLeft/(float)Sum_FullLeft);
  cpt_sat->SetBinContent(8,(float)cpt_FullRight/(float)Sum_FullRight);
  

  SaveData->cd();
  Compare_eloss_sigdigi->Write();
  Ih_VS_p->Write();
  Sigdigi_left_1->Write();
  Sigdigi_left_2->Write();
  Sigdigi_right_1->Write();
  Sigdigi_right_2->Write();
  Sigdigi_center_1->Write();
  Sigdigi_center_2->Write();
  Sigdigi_FullLeft->Write();
  Sigdigi_FullRight->Write();
  Sigdigi_FullLeft_channel->Write();
  Sigdigi_FullRight_channel->Write();
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

  cout<<"2 strips sat: "<<cpt_2sat<<endl;
  cout<<"left L: "<<cpt_left_1<<endl;
  cout<<"left R: "<<cpt_left_2<<endl;
  cout<<"right L: "<<cpt_right_1<<endl;
  cout<<"right R: "<<cpt_right_2<<endl;
  cout<<"center L: "<<cpt_center_1<<endl;
  cout<<"center R: "<<cpt_center_2<<endl;
  cout<<"FullLeft: "<<cpt_FullLeft<<endl;
  cout<<"FullRight: "<<cpt_FullRight<<endl;

  cout<<"Channel % 128 == 0 FullLeft: "<<cpt_channel128_FullLeft<<endl;
  cout<<"Channel+1 % 128 == 0 FullRight: "<<cpt_channel128_FullRight<<endl;
/*
  cout<<"         LEFT L"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_left_1_idem<<" / "<<cpt_left_1<<"  ("<<100*(double)cpt_left_1_idem/cpt_left_1<<" %)"<<endl;
  for (int i=0; i<cpt_cout_left_1.size(); i++) cout<<"cpt_cout_left_1["<<i<<"]  "<<cpt_cout_left_1[i]<<"   "<<100*(double)cpt_cout_left_1[i]/cpt_left_1_idem<<" %"<<endl;
  cout<<"         LEFT R"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_left_2_idem<<" / "<<cpt_left_2<<"  ("<<100*(double)cpt_left_2_idem/cpt_left_2<<" %)"<<endl;
  for (int i=0; i<cpt_cout_left_2.size(); i++) cout<<"cpt_cout_left_2["<<i<<"]  "<<cpt_cout_left_2[i]<<"   "<<100*(double)cpt_cout_left_2[i]/cpt_left_2_idem<<" %"<<endl;
  cout<<"         RIGHT L"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_right_1_idem<<" / "<<cpt_right_1<<"  ("<<100*(double)cpt_right_1_idem/cpt_right_1<<" %)"<<endl;
  for (int i=0; i<cpt_cout_right_1.size(); i++) cout<<"cpt_cout_right_1["<<i<<"]  "<<cpt_cout_right_1[i]<<"   "<<100*(double)cpt_cout_right_1[i]/cpt_right_1_idem<<" %"<<endl;
  cout<<"         RIGHT R"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_right_2_idem<<" / "<<cpt_right_2<<"  ("<<100*(double)cpt_right_2_idem/cpt_right_2<<" %)"<<endl;
  for (int i=0; i<cpt_cout_right_2.size(); i++) cout<<"cpt_cout_right_2["<<i<<"]  "<<cpt_cout_right_2[i]<<"   "<<100*(double)cpt_cout_right_2[i]/cpt_right_2_idem<<" %"<<endl;  
  cout<<"         CENTER L"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_center_1_idem<<" / "<<cpt_center_1<<"  ("<<100*(double)cpt_center_1_idem/cpt_center_1<<" %)"<<endl;
  for (int i=0; i<cpt_cout_center_1.size(); i++) cout<<"cpt_cout_center_1["<<i<<"]  "<<cpt_cout_center_1[i]<<"   "<<100*(double)cpt_cout_center_1[i]/cpt_center_1_idem<<" %"<<endl;
  cout<<"         CENTER R"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_center_2_idem<<" / "<<cpt_center_2<<"  ("<<100*(double)cpt_center_2_idem/cpt_center_2<<" %)"<<endl;
  for (int i=0; i<cpt_cout_center_2.size(); i++) cout<<"cpt_cout_center_2["<<i<<"]  "<<cpt_cout_center_2[i]<<"   "<<100*(double)cpt_cout_center_2[i]/cpt_center_2_idem<<" %"<<endl;
  cout<<"         FULL LEFT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_FullLeft_idem<<" / "<<cpt_FullLeft<<"  ("<<100*(double)cpt_FullLeft_idem/cpt_FullLeft<<" %)"<<endl;
  for (int i=0; i<cpt_cout_FullLeft.size(); i++) cout<<"cpt_cout_FullLeft["<<i<<"]  "<<cpt_cout_FullLeft[i]<<"   "<<100*(double)cpt_cout_FullLeft[i]/cpt_FullLeft_idem<<" %"<<endl;
  cout<<"         FULL RIGHT"<<endl;
  cout<<"nbr of clusters rejected  "<<cpt_FullRight_idem<<" / "<<cpt_FullRight<<"  ("<<100*(double)cpt_FullRight_idem/cpt_FullRight<<" %)"<<endl;
  for (int i=0; i<cpt_cout_FullRight.size(); i++) cout<<"cpt_cout_FullRight["<<i<<"]  "<<cpt_cout_FullRight[i]<<"   "<<100*(double)cpt_cout_FullRight[i]/cpt_FullRight_idem<<" %"<<endl;  
*/

  cout<<"         NBR OF CLUSTERS RECONSTRUCTED"<<endl;
  cout<<"nbr of clusters rejected  "<<(cpt_left_1_idem+cpt_left_2_idem+cpt_right_1_idem+cpt_right_2_idem+cpt_center_1_idem+cpt_center_2_idem+cpt_FullLeft_idem+cpt_FullRight_idem)<<" / "<<(cpt_left_1+cpt_left_2+cpt_right_1+cpt_right_2+cpt_center_1+cpt_center_2+cpt_FullLeft+cpt_FullRight)<<"  ("<<100*(double)(cpt_left_1_idem+cpt_left_2_idem+cpt_right_1_idem+cpt_right_2_idem+cpt_center_1_idem+cpt_center_2_idem+cpt_FullLeft_idem+cpt_FullRight_idem)/(cpt_left_1+cpt_left_2+cpt_right_1+cpt_right_2+cpt_center_1+cpt_center_2+cpt_FullLeft+cpt_FullRight)<<" %)"<<endl;
  for (int i=0; i<cpt_cout_all.size(); i++) cout<<"cpt_cout_all["<<i<<"]  "<<cpt_cout_all[i]<<endl;

}
