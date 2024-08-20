#define sat_observation_Wjets_cxx
#include "sat_observation_Wjets.h"
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
  else if(subdetid==4)  //TID
  {    
    if(((detid>>11)&0x3)==1) return 11;
    else if(((detid>>11)&0x3)==2) return 12;
    else if(((detid>>11)&0x3)==3) return 13;
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

std::vector<int> SaturationCorrection(const std::vector <int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  TMatrix A(N,N);

//---  que pour 1 max bien net
 if(Q.size()<2 || Q.size()>8){
        for (unsigned int i=0;i<Q.size();i++){
                QII.push_back((int) Q[i]);
        }
        return QII;
  }
 if(way){
          vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
          if(*mQ>253){
                 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
                 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
                     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
          }
      else{
          return Q; // no saturation --> no x-talk inversion
      }
  }
//---
 // do nothing else
 return Q;
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


void sat_observation_Wjets::Loop()
{
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  TFile* SaveData = new TFile("ROOT_SVG/Sat_Wjets_alltracks.root", "RECREATE");
  ifstream CheckDeadStrip("ROOT_SVG/CheckDeadStrip.txt");

  TH1F *sat_254_noclean=new TH1F("sat_254_noclean","sat_254_noclean",10000,0,1);
  TH1F *sat_255_noclean=new TH1F("sat_255_noclean","sat_255_noclean",10000,0,1);
  TH1F *sat_TOT_noclean=new TH1F("sat_TOT_noclean","sat_TOT_noclean",10000,0,1);

  TH1F *sat_254_clean2=new TH1F("sat_254_clean2","sat_254_clean2",10000,0,1);
  TH1F *sat_255_clean2=new TH1F("sat_255_clean2","sat_255_clean2",10000,0,1);
  TH1F *sat_TOT_clean2=new TH1F("sat_TOT_clean2","sat_TOT_clean2",10000,0,1);

  TH1F *sat_254_clean1=new TH1F("sat_254_clean1","sat_254_clean1",10000,0,1);
  TH1F *sat_255_clean1=new TH1F("sat_255_clean1","sat_255_clean1",10000,0,1);
  TH1F *sat_TOT_clean1=new TH1F("sat_TOT_clean1","sat_TOT_clean1",10000,0,1);

  TH1F *sat_254_newclean=new TH1F("sat_254_newclean","sat_254_newclean",10000,0,1);
  TH1F *sat_255_newclean=new TH1F("sat_255_newclean","sat_255_newclean",10000,0,1);
  TH1F *sat_TOT_newclean=new TH1F("sat_TOT_newclean","sat_TOT_newclean",10000,0,1);

  TH2F *sat_254_noclean_VS_p=new TH2F("sat_254_noclean_VS_p","sat_254_noclean_VS_p",300,0,300,100,0,1);
  TH2F *sat_255_noclean_VS_p=new TH2F("sat_255_noclean_VS_p","sat_255_noclean_VS_p",300,0,300,100,0,1);

  TH2F *sat_TOT_noclean_VS_p=new TH2F("sat_TOT_noclean_VS_p","sat_TOT_noclean_VS_p",300,0,300,100,0,1);
  TH2F *sat_TOT_noclean_VS_pt=new TH2F("sat_TOT_noclean_VS_pt","sat_TOT_noclean_VS_pt",300,0,300,100,0,1);
  TH2F *sat_TOT_noclean_VS_eta=new TH2F("sat_TOT_noclean_VS_eta","sat_TOT_noclean_VS_eta",100,-6,6,100,0,1);
  TH2F *sat_TOT_noclean_VS_p_COUNT=new TH2F("sat_TOT_noclean_VS_p_COUNT","sat_TOT_noclean_VS_p_COUNT",2000,0,2000,30,0,30);
  TH2F *sat_TOT_noclean_VS_pt_COUNT=new TH2F("sat_TOT_noclean_VS_pt_COUNT","sat_TOT_noclean_VS_pt_COUNT",300,0,300,30,0,30);
  TH2F *sat_TOT_noclean_VS_eta_COUNT=new TH2F("sat_TOT_noclean_VS_eta_COUNT","sat_TOT_noclean_VS_eta_COUNT",100,-6,6,30,0,30);
  

  TH2F *sat_254_clean2_VS_p=new TH2F("sat_254_clean2_VS_p","sat_254_clean2_VS_p",30,0,300,10000,0,1);
  TH2F *sat_255_clean2_VS_p=new TH2F("sat_255_clean2_VS_p","sat_255_clean2_VS_p",30,0,300,10000,0,1);
  TH2F *sat_TOT_clean2_VS_p=new TH2F("sat_TOT_clean2_VS_p","sat_TOT_clean2_VS_p",30,0,300,10000,0,1);

  TH2F *sat_254_clean1_VS_p=new TH2F("sat_254_clean1_VS_p","sat_254_clean1_VS_p",30,0,300,10000,0,1);
  TH2F *sat_255_clean1_VS_p=new TH2F("sat_255_clean1_VS_p","sat_255_clean1_VS_p",30,0,300,10000,0,1);
  TH2F *sat_TOT_clean1_VS_p=new TH2F("sat_TOT_clean1_VS_p","sat_TOT_clean1_VS_p",30,0,300,10000,0,1);

  TH2F *sat_254_newclean_VS_p=new TH2F("sat_254_newclean_VS_p","sat_254_newclean_VS_p",30,0,300,10000,0,1);
  TH2F *sat_255_newclean_VS_p=new TH2F("sat_255_newclean_VS_p","sat_255_newclean_VS_p",30,0,300,10000,0,1);
  TH2F *sat_TOT_newclean_VS_p=new TH2F("sat_TOT_newclean_VS_p","sat_TOT_newclean_VS_p",30,0,300,10000,0,1);

  TH1F *check_pt_matched_muon=new TH1F("check_pt_matched_muon","check_pt_matched_muon",1000,0,1000);
  TH1F *check_eta_matched_muon=new TH1F("check_eta_matched_muon","check_eta_matched_muon",100,-6,6);
  TH1F *check_phi_matched_muon=new TH1F("check_phi_matched_muon","check_phi_matched_muon",100,-7,7);

  TH2F *pathlength_VS_eta=new TH2F("pathlength_VS_eta","pathlength_VS_eta",100,0,0.1,100,-6,6);
  TH2F *clusterTrack_VS_eta=new TH2F("sat_clusterTrack_VS_eta","sat_clusterTrack_VS_eta",100,-6,6,30,0,30);
  
  TH2F *pt_VS_eta_sat=new TH2F("pt_VS_eta_sat","pt_VS_eta_sat",300,0,300,100,-6,6);
  TH2F *p_VS_eta_sat=new TH2F("p_VS_eta_sat","p_VS_eta_sat",2000,0,2000,100,-6,6);
  TH2F *pt_VS_eta_tot=new TH2F("pt_VS_eta_tot","pt_VS_eta_tot",300,0,300,100,-6,6);
  TH2F *p_VS_eta_tot=new TH2F("p_VS_eta_tot","p_VS_eta_tot",2000,0,2000,100,-6,6);
  
  TH2F *pt_VS_eta=new TH2F("pt_VS_eta","pt_VS_eta",300,0,300,100,-6,6);
  TH2F *p_VS_eta=new TH2F("p_VS_eta","p_VS_eta",2000,0,2000,100,-6,6);


  TH2F *Ih_VS_p_1=new TH2F("Ih_VS_p_1","Ih_VS_p_1",500,0,5,800,2,10);
  TH2F *Ih_VS_p_2=new TH2F("Ih_VS_p_2","Ih_VS_p_2",1000,0,100,800,2,10);
  TH2F *sat_cluster_1_VS_p=new TH2F("sat_cluster_1_VS_p","sat_cluster_1_VS_p",200,0,2000,10,0,10);
  TH2F *sat_cluster_1_VS_pt=new TH2F("sat_cluster_1_VS_pt","sat_cluster_1_VS_pt",300,0,300,10,0,10);
  TH2F *sat_cluster_2_VS_p=new TH2F("sat_cluster_2_VS_p","sat_cluster_2_VS_p",200,0,2000,10,0,10);
  TH2F *sat_cluster_2_VS_pt=new TH2F("sat_cluster_2_VS_pt","sat_cluster_2_VS_pt",300,0,300,10,0,10);

  TH1D *sat_TOT_VS_TECrings=new TH1D("sat_TOT_VS_TECrings","sat_TOT_VS_TECrings",20,0,20);
  TH1D *Cluster_layer_sat=new TH1D("Cluster_layer_sat","Cluster_layer_sat",35,0,35);
  TH1D *Cluster_layer_tot=new TH1D("Cluster_layer_tot","Cluster_layer_tot",35,0,35);
  TH1D *ClusterPerTrack_layer_tot=new TH1D("ClusterPerTrack_layer_tot","ClusterPerTrack_layer_tot",35,0,35);
  

  vector<int> NotDeadDetid; 
  double DeltaR=0;
  double DeltaR_val;
  int index;

  int cpt_ring1=0, cpt_ring2=0, cpt_ring3=0, cpt_ring4=0, cpt_ring5=0, cpt_ring6=0, cpt_ring7=0;
  int cpt_cluster_TEC1=0, cpt_cluster_TEC2=0, cpt_cluster_TEC3=0, cpt_cluster_TEC4=0, cpt_cluster_TEC5=0, cpt_cluster_TEC6=0, cpt_cluster_TEC7=0;
    
  int cpt_track=0;


  // Counter of all clusters, uncomment here and comment bellow
  // int cpt_sat_254_noclean=0, cpt_sat_255_noclean=0, cpt_sat_TOT_noclean=0, cpt_cluster_strip_noclean=0;
  // int cpt_sat_254_clean2=0, cpt_sat_255_clean2=0, cpt_sat_TOT_clean2=0, cpt_cluster_strip_clean2=0;
  // int cpt_sat_254_clean1=0, cpt_sat_255_clean1=0, cpt_sat_TOT_clean1=0, cpt_cluster_strip_clean1=0;
  // int cpt_sat_254_newclean=0, cpt_sat_255_newclean=0, cpt_sat_TOT_newclean=0, cpt_cluster_strip_newclean=0;

  nentries=100000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

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

  // Check dead strips
  string line;
  while (std::getline(CheckDeadStrip, line))
  {
    std::istringstream iss(line);
    int w_detid=-1, idx_dead_strip=-1;
    if (!(iss >> w_detid >> idx_dead_strip)) cout<<"ERROR: fill CheckDeadStrip"<<endl; // error
    
    for (int itr=0; itr<ntracks; itr++)
    { 
      if (track_nvalidhits[itr]>8 && track_qual[itr])
      {
        for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
        {
          if (dedx_isstrip[idedx])
          {
            if (w_detid==dedx_detid[idedx] && idx_dead_strip>=0  && strip_ampl[idx_dead_strip]>0 && idx_dead_strip<768)
            {
              NotDeadDetid.push_back(dedx_detid[idedx]);
            }
          }
        }
      }
    }
  }


  // Ih VS p
  for (int itr=0; itr<ntracks; itr++)
  {
    std::vector <int> ClusterCharge;
    bool ClusterSat = false;
    float Ih=0;

    for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
    {
      for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
      {
        ClusterCharge.push_back(strip_ampl[istrip]);
        if (strip_ampl[istrip]>=254) ClusterSat = true;
      }

      int ClusterChargeSUM = accumulate(ClusterCharge.begin(), ClusterCharge.end(), 0);
      Ih += 1./Power(ClusterChargeSUM,2);
    }

    Ih = 1./Sqrt(Ih/track_nhits[itr]);

    if (ClusterSat == false) Ih_VS_p_NoSat->Fill(track_p[itr], Ih, track_prescale[itr]);

    
    Ih_VS_p_SatNoCorr->Fill(track_p[itr], Ih, track_prescale[itr]);
    Ih_VS_p_SatCorr->Fill(track_p[itr], Ih, track_prescale[itr]);
  }



  // Clusters analysis
  if (track_nvalidhits[index]>8 && track_qual[index] && index>-1 && DeltaR_val<0.15)
  {
    check_pt_matched_muon->Fill(track_pt[index]);
    check_eta_matched_muon->Fill(track_eta[index]);
    check_phi_matched_muon->Fill(track_phi[index]);
  }

  for (int itr=0; itr<ntracks; itr++)
  {
    vector <float> charge_corr;
    vector <float> pathlength;
    vector <int> subdetId;
    vector <int> moduleGeometry;
    vector <bool> bool_cleaning;
    vector <bool> mustBeInside;  

    if (track_nvalidhits[itr]>8 && track_qual[itr]) // change "itr" by "index"
    {
      int cpt_sat_254_noclean=0, cpt_sat_255_noclean=0, cpt_sat_TOT_noclean=0, cpt_cluster_strip_noclean=0;
      int cpt_sat_254_clean2=0, cpt_sat_255_clean2=0, cpt_sat_TOT_clean2=0, cpt_cluster_strip_clean2=0;
      int cpt_sat_254_clean1=0, cpt_sat_255_clean1=0, cpt_sat_TOT_clean1=0, cpt_cluster_strip_clean1=0;
      int cpt_sat_254_newclean=0, cpt_sat_255_newclean=0, cpt_sat_TOT_newclean=0, cpt_cluster_strip_newclean=0;
      
      cpt_track++;

      for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
      {
        //Clusclean part
        vector <int> ampl;
        vector <int> ampl_corr;

        if (dedx_isstrip[idedx] && sclus_nstrip_corr[idedx]>1)
        {
          // Ih -------
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
            
            vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
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

          // Ih end ------

          // CLEANING
          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            ampl.push_back(strip_ampl[istrip]);
          }
          for (int istrip2=sclus_index_strip_corr[idedx]; istrip2<sclus_index_strip_corr[idedx]+sclus_nstrip_corr[idedx]; istrip2++)
          {
            ampl_corr.push_back(strip_ampl_corr[istrip2]);
          }

          bool clean_UNCHANGED=clusterCleaning_UNCHANGED(ampl,0,0,0.17,0.04,4);
          bool clean=clusterCleaning(ampl,0,0,0.17,0.04,4);
          bool clean2=clusterCleaning_UNCHANGED(ampl_corr,1,0,0.17,0.04,4);
          //bool clean_noDoublePic=clusterCleaning_UNCHANGED(ampl,0,0,10000,10000,10000);
          bool newclean=myNewCleaning(ampl,0,0);

          // -------
          
          bool sat_254_bool=false, sat_255_bool=false;

          // GENERATED SIGNAL
          // for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
          // {
          //   if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
          //   {
          //     if (sigdigi_strip_adc[isigdigi]>=254 && sigdigi_strip_adc[isigdigi]<512) sat_254_bool=true;
          //     if (sigdigi_strip_adc[isigdigi]>=512) sat_255_bool=true;
          //   }
          // }

          // RECONSTRUCTED SIGNAL (comment the other and uncomment here)
          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
          {
            if (strip_ampl[istrip]==254) sat_254_bool=true;
            if (strip_ampl[istrip]==255) sat_255_bool=true;
          }
          

          if (sat_254_bool)
          {
            cpt_sat_254_noclean++;
            if (clean2) cpt_sat_254_clean2++;
            if (clean_UNCHANGED) cpt_sat_254_clean1++;
            if (newclean) cpt_sat_254_newclean++;
          }
          if (sat_255_bool)
          {
            cpt_sat_255_noclean++;
            if (clean2) cpt_sat_255_clean2++;
            if (clean_UNCHANGED) cpt_sat_255_clean1++;
            if (newclean) cpt_sat_255_newclean++;
          }

          cpt_sat_TOT_noclean = cpt_sat_255_noclean + cpt_sat_254_noclean;
          if (clean2) cpt_sat_TOT_clean2 = cpt_sat_255_clean2 + cpt_sat_254_clean2;
          if (clean_UNCHANGED) cpt_sat_TOT_clean1 = cpt_sat_255_clean1 + cpt_sat_254_clean1;
          if (newclean) cpt_sat_TOT_newclean = cpt_sat_255_newclean + cpt_sat_254_newclean;
          
          cpt_cluster_strip_noclean++;
          if (clean2) cpt_cluster_strip_clean2++;
          if (clean_UNCHANGED) cpt_cluster_strip_clean1++;
          if (newclean) cpt_cluster_strip_newclean++;


          if (dedx_layer[idedx]>=5 && dedx_layer[idedx]<=10) pathlength_VS_eta->Fill(dedx_pathlength[idedx], track_eta[itr], track_prescale[itr]);

          // Clusters saturated layer by layer
          if (sat_254_bool || sat_255_bool) Cluster_layer_sat->Fill(dedx_layer[idedx]);
          Cluster_layer_tot->Fill(dedx_layer[idedx]);
          ClusterPerTrack_layer_tot->Fill(dedx_layer[idedx]);


          // Count the number of saturated clusters per ring
          if (dedx_subdetid[idedx]==6)
          {
            int ring = FindLayer(6,dedx_detid[idedx]);
            if (sat_254_bool || sat_255_bool)
            {
              Cluster_layer_sat->Fill(23+ring);
              if (ring==1) cpt_ring1++;
              else if (ring==2) cpt_ring2++;
              else if (ring==3) cpt_ring3++;
              else if (ring==4) cpt_ring4++;
              else if (ring==5) cpt_ring5++;
              else if (ring==6) cpt_ring6++;
              else if (ring==7) cpt_ring7++;
              else cout<<"ERROR: ring number greater than 7"<<endl;
            }

            Cluster_layer_tot->Fill(23+ring);
            ClusterPerTrack_layer_tot->Fill(23+ring);
            if (ring==1) cpt_cluster_TEC1++;
            else if (ring==2) cpt_cluster_TEC2++;
            else if (ring==3) cpt_cluster_TEC3++;
            else if (ring==4) cpt_cluster_TEC4++;
            else if (ring==5) cpt_cluster_TEC5++;
            else if (ring==6) cpt_cluster_TEC6++;
            else if (ring==7) cpt_cluster_TEC7++;
            else cout<<"ERROR: ring number greater than 7"<<endl;
          }


        } //end dedx_isstrip
      } //end idedx

      int nval=0;
      int nsat=0;
      double Ih= getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, NULL,2, 0., 0, nval, nsat);

      if (Ih >= -5.5*track_p[itr]+10 && Ih >= 4 && track_p[itr]<=2)
      {
        Ih_VS_p_1->Fill(track_p[itr], Ih, track_prescale[itr]);
        sat_cluster_1_VS_p->Fill(track_p[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
        sat_cluster_1_VS_pt->Fill(track_pt[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
      }
      if (track_p[itr]>2 && Ih <= 4)
      {
        Ih_VS_p_2->Fill(track_p[itr], Ih, track_prescale[itr]);
        sat_cluster_2_VS_p->Fill(track_p[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
        sat_cluster_2_VS_pt->Fill(track_pt[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
      }

      //if (cpt_sat_254_noclean!=0 || cpt_sat_254_clean2!=0 || cpt_sat_254_clean1!=0 || cpt_sat_254_newclean!=0) cout<<"ok"<<endl;
      
      sat_254_noclean->Fill((float)cpt_sat_254_noclean/(float)cpt_cluster_strip_noclean);
      sat_255_noclean->Fill((float)cpt_sat_255_noclean/(float)cpt_cluster_strip_noclean);
      sat_TOT_noclean->Fill((float)cpt_sat_TOT_noclean/(float)cpt_cluster_strip_noclean);

      sat_254_clean2->Fill((float)cpt_sat_254_clean2/(float)cpt_cluster_strip_clean2);
      sat_255_clean2->Fill((float)cpt_sat_255_clean2/(float)cpt_cluster_strip_clean2);
      sat_TOT_clean2->Fill((float)cpt_sat_TOT_clean2/(float)cpt_cluster_strip_clean2);

      sat_254_clean1->Fill((float)cpt_sat_254_clean1/(float)cpt_cluster_strip_clean1);
      sat_255_clean1->Fill((float)cpt_sat_255_clean1/(float)cpt_cluster_strip_clean1);
      sat_TOT_clean1->Fill((float)cpt_sat_TOT_clean1/(float)cpt_cluster_strip_clean1);

      sat_254_newclean->Fill((float)cpt_sat_254_newclean/(float)cpt_cluster_strip_newclean);
      sat_255_newclean->Fill((float)cpt_sat_255_newclean/(float)cpt_cluster_strip_newclean);
      sat_TOT_newclean->Fill((float)cpt_sat_TOT_newclean/(float)cpt_cluster_strip_newclean);


      // With weight
      sat_254_noclean_VS_p->Fill(track_p[itr], (float)cpt_sat_254_noclean/(float)cpt_cluster_strip_noclean, track_prescale[itr]);
      sat_255_noclean_VS_p->Fill(track_p[itr], (float)cpt_sat_255_noclean/(float)cpt_cluster_strip_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_p->Fill(track_p[itr], (float)cpt_sat_TOT_noclean/(float)cpt_cluster_strip_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_pt->Fill(track_pt[itr], (float)cpt_sat_TOT_noclean/(float)cpt_cluster_strip_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_eta->Fill(track_eta[itr], (float)cpt_sat_TOT_noclean/(float)cpt_cluster_strip_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_p_COUNT->Fill(track_p[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_pt_COUNT->Fill(track_pt[itr], cpt_sat_TOT_noclean, track_prescale[itr]);
      sat_TOT_noclean_VS_eta_COUNT->Fill(track_eta[itr], cpt_sat_TOT_noclean, track_prescale[itr]);

      clusterTrack_VS_eta->Fill(track_eta[itr], track_nhits[itr], track_prescale[itr]);
      
      //cout<<"bin: ["<<track_p[itr]<<" ; "<<track_eta[itr]<<"]   value: "<<cpt_sat_TOT_noclean<<endl;
      int bin_eta = (int)((track_eta[itr]+6)*(100./12.));
      int bin_p = (track_p[itr]-0)*(2000/2000);
      int bin_pt = (track_pt[itr]-0)*(300/300);
      pt_VS_eta_sat->SetBinContent(bin_pt,bin_eta, cpt_sat_TOT_noclean+pt_VS_eta_sat->GetBinContent(bin_pt, bin_eta));
      p_VS_eta_sat->SetBinContent(bin_p,bin_eta, cpt_sat_TOT_noclean+p_VS_eta_sat->GetBinContent(bin_p, bin_eta));
      pt_VS_eta_tot->SetBinContent(bin_pt,bin_eta, cpt_cluster_strip_noclean+pt_VS_eta_tot->GetBinContent(bin_pt, bin_eta));
      p_VS_eta_tot->SetBinContent(bin_p,bin_eta, cpt_cluster_strip_noclean+p_VS_eta_tot->GetBinContent(bin_p, bin_eta));
      
      pt_VS_eta->Fill(track_pt[itr], track_eta[itr], track_prescale[itr]);
      p_VS_eta->Fill(track_p[itr], track_eta[itr], track_prescale[itr]);
      

      // Without
      sat_254_clean2_VS_p->Fill(track_p[itr],(float)cpt_sat_254_clean2/(float)cpt_cluster_strip_clean2);
      sat_255_clean2_VS_p->Fill(track_p[itr],(float)cpt_sat_255_clean2/(float)cpt_cluster_strip_clean2);
      sat_TOT_clean2_VS_p->Fill(track_p[itr],(float)cpt_sat_TOT_clean2/(float)cpt_cluster_strip_clean2);

      sat_254_clean1_VS_p->Fill(track_p[itr],(float)cpt_sat_254_clean1/(float)cpt_cluster_strip_clean1);
      sat_255_clean1_VS_p->Fill(track_p[itr],(float)cpt_sat_255_clean1/(float)cpt_cluster_strip_clean1);
      sat_TOT_clean1_VS_p->Fill(track_p[itr],(float)cpt_sat_TOT_clean1/(float)cpt_cluster_strip_clean1);

      sat_254_newclean_VS_p->Fill(track_p[itr],(float)cpt_sat_254_newclean/(float)cpt_cluster_strip_newclean);
      sat_255_newclean_VS_p->Fill(track_p[itr],(float)cpt_sat_255_newclean/(float)cpt_cluster_strip_newclean);
      sat_TOT_newclean_VS_p->Fill(track_p[itr],(float)cpt_sat_TOT_newclean/(float)cpt_cluster_strip_newclean);


      // Saturated Clusters in TEC rings
      /*if (cpt_cluster_strip_noclean>0 && (cpt_ring1>0||cpt_ring2>0||cpt_ring3>0||cpt_ring4>0||cpt_ring5>0||cpt_ring6>0||cpt_ring7>0))
      {
        cpt_norm++;
        
        //if (jentry%1000==0) cout<<(float)cpt_ring1/cpt_cluster_strip_noclean<<" "<<(float)cpt_ring4/cpt_cluster_strip_noclean<<endl;
        sat_TOT_VS_TECrings->SetBinContent(1,(float)cpt_ring1/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(1));
        sat_TOT_VS_TECrings->SetBinContent(2,(float)cpt_ring2/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(2));
        sat_TOT_VS_TECrings->SetBinContent(3,(float)cpt_ring3/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(3));
        sat_TOT_VS_TECrings->SetBinContent(4,(float)cpt_ring4/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(4));
        sat_TOT_VS_TECrings->SetBinContent(5,(float)cpt_ring5/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(5));
        sat_TOT_VS_TECrings->SetBinContent(6,(float)cpt_ring6/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(6));
        sat_TOT_VS_TECrings->SetBinContent(7,(float)cpt_ring7/cpt_cluster_strip_noclean + sat_TOT_VS_TECrings->GetBinContent(7));
      }

      cpt_norm_cluster++;
      cluster_TOT_VS_TECrings->SetBinContent(1,cpt_cluster_TEC1 + cluster_TOT_VS_TECrings->GetBinContent(1));
      cluster_TOT_VS_TECrings->SetBinContent(2,cpt_cluster_TEC2 + cluster_TOT_VS_TECrings->GetBinContent(2));
      cluster_TOT_VS_TECrings->SetBinContent(3,cpt_cluster_TEC3 + cluster_TOT_VS_TECrings->GetBinContent(3));
      cluster_TOT_VS_TECrings->SetBinContent(4,cpt_cluster_TEC4 + cluster_TOT_VS_TECrings->GetBinContent(4));
      cluster_TOT_VS_TECrings->SetBinContent(5,cpt_cluster_TEC5 + cluster_TOT_VS_TECrings->GetBinContent(5));
      cluster_TOT_VS_TECrings->SetBinContent(6,cpt_cluster_TEC6 + cluster_TOT_VS_TECrings->GetBinContent(6));
      cluster_TOT_VS_TECrings->SetBinContent(7,cpt_cluster_TEC7 + cluster_TOT_VS_TECrings->GetBinContent(7));
      */

    }
  }

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<"%)"<<endl;

  } //end EVENTS
  
  std::sort(NotDeadDetid.begin(),NotDeadDetid.end());
  int cpt_NotDeadDetid=0;
  for (int i=0;i<NotDeadDetid.size()-1;i++)
  {
    if (NotDeadDetid[i]!=NotDeadDetid[i-1]) cpt_NotDeadDetid++;
  }
  if (NotDeadDetid[NotDeadDetid.size()-1]!=NotDeadDetid[NotDeadDetid.size()-2]) cpt_NotDeadDetid++;

  cout<<"NUMBER OF STRIPS ALIVE AMONG THE 7122 CONCERNED IN FL AND FR CLUSTERS: "<<cpt_NotDeadDetid<<" / 7122"<<endl;


  sat_TOT_VS_TECrings->SetBinContent(1,(float)cpt_ring1/cpt_cluster_TEC1);
  sat_TOT_VS_TECrings->SetBinContent(2,(float)cpt_ring2/cpt_cluster_TEC2);
  sat_TOT_VS_TECrings->SetBinContent(3,(float)cpt_ring3/cpt_cluster_TEC3);
  sat_TOT_VS_TECrings->SetBinContent(4,(float)cpt_ring4/cpt_cluster_TEC4);
  sat_TOT_VS_TECrings->SetBinContent(5,(float)cpt_ring5/cpt_cluster_TEC5);
  sat_TOT_VS_TECrings->SetBinContent(6,(float)cpt_ring6/cpt_cluster_TEC6);
  sat_TOT_VS_TECrings->SetBinContent(7,(float)cpt_ring7/cpt_cluster_TEC7);

  int sum_cpt_cluster_TEC=cpt_cluster_TEC1+cpt_cluster_TEC2+cpt_cluster_TEC3+cpt_cluster_TEC4+cpt_cluster_TEC5+cpt_cluster_TEC6+cpt_cluster_TEC7;

  sat_TOT_VS_TECrings->SetBinContent(11,(float)cpt_ring1/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(12,(float)cpt_ring2/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(13,(float)cpt_ring3/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(14,(float)cpt_ring4/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(15,(float)cpt_ring5/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(16,(float)cpt_ring6/sum_cpt_cluster_TEC);
  sat_TOT_VS_TECrings->SetBinContent(17,(float)cpt_ring7/sum_cpt_cluster_TEC);


  // cout<<"NUMBER OF SATURATED CLUSTERS / TRACK"<<endl;

  // cout<<"NO CLEANING"<<endl;
  // cout<<"Saturation 254: "<<cpt_sat_254_noclean<<" / "<<cpt_cluster_strip_noclean<<endl;
  // cout<<"Saturation 255: "<<cpt_sat_255_noclean<<" / "<<cpt_cluster_strip_noclean<<endl;
  // cout<<"Saturation 254 & 255: "<<cpt_sat_TOT_noclean<<" / "<<cpt_cluster_strip_noclean<<endl;
  // cout<<endl;
  // cout<<"CLEAN 2 (CROSSTALK INVERSION)"<<endl;
  // cout<<"Saturation 254: "<<cpt_sat_254_clean2<<" / "<<cpt_cluster_strip_clean2<<endl;
  // cout<<"Saturation 255: "<<cpt_sat_255_clean2<<" / "<<cpt_cluster_strip_clean2<<endl;
  // cout<<"Saturation 254 & 255: "<<cpt_sat_TOT_clean2<<" / "<<cpt_cluster_strip_clean2<<endl;
  // cout<<endl;
  // cout<<"CLEAN 1"<<endl;
  // cout<<"Saturation 254: "<<cpt_sat_254_clean1<<" / "<<cpt_cluster_strip_clean1<<endl;
  // cout<<"Saturation 255: "<<cpt_sat_255_clean1<<" / "<<cpt_cluster_strip_clean1<<endl;
  // cout<<"Saturation 254 & 255: "<<cpt_sat_TOT_clean1<<" / "<<cpt_cluster_strip_clean1<<endl;
  // cout<<endl;
  // cout<<"NEW CLEAN"<<endl;
  // cout<<"Saturation 254: "<<cpt_sat_254_newclean<<" / "<<cpt_cluster_strip_newclean<<endl;
  // cout<<"Saturation 255: "<<cpt_sat_255_newclean<<" / "<<cpt_cluster_strip_newclean<<endl;
  // cout<<"Saturation 254 & 255: "<<cpt_sat_TOT_newclean<<" / "<<cpt_cluster_strip_newclean<<endl;

  // cout<<"NO CLEANING"<<endl;
  // cout<<"Saturation 254: "<<100*sat_254_noclean->GetMean()<<" %"<<endl;
  // cout<<"Saturation 255: "<<100*sat_255_noclean->GetMean()<<" %"<<endl;
  // cout<<"Saturation 254 & 255: "<<100*sat_TOT_noclean->GetMean()<<" %"<<endl;
  // cout<<endl;
  // cout<<"CLEAN 2 (CROSSTALK INVERSION)"<<endl;
  // cout<<"Saturation 254: "<<100*sat_254_clean2->GetMean()<<" %"<<endl;
  // cout<<"Saturation 255: "<<100*sat_255_clean2->GetMean()<<" %"<<endl;
  // cout<<"Saturation 254 & 255: "<<100*sat_TOT_clean2->GetMean()<<" %"<<endl;
  // cout<<endl;
  // cout<<"CLEAN 1"<<endl;
  // cout<<"Saturation 254: "<<100*sat_254_clean1->GetMean()<<" %"<<endl;
  // cout<<"Saturation 255: "<<100*sat_255_clean1->GetMean()<<" %"<<endl;
  // cout<<"Saturation 254 & 255: "<<100*sat_TOT_clean1->GetMean()<<" %"<<endl;
  // cout<<endl;
  // cout<<"NEW CLEAN"<<endl;
  // cout<<"Saturation 254: "<<100*sat_254_newclean->GetMean()<<" %"<<endl;
  // cout<<"Saturation 255: "<<100*sat_255_newclean->GetMean()<<" %"<<endl;
  // cout<<"Saturation 254 & 255: "<<100*sat_TOT_newclean->GetMean()<<" %"<<endl;

  SaveData->cd();

  //sat_254_noclean->Write();
  //sat_255_noclean->Write();
  sat_TOT_noclean->Write();
  /*sat_254_clean2->Write();
  sat_255_clean2->Write();
  sat_TOT_clean2->Write();
  sat_254_clean1->Write();
  sat_255_clean1->Write();
  sat_TOT_clean1->Write();
  sat_254_newclean->Write();
  sat_255_newclean->Write();
  sat_TOT_newclean->Write();*/

  sat_254_noclean_VS_p->Write();
  sat_255_noclean_VS_p->Write();
  sat_TOT_noclean_VS_p->Write();
  sat_TOT_noclean_VS_pt->Write();
  sat_TOT_noclean_VS_eta->Write();
  sat_TOT_noclean_VS_p_COUNT->Write();
  sat_TOT_noclean_VS_pt_COUNT->Write();
  sat_TOT_noclean_VS_eta_COUNT->Write();
  /*sat_254_clean2_VS_p->Write();
  sat_255_clean2_VS_p->Write();
  sat_TOT_clean2_VS_p->Write();
  sat_254_clean1_VS_p->Write();
  sat_255_clean1_VS_p->Write();
  sat_TOT_clean1_VS_p->Write();
  sat_254_newclean_VS_p->Write();
  sat_255_newclean_VS_p->Write();
  sat_TOT_newclean_VS_p->Write();*/

  check_pt_matched_muon->Write();
  check_eta_matched_muon->Write();
  check_phi_matched_muon->Write();

  pathlength_VS_eta->Write();
  p_VS_eta_sat->Write();
  pt_VS_eta_sat->Write();
  p_VS_eta_tot->Write();
  pt_VS_eta_tot->Write();

  clusterTrack_VS_eta->Write();

  Ih_VS_p_1->Write();
  Ih_VS_p_2->Write();
  sat_cluster_1_VS_p->Write();
  sat_cluster_1_VS_pt->Write();
  sat_cluster_2_VS_p->Write();
  sat_cluster_2_VS_pt->Write();

  pt_VS_eta->Write();
  p_VS_eta->Write();

  sat_TOT_VS_TECrings->Write();
  Cluster_layer_sat->Write();
  Cluster_layer_tot->Write();
  ClusterPerTrack_layer_tot->Scale(1./cpt_track);
  ClusterPerTrack_layer_tot->Write();

}
