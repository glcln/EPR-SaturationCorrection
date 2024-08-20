#define Wjets_checkCorrection_cxx
#include "Wjets_checkCorrection.h"
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
  int QII = accumulate(Q.begin(), Q.end(), 0);
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

bool clusterCleaning_changed(std::vector<int> ampls,  int crosstalkInv, uint8_t * exitCode)
{

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
                //cout<<NofMax<<" 1"<<endl;
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                //cout<<NofMax<<" 2"<<endl;
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
        //if(shapecdtn==0) cout<<"NofMax: "<<NofMax<<endl;
        return shapecdtn;
      }

//      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        Float_t coeff1=1.7;     Float_t coeff2=2.0;
        Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) {shapecdtn=true; if (C_D==255) {cout<<"CD=255, pas normal"<<endl;}} else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;if (C_D==255) {cout<<"CD=255, pas normal 2"<<endl;}} else if (exitCode) *exitCode=3;

                              // debug
                                 //if (shapecdtn==0) cout<<"MaxOnStart"<<endl;

                              }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                 // debug
                                  //if (shapecdtn==0) cout<<"MaxOnEnd"<<endl;
                                }
                }

                if(MaxInMiddle==true){

                        //cout<<"in ?"<<endl;
                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*coeffn*C_M+coeff2*coeffnn*C_D+2*noise || C_M>=254){
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
                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D>=254)
                                           && C_Mnn<=coeff1*coeffn*C_Mn+coeff2*coeffnn*C_M+2*noise
                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise) {
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
   //if(shapecdtn==0) cout<<"Return 2, NofMax: "<<NofMax<<endl;
   return shapecdtn;
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



void Wjets_checkCorrection::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  TFile* SaveData = new TFile("ROOT_SVG/Wjets_checkCorrection.root", "RECREATE");


  TH2F *Ih_VS_p_NoSat=new TH2F("Ih_VS_p_NoSat","Ih_VS_p_NoSat",300,0,15,300,0,30);
  TH2F *Ih_VS_p_SatNoCorr=new TH2F("Ih_VS_p_SatNoCorr","Ih_VS_p_SatNoCorr",300,0,15,300,0,30);
  TH2F *Ih_VS_p_SatCorr=new TH2F("Ih_VS_p_SatCorr","Ih_VS_p_SatCorr",300,0,15,300,0,30);
  TH2F *Ih_VS_p_OldCorr=new TH2F("Ih_VS_p_OldCorr","Ih_VS_p_OldCorr",300,0,15,300,0,30);

  int cpt_cluster_strip = 0;
  int cpt_sat_1=0, cpt_sat_2=0, cpt_sat_3=0, cpt_sat_4=0, cpt_sat_5=0, cpt_sat_more5=0;


  //nentries=1000000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    for (int itr=0; itr<ntracks; itr++)
    {
      float Ih = 0, Ih_nocorr = 0, Ih_oldcorr = 0;
      int cpt_isstrip = 0, cpt_isnotsat = 0;
      bool selection = true;
      bool selection2024 = true;
      vector <float> Ih_nocorr_charge_corr, Ih_charge_corr, Ih_oldcorr_charge_corr;
      vector <float> Ih_nocorr_pathlength, Ih_pathlength, Ih_oldcorr_pathlength;
      vector <int> Ih_nocorr_subdetId, Ih_subdetId, Ih_oldcorr_subdetId;
      vector <int> Ih_nocorr_moduleGeometry, Ih_moduleGeometry, Ih_oldcorr_moduleGeometry;
      vector <bool> Ih_nocorr_bool_cleaning, Ih_bool_cleaning, Ih_oldcorr_bool_cleaning;
      vector <bool> Ih_nocorr_mustBeInside, Ih_mustBeInside, Ih_oldcorr_mustBeInside;

      if (track_nvalidhits[itr]<8) selection=false;
      if (track_npixhits[itr]<2) selection=false;
      if (track_validfraction[itr]<0.8) selection=false;
      if (track_chi2[itr]>5) selection=false;
      if (abs(track_dxy[itr])>0.5) selection=false;
      if (abs(track_dz[itr])>0.5) selection=false;
      if (track_nvalidhits[itr]<=8) selection=false;
      if (!track_qual[itr]) selection=false;


      if (!track_qual[itr]) selection2024=false;
			if (track_pt[itr] < 0.5) selection2024 = false;
			if (track_nvalidhits[itr]<8) selection2024=false;
			if (track_npixhits[itr]<2) selection2024=false;
			if (track_validfraction[itr]<0.8) selection2024=false;
      if (abs(track_dxy[itr])>0.5) selection2024=false;
			if (abs(track_dz[itr])>0.5) selection2024=false;
      if (track_pterr[itr]/track_pt[itr]>0.25) selection2024=false;
			if (track_chi2[itr]>5) selection2024=false;
      for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
      {
        if (dedx_isstrip[idedx])
        {
          std::vector <int> ClusterCharge;
          for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++) ClusterCharge.push_back(strip_ampl[istrip]);
          if (!ClusterCharge.empty()) cpt_cluster_strip++;
        }
      }
      if (cpt_cluster_strip<10) selection2024=false;

      
      if (selection2024) //selection)
      {

        for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
        {
          if (dedx_isstrip[idedx])
          {
            cpt_isstrip++;

                // LAYER
            int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

                // SETUP
            std::vector <int> ClusterCharge;
            std::vector <int> ClusterStrip;
            int cpt_sat = 0;
            bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
            int Nleft = 0, Nright = 0;
            int Qcorr = 0;
            int Qoldcorr = 0;

            for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
            {
              ClusterCharge.push_back(strip_ampl[istrip]);
              ClusterStrip.push_back(strip_channel[istrip]);
              if (strip_ampl[istrip]>=254) cpt_sat++;
            }

                // Saturated clusters stats
            int cpt_sat_insidecluster = ConsecutiveSaturatedStrips(ClusterCharge);
            cpt_cluster_strip++;
            if (cpt_sat_insidecluster==1) cpt_sat_1++;
            if (cpt_sat_insidecluster==2) cpt_sat_2++;
            if (cpt_sat_insidecluster==3) cpt_sat_3++;
            if (cpt_sat_insidecluster==4) cpt_sat_4++;
            if (cpt_sat_insidecluster==5) cpt_sat_5++;
            if (cpt_sat_insidecluster>5) cpt_sat_more5++;

            
            if (!ClusterCharge.empty())
            {
              int value_max = *max_element(ClusterCharge.begin(), ClusterCharge.end());
              int i_max = getIndex(ClusterCharge,value_max);
              int sum_reco = accumulate(ClusterCharge.begin(), ClusterCharge.end(), 0);
              if (ClusterCharge.size()>=3) {Nleft = ClusterCharge[i_max-1]; Nright = ClusterCharge[i_max+1];}

                  // SHAPE
              if (ClusterCharge.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterCharge.size()-1 && Nleft>1.1*Nright) left=true;
              if (ClusterCharge.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterCharge.size()-1 && Nleft<0.9*Nright) right=true;
              if (ClusterCharge.size()>=3 && cpt_sat==1 && i_max>0 && i_max<ClusterCharge.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
              if (ClusterCharge.size()>=2 && cpt_sat==1 && i_max==0) FullLeft=true;
              if (ClusterCharge.size()>=2 && cpt_sat==1 && i_max==ClusterCharge.size()-1) FullRight=true;


                  // CORRECTION 
              Qoldcorr = Correction_OLD(ClusterCharge);
              if (center) Qcorr = Correction_wNeighbours(ClusterCharge, layer, "ROOT_SVG/Template_correction/Template_CENTER.txt", true);
              if (left || right) Qcorr = Correction_wNeighbours(ClusterCharge, layer, "ROOT_SVG/Template_correction/Template_LEFTRIGHT.txt", false);
              if (FullLeft || FullRight) Qcorr = Correction_FL_FR_xtalk(ClusterCharge, layer, "ROOT_SVG/Template_correction/Template_FLFR.txt");
              if (Qcorr < 254 || (!left && !right && !center && !FullLeft && !FullRight))
              {
                Qcorr = sum_reco;
                cpt_isnotsat++;
              }

              if (Qcorr >= 254)
              {
                for (int i=0; i<ClusterCharge.size(); i++)
                {
                  if (i != i_max) Qcorr += ClusterCharge[i];
                }
              }

              if (dedx_insideTkMod[idedx])
              {
                Ih_charge_corr.push_back(Qcorr);
                Ih_pathlength.push_back(dedx_pathlength[idedx]);
                Ih_subdetId.push_back(dedx_subdetid[idedx]);
                Ih_moduleGeometry.push_back(dedx_modulgeom[idedx]);
                Ih_mustBeInside.push_back(dedx_insideTkMod[idedx]);
                Ih_bool_cleaning.push_back(clusterCleaning_changed(ClusterCharge,0,0));

                Ih_nocorr_charge_corr.push_back(sum_reco);
                Ih_nocorr_pathlength.push_back(dedx_pathlength[idedx]);
                Ih_nocorr_subdetId.push_back(dedx_subdetid[idedx]);
                Ih_nocorr_moduleGeometry.push_back(dedx_modulgeom[idedx]);
                Ih_nocorr_mustBeInside.push_back(dedx_insideTkMod[idedx]);
                Ih_nocorr_bool_cleaning.push_back(clusterCleaning_changed(ClusterCharge,0,0));

                Ih_oldcorr_charge_corr.push_back(Qoldcorr);
                Ih_oldcorr_pathlength.push_back(dedx_pathlength[idedx]);
                Ih_oldcorr_subdetId.push_back(dedx_subdetid[idedx]);
                Ih_oldcorr_moduleGeometry.push_back(dedx_modulgeom[idedx]);
                Ih_oldcorr_mustBeInside.push_back(dedx_insideTkMod[idedx]);
                Ih_oldcorr_bool_cleaning.push_back(clusterCleaning_changed(ClusterCharge,0,0));
              }

            } // end check ClusterCharge.empty()
          } // end check dedx_isstrip
        } // end dedx loop

        int trash;
        Ih = getdEdX(Ih_charge_corr, Ih_pathlength, Ih_subdetId, Ih_moduleGeometry, Ih_bool_cleaning, Ih_mustBeInside, NULL,2, 0., 0, trash, trash);
        Ih_nocorr = getdEdX(Ih_nocorr_charge_corr, Ih_nocorr_pathlength, Ih_nocorr_subdetId, Ih_nocorr_moduleGeometry, Ih_nocorr_bool_cleaning, Ih_nocorr_mustBeInside, NULL,2, 0., 0, trash, trash);
        Ih_oldcorr = getdEdX(Ih_oldcorr_charge_corr, Ih_oldcorr_pathlength, Ih_oldcorr_subdetId, Ih_oldcorr_moduleGeometry, Ih_oldcorr_bool_cleaning, Ih_oldcorr_mustBeInside, NULL,2, 0., 0, trash, trash);

        if (cpt_isnotsat == cpt_isstrip) Ih_VS_p_NoSat->Fill(track_p[itr], Ih, track_prescale[itr]);
        else
        {
          Ih_VS_p_SatCorr->Fill(track_p[itr], Ih, track_prescale[itr]);
          Ih_VS_p_SatNoCorr->Fill(track_p[itr], Ih_nocorr, track_prescale[itr]);
          Ih_VS_p_OldCorr->Fill(track_p[itr], Ih_oldcorr, track_prescale[itr]);
        }

      } // end selection
    
    } // end track loop

    if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100*(float)jentry/(float)nentries<<" %)"<<endl;
  }

  SaveData->cd();

  Ih_VS_p_NoSat->Write();
  Ih_VS_p_SatNoCorr->Write();
  Ih_VS_p_SatCorr->Write();
  Ih_VS_p_OldCorr->Write();

  SaveData->Close();

  cout << " --------- Saturated clusters --------- " << endl;
  int cpt = cpt_sat_1 + cpt_sat_2 + cpt_sat_3 + cpt_sat_4 + cpt_sat_5 + cpt_sat_more5;
  cout << "Total cluster in strip: " << cpt_cluster_strip << endl;
  cout << "Total saturated clusters: " << cpt << endl;
  cout << "1 saturated strip: " << cpt_sat_1 << " / " << cpt << " (" << 100*(float)cpt_sat_1/(float)cpt << " %)" << endl;
  cout << "2 consecutives saturated strips: " << cpt_sat_2 << " / " << cpt << " (" << 100*(float)cpt_sat_2/(float)cpt << " %)" << endl;
  cout << "3 consecutives saturated strips: " << cpt_sat_3 << " / " << cpt << " (" << 100*(float)cpt_sat_3/(float)cpt << " %)" << endl;
  cout << "4 consecutives saturated strips: " << cpt_sat_4 << " / " << cpt << " (" << 100*(float)cpt_sat_4/(float)cpt << " %)" << endl;
  cout << "5 consecutives saturated strips: " << cpt_sat_5 << " / " << cpt << " (" << 100*(float)cpt_sat_5/(float)cpt << " %)" << endl;
  cout << "More than 5 consecutives saturated strips: " << cpt_sat_more5 << " / " << cpt << " (" << 100*(float)cpt_sat_more5/(float)cpt << " %)" << endl;
}