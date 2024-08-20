#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <vector>

using namespace TMath;
using namespace std;

void histo_Wjets_reconstructed()
{
    TFile *Sat_Wjets = new TFile("ROOT_SVG/Sat_Wjets_reconstructed.root", "READ");

    TH2F *sat_TOT_noclean_VS_p = (TH2F*)Sat_Wjets->Get("sat_TOT_noclean_VS_p");
    TH2F *sat_TOT_clean2_VS_p = (TH2F*)Sat_Wjets->Get("sat_TOT_clean2_VS_p");
    TH2F *sat_TOT_clean1_VS_p = (TH2F*)Sat_Wjets->Get("sat_TOT_clean1_VS_p");
    TH2F *sat_TOT_newclean_VS_p = (TH2F*)Sat_Wjets->Get("sat_TOT_newclean_VS_p");
    TProfile *PROF_sat_TOT_noclean_VS_p=sat_TOT_noclean_VS_p->ProfileX();
    TProfile *PROF_sat_TOT_clean2_VS_p=sat_TOT_clean2_VS_p->ProfileX();
    TProfile *PROF_sat_TOT_clean1_VS_p=sat_TOT_clean1_VS_p->ProfileX();
    TProfile *PROF_sat_TOT_newclean_VS_p=sat_TOT_newclean_VS_p->ProfileX();

    TCanvas *c_fractionVSp=new TCanvas("c_fractionVSp","c_fractionVSp",1100,700);
    c_fractionVSp->cd();
    gStyle->SetOptStat(0);

    PROF_sat_TOT_noclean_VS_p->SetLineColor(kBlack);
    PROF_sat_TOT_noclean_VS_p->SetLineWidth(2);

    PROF_sat_TOT_clean2_VS_p->SetLineColor(kBlue);
    PROF_sat_TOT_clean2_VS_p->SetLineWidth(2);
    
    PROF_sat_TOT_clean1_VS_p->SetLineColor(kRed);
    PROF_sat_TOT_clean1_VS_p->SetLineWidth(2);

    PROF_sat_TOT_newclean_VS_p->SetLineColor(kOrange);
    PROF_sat_TOT_newclean_VS_p->SetLineWidth(2);
    
    PROF_sat_TOT_noclean_VS_p->Draw();
    PROF_sat_TOT_clean2_VS_p->Draw("SAME");
    PROF_sat_TOT_clean1_VS_p->Draw("SAME");
    PROF_sat_TOT_newclean_VS_p->Draw("SAME");
    PROF_sat_TOT_newclean_VS_p->GetYaxis()->SetTitle("Nbr saturated clusters / tracks");
    PROF_sat_TOT_noclean_VS_p->GetXaxis()->SetTitle("p");

    auto l1= new TLegend(0.7,0.7,0.89,0.89);
    l1->AddEntry(PROF_sat_TOT_noclean_VS_p,"no clean","l");
    l1->AddEntry(PROF_sat_TOT_clean2_VS_p,"clean2","l");
    l1->AddEntry(PROF_sat_TOT_clean1_VS_p,"clean1","l");
    l1->AddEntry(PROF_sat_TOT_newclean_VS_p,"new clean","l");
    l1->SetBorderSize(0);
    l1->Draw();
}

void histo_Wjets_alltracks()
{
    TFile *Sat_Wjets_alltracks = new TFile("ROOT_SVG/Sat_Wjets_alltracks.root", "READ");
    
    TH2F *sat_TOT_noclean_VS_p = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_p");
    TH2F *sat_TOT_noclean_VS_p_COUNT = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_p_COUNT");
    TH2F *sat_TOT_noclean_VS_pt = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_pt");
    TH2F *sat_TOT_noclean_VS_pt_COUNT = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_pt_COUNT");
    TH2F *sat_TOT_noclean_VS_eta = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_eta");
    TH2F *sat_TOT_noclean_VS_eta_COUNT = (TH2F*)Sat_Wjets_alltracks->Get("sat_TOT_noclean_VS_eta_COUNT");
    TH2F *sat_clusterTrack_VS_eta = (TH2F*)Sat_Wjets_alltracks->Get("sat_clusterTrack_VS_eta");
    TH2F *pathlength_VS_eta = (TH2F*)Sat_Wjets_alltracks->Get("pathlength_VS_eta");
    TH2F *p_VS_eta_sat = (TH2F*)Sat_Wjets_alltracks->Get("p_VS_eta_sat");
    TH2F *pt_VS_eta_sat = (TH2F*)Sat_Wjets_alltracks->Get("pt_VS_eta_sat");
    TH2F *p_VS_eta_tot = (TH2F*)Sat_Wjets_alltracks->Get("p_VS_eta_tot");
    TH2F *pt_VS_eta_tot = (TH2F*)Sat_Wjets_alltracks->Get("pt_VS_eta_tot");
    TH1D *sat_TOT_VS_TECrings = (TH1D*)Sat_Wjets_alltracks->Get("sat_TOT_VS_TECrings");

    TH1D *Cluster_layer_sat = (TH1D*)Sat_Wjets_alltracks->Get("Cluster_layer_sat");
    TH1D *Cluster_layer_tot = (TH1D*)Sat_Wjets_alltracks->Get("Cluster_layer_tot");
    TH1D *ClusterPerTrack_layer_tot = (TH1D*)Sat_Wjets_alltracks->Get("ClusterPerTrack_layer_tot");

    TH2F *p_VS_eta = (TH2F*)Sat_Wjets_alltracks->Get("p_VS_eta");
    TH2F *pt_VS_eta = (TH2F*)Sat_Wjets_alltracks->Get("pt_VS_eta");
 
    
    
    sat_TOT_noclean_VS_p->Sumw2();
    sat_TOT_noclean_VS_p_COUNT->Sumw2();
    sat_TOT_noclean_VS_pt->Sumw2();
    sat_TOT_noclean_VS_pt_COUNT->Sumw2();
    sat_TOT_noclean_VS_eta->Sumw2();
    sat_TOT_noclean_VS_eta_COUNT->Sumw2();
    sat_clusterTrack_VS_eta->Sumw2();
    pathlength_VS_eta->Sumw2();
    p_VS_eta_sat->Sumw2();
    pt_VS_eta_sat->Sumw2();
    p_VS_eta_tot->Sumw2();
    pt_VS_eta_tot->Sumw2();
    sat_TOT_VS_TECrings->Sumw2();
    p_VS_eta->Sumw2();
    pt_VS_eta->Sumw2();
    Cluster_layer_sat->Sumw2();
    Cluster_layer_tot->Sumw2();
    ClusterPerTrack_layer_tot->Sumw2();

    TH1D *PROJECTION_p_VS_eta_sat=p_VS_eta_sat->ProjectionX();
    TH1D *PROJECTION_pt_VS_eta_sat=pt_VS_eta_sat->ProjectionX();
    TH1D *PROJECTION_p_VS_eta_tot=p_VS_eta_tot->ProjectionX();
    TH1D *PROJECTION_pt_VS_eta_tot=pt_VS_eta_tot->ProjectionX();

    p_VS_eta_sat->Divide(p_VS_eta_tot);
    pt_VS_eta_sat->Divide(pt_VS_eta_tot);
    Cluster_layer_sat->Divide(Cluster_layer_tot);
    PROJECTION_p_VS_eta_sat->Divide(PROJECTION_p_VS_eta_tot);
    PROJECTION_pt_VS_eta_sat->Divide(PROJECTION_pt_VS_eta_tot);
    

    TH1D *PROJECTION_alltracks_noclean_VS_p=sat_TOT_noclean_VS_p->ProjectionX();
    TH1D *PROJECTION_alltracks_noclean_VS_pt=sat_TOT_noclean_VS_pt->ProjectionX();
    TProfile *PROF_alltracks_noclean_VS_p=sat_TOT_noclean_VS_p->ProfileX();
    TProfile *PROF_alltracks_noclean_VS_p_COUNT=sat_TOT_noclean_VS_p_COUNT->ProfileX();
    TProfile *PROF_alltracks_noclean_VS_pt=sat_TOT_noclean_VS_pt->ProfileX();
    TProfile *PROF_alltracks_noclean_VS_pt_COUNT=sat_TOT_noclean_VS_pt_COUNT->ProfileX();
    TProfile *PROF_alltracks_noclean_VS_eta=sat_TOT_noclean_VS_eta->ProfileX();
    TProfile *PROF_alltracks_noclean_VS_eta_COUNT=sat_TOT_noclean_VS_eta_COUNT->ProfileX();

    TH2F *Ih_VS_p_1 = (TH2F*)Sat_Wjets_alltracks->Get("Ih_VS_p_1");
    TH2F *Ih_VS_p_2 = (TH2F*)Sat_Wjets_alltracks->Get("Ih_VS_p_2");
    TH2F *sat_cluster_1_VS_p = (TH2F*)Sat_Wjets_alltracks->Get("sat_cluster_1_VS_p");
    TH2F *sat_cluster_1_VS_pt = (TH2F*)Sat_Wjets_alltracks->Get("sat_cluster_1_VS_pt");
    TH2F *sat_cluster_2_VS_p = (TH2F*)Sat_Wjets_alltracks->Get("sat_cluster_2_VS_p");
    TH2F *sat_cluster_2_VS_pt = (TH2F*)Sat_Wjets_alltracks->Get("sat_cluster_2_VS_pt");
    sat_cluster_1_VS_p->Sumw2();
    sat_cluster_1_VS_pt->Sumw2();
    sat_cluster_2_VS_p->Sumw2();
    sat_cluster_2_VS_pt->Sumw2();
    TProfile *PROF_sat_cluster_1_VS_p=sat_cluster_1_VS_p->ProfileX();
    TProfile *PROF_sat_cluster_1_VS_pt=sat_cluster_1_VS_pt->ProfileX();
    TProfile *PROF_sat_cluster_2_VS_p=sat_cluster_2_VS_p->ProfileX();
    TProfile *PROF_sat_cluster_2_VS_pt=sat_cluster_2_VS_pt->ProfileX();
    TProfile *PROF_clusterTrack_VS_eta=sat_clusterTrack_VS_eta->ProfileX();
    PROF_sat_cluster_1_VS_p->SetLineColor(kBlue);
    PROF_sat_cluster_1_VS_p->SetLineWidth(2);
    PROF_sat_cluster_1_VS_p->GetXaxis()->SetTitle("p");
    PROF_sat_cluster_1_VS_p->GetYaxis()->SetTitle("Nbr of saturated clusters");
    PROF_sat_cluster_2_VS_p->SetLineColor(kRed);
    PROF_sat_cluster_2_VS_p->SetLineWidth(2);
    PROF_sat_cluster_1_VS_pt->SetLineColor(kBlue);
    PROF_sat_cluster_1_VS_pt->SetLineWidth(2);
    PROF_sat_cluster_1_VS_pt->GetXaxis()->SetTitle("p_{T}");
    PROF_sat_cluster_1_VS_pt->GetYaxis()->SetTitle("Nbr of saturated clusters");
    PROF_sat_cluster_2_VS_pt->SetLineColor(kRed);
    PROF_sat_cluster_2_VS_pt->SetLineWidth(2);
    sat_TOT_VS_TECrings->SetLineColor(kBlack);
    sat_TOT_VS_TECrings->SetLineWidth(2);
    sat_TOT_VS_TECrings->GetXaxis()->SetTitle("TEC ring");
    sat_TOT_VS_TECrings->GetYaxis()->SetTitle("Nbr of clusters / tracks");
    sat_TOT_VS_TECrings->SetTitle("No cleaning");
    Cluster_layer_sat->SetLineColor(kBlack);
    Cluster_layer_sat->SetLineWidth(2);
    Cluster_layer_sat->GetXaxis()->SetTitle("Tracker layer");
    Cluster_layer_sat->GetYaxis()->SetTitle("# saturated clusters in layer i / # clusters in layer i");
    Cluster_layer_sat->SetTitle("No cleaning");
    ClusterPerTrack_layer_tot->SetLineColor(kBlack);
    ClusterPerTrack_layer_tot->SetLineWidth(2);
    ClusterPerTrack_layer_tot->GetXaxis()->SetTitle("Tracker layer");
    ClusterPerTrack_layer_tot->GetYaxis()->SetTitle("# clusters in layer i / track");
    ClusterPerTrack_layer_tot->SetTitle("No cleaning");
    p_VS_eta->GetXaxis()->SetTitle("p");
    p_VS_eta->GetYaxis()->SetTitle("#eta");
    p_VS_eta->SetTitle("No cleaning");
    pt_VS_eta->GetXaxis()->SetTitle("p_{T}");
    pt_VS_eta->GetYaxis()->SetTitle("#eta");
    pt_VS_eta->SetTitle("No cleaning");
    PROJECTION_p_VS_eta_sat->GetXaxis()->SetTitle("p");
    PROJECTION_p_VS_eta_sat->GetYaxis()->SetTitle("Nbr of saturated clusters / track");
    PROJECTION_p_VS_eta_sat->SetTitle("No cleaning");
    PROJECTION_pt_VS_eta_sat->GetXaxis()->SetTitle("p_{T}");
    PROJECTION_pt_VS_eta_sat->GetYaxis()->SetTitle("Nbr of saturated clusters / track");
    PROJECTION_pt_VS_eta_sat->SetTitle("No cleaning");

    
    
    /*TCanvas *c_fractionVSp_alltracks_proj=new TCanvas("c_fractionVSp_alltracks_proj","c_fractionVSp_alltracks_proj",1100,700);
    c_fractionVSp_alltracks_proj->cd();
    gStyle->SetOptStat(0);

    PROJECTION_alltracks_noclean_VS_p->SetLineColor(kBlack);
    PROJECTION_alltracks_noclean_VS_p->SetLineWidth(2);

    PROJECTION_alltracks_noclean_VS_p->Draw();
    PROJECTION_alltracks_noclean_VS_p->GetYaxis()->SetTitle("Nbr saturated clusters / tracks * bin in p");
    PROJECTION_alltracks_noclean_VS_p->GetXaxis()->SetTitle("p");
    //PROJECTION_alltracks_noclean_VS_p->SetAxisRange(0,50,"X");
    PROJECTION_alltracks_noclean_VS_p->SetTitle("No cleaning - ProjectionX");
    c_fractionVSp_alltracks_proj->Update();

    
    TCanvas *c_fractionVSpt_alltracks_proj=new TCanvas("c_fractionVSpt_alltracks_proj","c_fractionVSpt_alltracks_proj",1100,700);
    c_fractionVSpt_alltracks_proj->cd();
    gStyle->SetOptStat(0);

    PROJECTION_alltracks_noclean_VS_pt->SetLineColor(kBlack);
    PROJECTION_alltracks_noclean_VS_pt->SetLineWidth(2);

    PROJECTION_alltracks_noclean_VS_pt->Draw();
    PROJECTION_alltracks_noclean_VS_pt->GetYaxis()->SetTitle("Nbr saturated clusters / tracks * bin en p_{T}");
    PROJECTION_alltracks_noclean_VS_pt->GetXaxis()->SetTitle("p_{T}");
    //PROJECTION_alltracks_noclean_VS_pt->SetAxisRange(0,50,"X");
    PROJECTION_alltracks_noclean_VS_pt->SetTitle("No cleaning - ProjectionX");
    c_fractionVSpt_alltracks_proj->Update();


    TCanvas *c_fractionVSp_alltracks=new TCanvas("c_fractionVSp_alltracks","c_fractionVSp_alltracks",1100,700);
    c_fractionVSp_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_p->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_p->SetLineWidth(2);

    PROF_alltracks_noclean_VS_p->Draw();
    PROF_alltracks_noclean_VS_p->GetYaxis()->SetTitle("Nbr saturated clusters / tracks");
    PROF_alltracks_noclean_VS_p->GetXaxis()->SetTitle("p");
    //PROF_alltracks_noclean_VS_p->SetAxisRange(0,50,"X");
    PROF_alltracks_noclean_VS_p->SetTitle("No cleaning");
    c_fractionVSp_alltracks->Update();


    TCanvas *c_fractionVSp_COUNT_alltracks=new TCanvas("c_fractionVSp_COUNT_alltracks","c_fractionVSp_COUNT_alltracks",1100,700);
    c_fractionVSp_COUNT_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_p_COUNT->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_p_COUNT->SetLineWidth(2);

    PROF_alltracks_noclean_VS_p_COUNT->Draw();
    PROF_alltracks_noclean_VS_p_COUNT->GetYaxis()->SetTitle("Nbr saturated clusters");
    PROF_alltracks_noclean_VS_p_COUNT->GetXaxis()->SetTitle("p");
    //PROF_alltracks_noclean_VS_p_COUNT->SetAxisRange(0,50,"X");
    PROF_alltracks_noclean_VS_p_COUNT->SetTitle("No cleaning");
    c_fractionVSp_COUNT_alltracks->Update();



    TCanvas *c_fractionVSpt_alltracks=new TCanvas("c_fractionVSpt_alltracks","c_fractionVSpt_alltracks",1100,700);
    c_fractionVSpt_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_pt->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_pt->SetLineWidth(2);

    PROF_alltracks_noclean_VS_pt->Draw();
    PROF_alltracks_noclean_VS_pt->GetYaxis()->SetTitle("Nbr saturated clusters / tracks");
    PROF_alltracks_noclean_VS_pt->GetXaxis()->SetTitle("p_{T}");
    //PROF_alltracks_noclean_VS_pt->SetAxisRange(0,50,"X");
    PROF_alltracks_noclean_VS_pt->SetTitle("No cleaning");
    c_fractionVSpt_alltracks->Update();



    TCanvas *c_fractionVSpt_COUNT_alltracks=new TCanvas("c_fractionVSpt_COUNT_alltracks","c_fractionVSpt_COUNT_alltracks",1100,700);
    c_fractionVSpt_COUNT_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_pt_COUNT->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_pt_COUNT->SetLineWidth(2);

    PROF_alltracks_noclean_VS_pt_COUNT->Draw();
    PROF_alltracks_noclean_VS_pt_COUNT->GetYaxis()->SetTitle("Nbr saturated clusters");
    PROF_alltracks_noclean_VS_pt_COUNT->GetXaxis()->SetTitle("p_{T}");
    //PROF_alltracks_noclean_VS_pt_COUNT->SetAxisRange(0,50,"X");
    PROF_alltracks_noclean_VS_pt_COUNT->SetTitle("No cleaning");
    c_fractionVSpt_COUNT_alltracks->Update();



    TCanvas *c_fractionVSeta_alltracks=new TCanvas("c_fractionVSeta_alltracks","c_fractionVSeta_alltracks",1100,700);
    c_fractionVSeta_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_eta->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_eta->SetLineWidth(2);

    PROF_alltracks_noclean_VS_eta->Draw();
    PROF_alltracks_noclean_VS_eta->GetYaxis()->SetTitle("Nbr saturated clusters / tracks");
    PROF_alltracks_noclean_VS_eta->GetXaxis()->SetTitle("#eta");
    PROF_alltracks_noclean_VS_eta->SetTitle("No cleaning");
    c_fractionVSeta_alltracks->Update();



    TCanvas *c_fractionVSeta_COUNT_alltracks=new TCanvas("c_fractionVSeta_COUNT_alltracks","c_fractionVSeta_COUNT_alltracks",1100,700);
    c_fractionVSeta_COUNT_alltracks->cd();
    gStyle->SetOptStat(0);

    PROF_alltracks_noclean_VS_eta_COUNT->SetLineColor(kBlack);
    PROF_alltracks_noclean_VS_eta_COUNT->SetLineWidth(2);

    PROF_alltracks_noclean_VS_eta_COUNT->Draw();
    PROF_alltracks_noclean_VS_eta_COUNT->GetYaxis()->SetTitle("Nbr saturated clusters");
    PROF_alltracks_noclean_VS_eta_COUNT->GetXaxis()->SetTitle("#eta");
    PROF_alltracks_noclean_VS_eta_COUNT->SetTitle("No cleaning");
    c_fractionVSeta_COUNT_alltracks->Update();


    TCanvas *c_clusterTrack_VS_eta=new TCanvas("c_clusterTrack_VS_eta","c_clusterTrack_VS_eta",1100,700);
    c_clusterTrack_VS_eta->cd();
    gStyle->SetOptStat(0);

    PROF_clusterTrack_VS_eta->SetLineColor(kBlack);
    PROF_clusterTrack_VS_eta->SetLineWidth(2);

    PROF_clusterTrack_VS_eta->Draw();
    PROF_clusterTrack_VS_eta->GetYaxis()->SetTitle("Nbr clusters");
    PROF_clusterTrack_VS_eta->GetXaxis()->SetTitle("#eta");
    PROF_clusterTrack_VS_eta->SetTitle("No cleaning");
    c_clusterTrack_VS_eta->Update();


    TCanvas *c_pathlength_VS_eta=new TCanvas("c_pathlength_VS_eta","c_pathlength_VS_eta",1100,700);
    c_pathlength_VS_eta->cd();
    gStyle->SetOptStat(0);
    pathlength_VS_eta->Draw("colz");
    pathlength_VS_eta->GetXaxis()->SetTitle("pathlength in TOB");
    pathlength_VS_eta->GetYaxis()->SetTitle("#eta");
    pathlength_VS_eta->SetTitle("No cleaning");
    c_pathlength_VS_eta->Update();
    
    
    
    TCanvas *c_Ih_VS_p_1=new TCanvas("c_Ih_VS_p_1","c_Ih_VS_p_1",1100,700);
    c_Ih_VS_p_1->cd();
    gStyle->SetOptStat(0);
    Ih_VS_p_1->Draw("colz");
    Ih_VS_p_1->GetXaxis()->SetTitle("p");
    Ih_VS_p_1->GetYaxis()->SetTitle("I_{h}");
    Ih_VS_p_1->SetTitle("No cleaning");


    TCanvas *c_Ih_VS_p_2=new TCanvas("c_Ih_VS_p_2","c_Ih_VS_p_2",1100,700);
    c_Ih_VS_p_2->cd();
    gStyle->SetOptStat(0);
    Ih_VS_p_2->Draw("colz");
    Ih_VS_p_2->GetXaxis()->SetTitle("p");
    Ih_VS_p_2->GetYaxis()->SetTitle("I_{h}");
    Ih_VS_p_2->SetTitle("No cleaning");
 
    
    TCanvas *c_p_VS_eta=new TCanvas("c_p_VS_eta","c_p_VS_eta",1100,700);
    c_p_VS_eta->cd();
    gStyle->SetOptStat(0);
    p_VS_eta_sat->Draw("colz");
    p_VS_eta_sat->GetXaxis()->SetTitle("p");
    p_VS_eta_sat->GetYaxis()->SetTitle("#eta");
    p_VS_eta_sat->SetTitle("No cleaning");
    c_p_VS_eta->Update();


    TCanvas *c_pt_VS_eta=new TCanvas("c_pt_VS_eta","c_pt_VS_eta",1100,700);
    c_pt_VS_eta->cd();
    gStyle->SetOptStat(0);
    pt_VS_eta_sat->Draw("colz");
    pt_VS_eta_sat->GetXaxis()->SetTitle("p_{T}");
    pt_VS_eta_sat->GetYaxis()->SetTitle("#eta");
    pt_VS_eta_sat->SetTitle("No cleaning");
    c_pt_VS_eta->Update();


    TCanvas *c_sat_VS_p=new TCanvas("c_sat_VS_p","c_sat_VS_p",1100,700);
    c_sat_VS_p->cd();
    gStyle->SetOptStat(0);
    PROF_sat_cluster_1_VS_p->Draw("hist");
    PROF_alltracks_noclean_VS_p_COUNT->Draw("hist SAME");
    PROF_sat_cluster_2_VS_p->Draw("hist SAME");
    PROF_sat_cluster_2_VS_p->GetYaxis()->SetTitle("Nbr saturated clusters");
    PROF_sat_cluster_2_VS_p->GetXaxis()->SetTitle("p");
    auto ld= new TLegend(0.6,0.35,0.85,0.55);
    ld->AddEntry(PROF_sat_cluster_1_VS_p,"Zone 1: I_{h} #geq -6p+10 && I_{h} #geq 4 && p #leq 2","l");
    ld->AddEntry(PROF_sat_cluster_2_VS_p,"Zone 2: p > 2 && I_{h} < 4","l");
    ld->AddEntry(PROF_alltracks_noclean_VS_p_COUNT,"all p","l");
    ld->SetBorderSize(0);
    ld->Draw();


    TCanvas *c_sat_VS_pt=new TCanvas("c_sat_VS_pt","c_sat_VS_pt",1100,700);
    c_sat_VS_pt->cd();
    gStyle->SetOptStat(0);
    PROF_sat_cluster_1_VS_pt->Draw("hist");
    PROF_alltracks_noclean_VS_pt_COUNT->Draw("hist SAME");
    PROF_sat_cluster_2_VS_pt->Draw("hist SAME");
    PROF_sat_cluster_2_VS_pt->GetYaxis()->SetTitle("Nbr saturated clusters");
    PROF_sat_cluster_2_VS_pt->GetXaxis()->SetTitle("p_{T}");
    auto le= new TLegend(0.6,0.35,0.85,0.55);
    le->AddEntry(PROF_sat_cluster_1_VS_p,"Zone 1: I_{h} #geq -6p+10 && I_{h} #geq 4 && p #leq 2","l");
    le->AddEntry(PROF_sat_cluster_2_VS_p,"Zone 2: p > 2 && I_{h} < 4","l");
    le->AddEntry(PROF_alltracks_noclean_VS_pt_COUNT,"all p_{T}","l");
    le->SetBorderSize(0);
    le->Draw();*/

    TCanvas *c_sat_TECring=new TCanvas("c_sat_TECring","c_sat_TECring",1100,700);
    c_sat_TECring->cd();
    gStyle->SetOptStat(0);
    sat_TOT_VS_TECrings->Draw("hist");


    TCanvas *c_sat_layer=new TCanvas("c_sat_layer","c_sat_layer",1100,700);
    c_sat_layer->cd();
    gStyle->SetOptStat(0);
    Cluster_layer_sat->Draw("E0");

    TCanvas *c_cluster_layer=new TCanvas("c_cluster_layer","c_cluster_layer",1100,700);
    c_cluster_layer->cd();
    gStyle->SetOptStat(0);
    ClusterPerTrack_layer_tot->Draw("E0");

    TCanvas *c_pVSeta=new TCanvas("c_pVSeta","c_pVSeta",1100,700);
    c_pVSeta->cd();
    gStyle->SetOptStat(0);
    p_VS_eta->Draw("colz");

    TCanvas *c_ptVSeta=new TCanvas("c_ptVSeta","c_ptVSeta",1100,700);
    c_ptVSeta->cd();
    gStyle->SetOptStat(0);
    pt_VS_eta->Draw("colz");

    TCanvas *c_Proj_p=new TCanvas("c_Proj_p","c_Proj_p",1100,700);
    c_Proj_p->cd();
    gStyle->SetOptStat(0);
    PROJECTION_p_VS_eta_sat->Draw("E0");

    TCanvas *c_Proj_pt=new TCanvas("c_Proj_pt","c_Proj_pt",1100,700);
    c_Proj_pt->cd();
    gStyle->SetOptStat(0);
    PROJECTION_pt_VS_eta_sat->Draw("E0");

}

void histo_MuonPU_1strip()
{
    TFile *Sat_MuonPU_alltracks = new TFile("ROOT_SVG/Sat_MuonPU_alltracks.root", "READ");

    TH1F *Sigdigi_left = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_left");
    TH1F *Sigdigi_right = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_right");
    TH1F *Sigdigi_center = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_center");
    TH1F *Sigdigi_FullLeft = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft");
    TH1F *Sigdigi_FullRight = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight");
    TH1F *Sigdigi_FullLeft_recons = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft_recons");
    TH1F *Sigdigi_FullRight_recons = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight_recons");
    TH1F *Sigdigi_FullLeft_noBorder = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft_noBorder");
    TH1F *Sigdigi_FullRight_noBorder = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight_noBorder");
    TH1F *Sigdigi_FullLeft_Border = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft_Border");
    TH1F *Sigdigi_FullRight_Border = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight_Border");
    TH1F *Sigdigi_FullLeft_layer = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft_layer");
    TH1F *Sigdigi_FullRight_layer = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight_layer");
    TH1F *Sigdigi_layer_maxsat = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_layer_maxsat");
    
    TH1F *cpt_sat = (TH1F*)Sat_MuonPU_alltracks->Get("cpt_sat");
    TH1F *Compare_eloss_over_sigdigi = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_eloss_over_sigdigi");
    TH1F *Compare_eloss_over_sigdigi_nosat = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_eloss_over_sigdigi_nosat");

    TH1F *Compare_QII_Over_SigDigi = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi");
    TH1F *Compare_QII_Over_SigDigi_L = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_L");
    TH1F *Compare_QII_Over_SigDigi_R = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_R");
    TH1F *Compare_QII_Over_SigDigi_C = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_C");
    TH1F *Compare_QII_Over_SigDigi_FL = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_FL");
    TH1F *Compare_QII_Over_SigDigi_FR = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_FR");

    TH1F *Compare_QII_Minus_SigDigi = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi");
    TH1F *Compare_QII_Minus_SigDigi_L = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_L");
    TH1F *Compare_QII_Minus_SigDigi_R = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_R");
    TH1F *Compare_QII_Minus_SigDigi_C = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_C");
    TH1F *Compare_QII_Minus_SigDigi_FL = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_FL");
    TH1F *Compare_QII_Minus_SigDigi_FR = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_FR");

    TH1F *check_L_Original = (TH1F*)Sat_MuonPU_alltracks->Get("check_L_Original");
    TH1F *check_R_Original = (TH1F*)Sat_MuonPU_alltracks->Get("check_R_Original");
    TH1F *check_C_Original = (TH1F*)Sat_MuonPU_alltracks->Get("check_C_Original");
    TH1F *check_FL_Original = (TH1F*)Sat_MuonPU_alltracks->Get("check_FL_Original");
    TH1F *check_FR_Original = (TH1F*)Sat_MuonPU_alltracks->Get("check_FR_Original");

    TH2F *Ih_VS_p_1 = (TH2F*)Sat_MuonPU_alltracks->Get("Ih_VS_p_1");
    TH2F *Ih_VS_p_2 = (TH2F*)Sat_MuonPU_alltracks->Get("Ih_VS_p_2");

    TH2F *FirstNeighbour_VS_Barycenter = (TH2F*)Sat_MuonPU_alltracks->Get("FirstNeighbour_VS_Barycenter");
    TH2F *MaxNeighbourRatio_VS_Barycenter = (TH2F*)Sat_MuonPU_alltracks->Get("MaxNeighbourRatio_VS_Barycenter");
    TH2F *MaxOverMaxNeighbour_VS_MaxNeihbourRatio = (TH2F*)Sat_MuonPU_alltracks->Get("MaxOverMaxNeighbour_VS_MaxNeihbourRatio");
    TH2F *TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio = (TH2F*)Sat_MuonPU_alltracks->Get("TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio");
    TH2F *TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio = (TH2F*)Sat_MuonPU_alltracks->Get("TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio");
    TH2F *TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio = (TH2F*)Sat_MuonPU_alltracks->Get("TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio");
    TH2F *TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio = (TH2F*)Sat_MuonPU_alltracks->Get("TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio");
    TH2F *MaxOverMaxNeighbour_VS_Barycenter = (TH2F*)Sat_MuonPU_alltracks->Get("MaxOverMaxNeighbour_VS_Barycenter");

    TH1F *StripSat_Over5_layer = (TH1F*)Sat_MuonPU_alltracks->Get("StripSat_Over5_layer");
    TH2F *StripSat_Over5_layer_eta = (TH2F*)Sat_MuonPU_alltracks->Get("StripSat_Over5_layer_eta");
    TH1F *ClusterUnderCorrection_layer_B = (TH1F*)Sat_MuonPU_alltracks->Get("ClusterUnderCorrection_layer_B");
    TH1F *ClusterUnderCorrection_layer_BandN = (TH1F*)Sat_MuonPU_alltracks->Get("ClusterUnderCorrection_layer_BandN");

    TH1F *ClusterUnderCorrection_layer = (TH1F*)Sat_MuonPU_alltracks->Get("ClusterUnderCorrection_layer");
    TH1F *ClusterUnderCorrection_layer_all = (TH1F*)Sat_MuonPU_alltracks->Get("ClusterUnderCorrection_layer_all");
    TH2F *ClusterUnderCorrection_sizeVSeta = (TH2F*)Sat_MuonPU_alltracks->Get("ClusterUnderCorrection_sizeVSeta");
    TH2F *ClusterGoodCorrection_sizeVSeta = (TH2F*)Sat_MuonPU_alltracks->Get("ClusterGoodCorrection_sizeVSeta");

    TH2F *MaxcorrOverMax_VS_eta_TIBL2 = (TH2F*)Sat_MuonPU_alltracks->Get("MaxcorrOverMax_VS_eta_TIBL2");
    TH2F *MaxcorrOverMax_VS_shape_TIDR2 = (TH2F*)Sat_MuonPU_alltracks->Get("MaxcorrOverMax_VS_shape_TIDR2");
    ClusterUnderCorrection_layer->Sumw2();
    ClusterUnderCorrection_layer_all->Sumw2();
    MaxcorrOverMax_VS_eta_TIBL2->Sumw2();
    MaxcorrOverMax_VS_shape_TIDR2->Sumw2();

    std::vector <TH2D*> histoVector_IonizationNeighbour;
    std::vector <TH2D*> histoVector_MaxNeighbourBarycenter;
    std::vector <TH2D*> histoVector_MaxVmaxBarycenter;
    for (int i=0; i<=29; i++)
    {
        TString histName = Form("MaxOverMaxNeighbour_VS_MaxNeihbourRatio_layer_%d", i);
        TH2D *histo = (TH2D*)Sat_MuonPU_alltracks->Get(histName);
        histoVector_IonizationNeighbour.push_back(histo);
    }
    for (int i=0; i<=29; i++)
    {
        TString histName = Form("MaxNeighbourRatio_VS_Barycenter_layer_%d", i);
        TH2D *histo = (TH2D*)Sat_MuonPU_alltracks->Get(histName);
        histoVector_MaxNeighbourBarycenter.push_back(histo);
    }
    for (int i=0; i<=20; i++)
    {
        TString histName = Form("MaxOverMaxNeighbour_VS_Barycenter_layer_%d", i);
        TH2D *histo = (TH2D*)Sat_MuonPU_alltracks->Get(histName);
        histoVector_MaxVmaxBarycenter.push_back(histo);
    }
    


    vector <TLine*> line;
    for (int i=1;i<=9;i++)
    {
        TLine *l1 = new TLine(0, 255*cpt_sat->GetBinContent(i), 50, 255*cpt_sat->GetBinContent(i));
        l1->SetLineColor(kRed);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        line.push_back(l1);
    }

    TLine *l2 = new TLine(-3, 1.1, 3, 1.1);
    TLine *l4= new TLine(-3, 0.909, 3, 0.909);
    l2->SetLineColor(kRed);
    l2->SetLineStyle(2);
    l2->SetLineWidth(2);
    l4->SetLineColor(kRed);
    l4->SetLineStyle(2);
    l4->SetLineWidth(2);
    line.push_back(l2);
    line.push_back(l4);

    auto l3= new TLegend(0.7,0.7,0.83,0.83);
    l3->AddEntry(line[0],"255 ADC","l");
    l3->SetBorderSize(0);

    Sigdigi_left->SetLineColor(kBlack);
    Sigdigi_left->SetLineWidth(2);
    Sigdigi_left->SetTitle("1st left neighbor/1st right neighbor > 10\%");
    Sigdigi_right->SetLineColor(kBlack);
    Sigdigi_right->SetLineWidth(2);
    Sigdigi_right->SetTitle("1st right neighbor/1st left neighbor > 10\%");
    Sigdigi_center->SetLineColor(kBlack);
    Sigdigi_center->SetLineWidth(2);
    Sigdigi_center->SetTitle("|1- 1st left neighbor/1st right neighbor|<10%");
    Sigdigi_FullLeft->SetLineColor(kBlack);
    Sigdigi_FullLeft->SetLineWidth(2);
    Sigdigi_FullLeft->SetTitle("Max on left");
    Sigdigi_FullRight->SetLineColor(kBlack);
    Sigdigi_FullRight->SetLineWidth(2);
    Sigdigi_FullRight->SetTitle("Max on right");
    Sigdigi_FullLeft_recons->SetLineColor(kGreen+1);
    Sigdigi_FullLeft_recons->SetLineWidth(2);
    Sigdigi_FullLeft_recons->SetTitle("Max on left");
    Sigdigi_FullRight_recons->SetLineColor(kGreen+1);
    Sigdigi_FullRight_recons->SetLineWidth(2);
    Sigdigi_FullRight_recons->SetTitle("Max on right");
    Compare_eloss_over_sigdigi->SetLineColor(kBlack);
    Compare_eloss_over_sigdigi->SetLineWidth(2);
    Compare_eloss_over_sigdigi->SetTitle("Qloss/Sum_sigdigi");
    Compare_eloss_over_sigdigi_nosat->SetLineColor(kRed);
    Compare_eloss_over_sigdigi_nosat->SetLineWidth(2);
    Compare_eloss_over_sigdigi_nosat->SetTitle("Qloss/Sum_sigdigi");
    Sigdigi_FullLeft_noBorder->SetLineColor(kBlack);
    Sigdigi_FullLeft_noBorder->SetLineWidth(2);
    Sigdigi_FullLeft_noBorder->SetTitle("Max on left - not on border");
    Sigdigi_FullRight_noBorder->SetLineColor(kBlack);
    Sigdigi_FullRight_noBorder->SetLineWidth(2);
    Sigdigi_FullRight_noBorder->SetTitle("Max on right - not on border");
    Sigdigi_FullLeft_Border->SetLineColor(kBlack);
    Sigdigi_FullLeft_Border->SetLineWidth(2);
    Sigdigi_FullLeft_Border->SetTitle("Max on left - on border");
    Sigdigi_FullRight_Border->SetLineColor(kBlack);
    Sigdigi_FullRight_Border->SetLineWidth(2);
    Sigdigi_FullRight_Border->SetTitle("Max on right - on border");
    Sigdigi_FullLeft_layer->SetLineColor(kBlue);
    Sigdigi_FullLeft_layer->SetLineWidth(2);
    Sigdigi_FullLeft_layer->Scale(1./Sigdigi_FullLeft_layer->Integral());
    Sigdigi_FullRight_layer->SetLineColor(kRed);
    Sigdigi_FullRight_layer->SetLineWidth(2);
    Sigdigi_FullRight_layer->Scale(1./Sigdigi_FullRight_layer->Integral());
    Sigdigi_layer_maxsat->SetLineColor(kBlack);
    Sigdigi_layer_maxsat->SetLineWidth(2);
    Sigdigi_layer_maxsat->Scale(1./Sigdigi_layer_maxsat->Integral());
    Sigdigi_layer_maxsat->SetTitle("Layer of max saturated");
    Compare_QII_Over_SigDigi->SetLineColor(kBlack);
    Compare_QII_Over_SigDigi->SetLineWidth(2);
    Compare_QII_Over_SigDigi->SetTitle("QII/Sum_sigdigi");
    Compare_QII_Over_SigDigi_L->SetLineColor(kBlack);
    Compare_QII_Over_SigDigi_L->SetLineWidth(2);
    Compare_QII_Over_SigDigi_L->Scale(1./Compare_QII_Over_SigDigi_L->Integral());
    Compare_QII_Over_SigDigi_R->SetLineColor(kRed);
    Compare_QII_Over_SigDigi_R->SetLineWidth(2);
    Compare_QII_Over_SigDigi_R->Scale(1./Compare_QII_Over_SigDigi_R->Integral());
    Compare_QII_Over_SigDigi_C->SetLineColor(kOrange);
    Compare_QII_Over_SigDigi_C->SetLineWidth(2);
    Compare_QII_Over_SigDigi_C->Scale(1./Compare_QII_Over_SigDigi_C->Integral());
    Compare_QII_Over_SigDigi_FL->SetLineColor(kGreen);
    Compare_QII_Over_SigDigi_FL->SetLineWidth(2);
    Compare_QII_Over_SigDigi_FL->Scale(1./Compare_QII_Over_SigDigi_FL->Integral());
    Compare_QII_Over_SigDigi_FR->SetLineColor(kBlue);
    Compare_QII_Over_SigDigi_FR->SetLineWidth(2);
    Compare_QII_Over_SigDigi_FR->Scale(1./Compare_QII_Over_SigDigi_FR->Integral());
    
    Compare_QII_Minus_SigDigi->SetLineColor(kBlack);
    Compare_QII_Minus_SigDigi->SetLineWidth(2);
    Compare_QII_Minus_SigDigi->SetTitle("QII - Sum_sigdigi");
    Compare_QII_Minus_SigDigi_L->SetLineColor(kBlack);
    Compare_QII_Minus_SigDigi_L->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_L->Scale(1./Compare_QII_Minus_SigDigi_L->Integral());
    Compare_QII_Minus_SigDigi_R->SetLineColor(kRed);
    Compare_QII_Minus_SigDigi_R->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_R->Scale(1./Compare_QII_Minus_SigDigi_R->Integral());
    Compare_QII_Minus_SigDigi_C->SetLineColor(kOrange);
    Compare_QII_Minus_SigDigi_C->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_C->Scale(1./Compare_QII_Minus_SigDigi_C->Integral());
    Compare_QII_Minus_SigDigi_FL->SetLineColor(kGreen);
    Compare_QII_Minus_SigDigi_FL->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_FL->Scale(1./Compare_QII_Minus_SigDigi_FL->Integral());
    Compare_QII_Minus_SigDigi_FR->SetLineColor(kBlue);
    Compare_QII_Minus_SigDigi_FR->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_FR->Scale(1./Compare_QII_Minus_SigDigi_FR->Integral());

    check_L_Original->SetLineColor(kBlack);
    check_L_Original->SetLineWidth(2);
    check_L_Original->Scale(1./check_L_Original->Integral());
    check_R_Original->SetLineColor(kRed);
    check_R_Original->SetLineWidth(2);
    check_R_Original->Scale(1./check_R_Original->Integral());
    check_C_Original->SetLineColor(kOrange);
    check_C_Original->SetLineWidth(2);
    check_C_Original->Scale(1./check_C_Original->Integral());
    check_FL_Original->SetLineColor(kGreen);
    check_FL_Original->SetLineWidth(2);
    check_FL_Original->Scale(1./check_FL_Original->Integral());
    check_FR_Original->SetLineColor(kBlue);
    check_FR_Original->SetLineWidth(2);
    check_FR_Original->Scale(1./check_FR_Original->Integral());

    FirstNeighbour_VS_Barycenter->GetXaxis()->SetTitle("i_{max} - barycenter");
    FirstNeighbour_VS_Barycenter->GetYaxis()->SetTitle("1^{st} Left neighbour / 1^{st} Right neighbour");
    FirstNeighbour_VS_Barycenter->SetTitle("No cleaning");
    MaxNeighbourRatio_VS_Barycenter->GetXaxis()->SetTitle("i_{max} - barycenter");
    MaxNeighbourRatio_VS_Barycenter->GetYaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    MaxNeighbourRatio_VS_Barycenter->SetTitle("No cleaning");
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->SetTitle("No cleaning");
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->SetTitle("No cleaning - TIB");
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->SetTitle("No cleaning - TOB");
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->SetTitle("No cleaning - TID");
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->SetTitle("No cleaning - TEC");
    MaxOverMaxNeighbour_VS_Barycenter->GetXaxis()->SetTitle("i_{max} - barycenter");
    MaxOverMaxNeighbour_VS_Barycenter->GetYaxis()->SetTitle("Max / V_{max}");
    MaxOverMaxNeighbour_VS_Barycenter->SetTitle("No cleaning");
    StripSat_Over5_layer->GetXaxis()->SetTitle("layer");
    StripSat_Over5_layer->GetYaxis()->SetTitle("# Clusters with more than 5 sat strips");
    StripSat_Over5_layer->SetTitle("No cleaning");
    StripSat_Over5_layer_eta->GetXaxis()->SetTitle("layer");
    StripSat_Over5_layer_eta->GetYaxis()->SetTitle("#eta");
    StripSat_Over5_layer_eta->SetTitle("No cleaning");
    ClusterUnderCorrection_layer->GetXaxis()->SetTitle("layer");
    ClusterUnderCorrection_layer->GetYaxis()->SetTitle("# Clusters undercorrected");
    ClusterUnderCorrection_layer->SetTitle("No cleaning - Neighbours correction");
    ClusterUnderCorrection_layer->SetLineColor(kBlack);
    ClusterUnderCorrection_layer->SetLineWidth(2);
    ClusterUnderCorrection_layer_B->GetXaxis()->SetTitle("layer");
    ClusterUnderCorrection_layer_B->GetYaxis()->SetTitle("# Clusters undercorrected");
    ClusterUnderCorrection_layer_B->SetTitle("No cleaning - Barycenter correction");
    ClusterUnderCorrection_layer_B->SetLineColor(kRed);
    ClusterUnderCorrection_layer_B->SetLineWidth(2);
    ClusterUnderCorrection_layer_BandN->GetXaxis()->SetTitle("layer");
    ClusterUnderCorrection_layer_BandN->GetYaxis()->SetTitle("# Clusters undercorrected");
    ClusterUnderCorrection_layer_BandN->SetTitle("No cleaning - Barycenter AND Neighbours correction");
    ClusterUnderCorrection_layer_BandN->SetLineColor(kBlue);
    ClusterUnderCorrection_layer_BandN->SetLineWidth(2);
    ClusterUnderCorrection_sizeVSeta->GetXaxis()->SetTitle("Cluster size");
    ClusterUnderCorrection_sizeVSeta->GetYaxis()->SetTitle("#eta");
    ClusterGoodCorrection_sizeVSeta->SetTitle("No cleaning - cluster well corrected");
    ClusterGoodCorrection_sizeVSeta->GetXaxis()->SetTitle("Cluster size");
    ClusterGoodCorrection_sizeVSeta->GetYaxis()->SetTitle("#eta");
    ClusterUnderCorrection_sizeVSeta->SetTitle("No cleaning - cluster miscorrected");
    ClusterUnderCorrection_layer->GetXaxis()->SetTitle("layer");
    ClusterUnderCorrection_layer->GetYaxis()->SetTitle("f#frac{# miscorrected clusters in the layer}{All corrected clusters in the layer}");
    ClusterUnderCorrection_layer->SetTitle("No cleaning - Neighbours correction");
    ClusterUnderCorrection_layer->SetLineColor(kBlack);
    ClusterUnderCorrection_layer->SetLineWidth(2);
    ClusterUnderCorrection_layer_all->GetXaxis()->SetTitle("layer");
    ClusterUnderCorrection_layer_all->GetYaxis()->SetTitle("# Clusters undercorrected");
    ClusterUnderCorrection_layer_all->SetTitle("No cleaning");
    //MaxcorrOverMax_VS_shape_TIBL2->GetXaxis()->SetTitle("");
    //MaxcorrOverMax_VS_shape_TIBL2->GetYaxis()->SetTitle("Max_{corr}/Max");
    //MaxcorrOverMax_VS_shape_TIBL2->SetTitle("No cleaning - TIB L2");
    MaxcorrOverMax_VS_shape_TIDR2->GetXaxis()->SetTitle("");
    MaxcorrOverMax_VS_shape_TIDR2->GetYaxis()->SetTitle("Max_{corr}/Max");
    MaxcorrOverMax_VS_shape_TIDR2->SetTitle("No cleaning - TID R2");

/*  
    TCanvas *c_left=new TCanvas("c_left","c_left",1100,700);
    c_left->cd();
    gStyle->SetOptStat(0);
    Sigdigi_left->Draw("hist text");
    line[0]->Draw("same");
    l3->Draw();

    TCanvas *c_right=new TCanvas("c_right","c_right",1100,700);
    c_right->cd();
    gStyle->SetOptStat(0);
    Sigdigi_right->Draw("hist text");
    line[1]->Draw("same");
    l3->Draw();

    TCanvas *c_center=new TCanvas("c_center","c_center",1100,700);
    c_center->cd();
    gStyle->SetOptStat(0);
    Sigdigi_center->Draw("hist text");
    line[2]->Draw("same");
    l3->Draw();

    TCanvas *c_FullLeft=new TCanvas("c_FullLeft","c_FullLeft",1100,700);
    c_FullLeft->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullLeft->Draw("hist text");
    Sigdigi_FullLeft_recons->Draw("hist same");
    line[3]->Draw("same");
    l3->Draw();

    TCanvas *c_FullRight=new TCanvas("c_FullRight","c_FullRight",1100,700);
    c_FullRight->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullRight->Draw("hist text");
    Sigdigi_FullRight_recons->Draw("hist same");
    line[4]->Draw("same");
    l3->Draw();


    TCanvas *c_NeighbourBarycenter=new TCanvas("c_NeighbourBarycenter","c_NeighbourBarycenter",1100,700);
    c_NeighbourBarycenter->cd();
    gStyle->SetOptStat(0);
    FirstNeighbour_VS_Barycenter->GetZaxis()->SetRangeUser(0,250);
    FirstNeighbour_VS_Barycenter->Draw("colz");
    line[9]->Draw("same");
    line[10]->Draw("same");

    
    TCanvas *c_MaxNeighbourBarycenter=new TCanvas("c_MaxNeighbourBarycenter","c_MaxNeighbourBarycenter",1100,700);
    c_MaxNeighbourBarycenter->cd();
    gStyle->SetOptStat(0);
    MaxNeighbourRatio_VS_Barycenter->GetZaxis()->SetRangeUser(0,250);
    MaxNeighbourRatio_VS_Barycenter->GetYaxis()->SetRangeUser(0.7,10);
    MaxNeighbourRatio_VS_Barycenter->Draw("colz");


    TCanvas *c_MaxOverNeighbour_VS_NeighbourRatio=new TCanvas("c_MaxOverNeighbour_VS_NeighbourRatio","c_MaxOverNeighbour_VS_NeighbourRatio",1100,700);
    c_MaxOverNeighbour_VS_NeighbourRatio->cd();
    gStyle->SetOptStat(0);
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetZaxis()->SetRangeUser(0,250);
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetRangeUser(0.7,20);
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetRangeUser(0.7,10);
    MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Draw("colz");

    TCanvas *c_TIB_MaxOverNeighbour_VS_NeighbourRatio=new TCanvas("c_TIB_MaxOverNeighbour_VS_NeighbourRatio","c_TIB_MaxOverNeighbour_VS_NeighbourRatio",1100,700);
    c_TIB_MaxOverNeighbour_VS_NeighbourRatio->cd();
    gStyle->SetOptStat(0);
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetZaxis()->SetRangeUser(0,250);
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetRangeUser(0.7,20);
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetRangeUser(0.7,10);
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Draw("colz");
    TF1 *fit_TIB = new TF1("fit_TIB","[0]+[1]/x",0,10);
    TIB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fit("fit_TIB");

    TCanvas *c_TOB_MaxOverNeighbour_VS_NeighbourRatio=new TCanvas("c_TOB_MaxOverNeighbour_VS_NeighbourRatio","c_TOB_MaxOverNeighbour_VS_NeighbourRatio",1100,700);
    c_TOB_MaxOverNeighbour_VS_NeighbourRatio->cd();
    gStyle->SetOptStat(0);
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetZaxis()->SetRangeUser(0,250);
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetRangeUser(0.7,20);
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetRangeUser(0.7,10);
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Draw("colz");
    TF1 *fit_TOB = new TF1("fit_TOB","[0]+[1]/x",0,10);
    TOB_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fit("fit_TOB");
    
    TCanvas *c_TID_MaxOverNeighbour_VS_NeighbourRatio=new TCanvas("c_TID_MaxOverNeighbour_VS_NeighbourRatio","c_TID_MaxOverNeighbour_VS_NeighbourRatio",1100,700);
    c_TID_MaxOverNeighbour_VS_NeighbourRatio->cd();
    gStyle->SetOptStat(0);
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetZaxis()->SetRangeUser(0,250);
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetRangeUser(0.7,20);
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetRangeUser(0.7,10);
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Draw("colz");
    TF1 *fit_TID = new TF1("fit_TID","[0]+[1]/x",0,10);
    TID_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fit("fit_TID");
    
    TCanvas *c_TEC_MaxOverNeighbour_VS_NeighbourRatio=new TCanvas("c_TEC_MaxOverNeighbour_VS_NeighbourRatio","c_TEC_MaxOverNeighbour_VS_NeighbourRatio",1100,700);
    c_TEC_MaxOverNeighbour_VS_NeighbourRatio->cd();
    gStyle->SetOptStat(0);
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetZaxis()->SetRangeUser(0,250);
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetYaxis()->SetRangeUser(0.7,20);
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->GetXaxis()->SetRangeUser(0.7,10);
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Draw("colz");
    TF1 *fit_TEC = new TF1("fit_TEC","[0]+[1]/x",0,10);
    TEC_MaxOverMaxNeighbour_VS_MaxNeihbourRatio->Fit("fit_TEC");
*/

/*  TCanvas *c_trash = new TCanvas("c_trash", "c_trash", 100, 100);
    c_trash->cd();
    std::vector <float> p0;
    std::vector <float> p1;
    std::vector <TString> eq_fit;
    std::vector <TF1*> func_fit;
    for (int i=0; i<histoVector_IonizationNeighbour.size(); i++)
    {
        TF1 *fit = new TF1("fit","[0]+[1]/x",0,10);
        histoVector_IonizationNeighbour[i]->Fit("fit","REM");

        p0.push_back(fit->GetParameter(0));
        p1.push_back(fit->GetParameter(1));
        eq_fit.push_back(Form("%f + %f/x", p0[i], p1[i]));

        TF1 *func = new TF1("func", eq_fit[i], 0, 10);
        func->SetLineColor(i);
        if (i>=5 && i<=10) func->SetLineColor(i-4);
        if (i>=11 && i<=13) func->SetLineColor(i-10);
        if (i>=14 && i<=22) func->SetLineColor(i-13);
        if (i>=23 && i<=29) func->SetLineColor(i-22);
        func->SetLineWidth(2);
        func_fit.push_back(func);
    }
    for (int i=0; i<eq_fit.size(); i++) cout<< "Fit of histoVector_IonizationNeighbour["<<i<<"]  : y= "<<eq_fit[i]<<endl;

    func_fit[0]->SetTitle("No cleaning");
    func_fit[0]->GetXaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    func_fit[0]->GetYaxis()->SetTitle("Max / 1^{st} Max neighbour");
    func_fit[0]->GetXaxis()->SetRangeUser(0,10);
    func_fit[0]->GetYaxis()->SetRangeUser(0,20);
    
    TCanvas *c_fit1_TIB = new TCanvas("c_fit1_TIB", "c_fit1_TIB", 1100, 700);
    c_fit1_TIB->cd();
    gStyle->SetOptStat(0);
    func_fit[0]->Draw();
    for (int i=1; i<=4; i++) func_fit[i]->Draw("same");
    auto l_TIB= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=1; i<=4; i++) l_TIB->AddEntry(func_fit[i],Form("TIB layer %d",i),"l");
    l_TIB->SetBorderSize(0);
    l_TIB->Draw();

    TCanvas *c_fit1_TOB = new TCanvas("c_fit1_TOB", "c_fit1_TOB", 1100, 700);
    c_fit1_TOB->cd();
    gStyle->SetOptStat(0);
    func_fit[0]->Draw();
    for (int i=5; i<=10; i++) func_fit[i]->Draw("same");
    auto l_TOB= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=5; i<=10; i++) l_TOB->AddEntry(func_fit[i],Form("TOB layer %d",i-4),"l");
    l_TOB->SetBorderSize(0);
    l_TOB->Draw();

    TCanvas *c_fit1_TID = new TCanvas("c_fit1_TID", "c_fit1_TID", 1100, 700);
    c_fit1_TID->cd();
    gStyle->SetOptStat(0);
    func_fit[0]->Draw();
    for (int i=11; i<=13; i++) func_fit[i]->Draw("same");
    auto l_TID= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=11; i<=13; i++) l_TID->AddEntry(func_fit[i],Form("TID ring %d",i-10),"l");
    l_TID->SetBorderSize(0);
    l_TID->Draw();

    TCanvas *c_fit1_TEC_wheel = new TCanvas("c_fit1_TEC_wheel", "c_fit1_TEC_wheel", 1100, 700);
    c_fit1_TEC_wheel->cd();
    gStyle->SetOptStat(0);
    func_fit[0]->Draw();
    for (int i=14; i<=22; i++) func_fit[i]->Draw("same");
    auto l_TEC= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=14; i<=22; i++) l_TEC->AddEntry(func_fit[i],Form("TEC wheel %d",i-13),"l");
    l_TEC->SetBorderSize(0);
    l_TEC->Draw();

    TCanvas *c_fit1_TEC_ring = new TCanvas("c_fit1_TEC_ring", "c_fit1_TEC_ring", 1100, 700);
    c_fit1_TEC_ring->cd();
    gStyle->SetOptStat(0);
    func_fit[0]->Draw();
    for (int i=23; i<=29; i++) func_fit[i]->Draw("same");
    auto l_TEC_ring= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=23; i<=29; i++) l_TEC_ring->AddEntry(func_fit[i],Form("TEC ring %d",i-22),"l");
    l_TEC_ring->SetBorderSize(0);
    l_TEC_ring->Draw();


    TCanvas *c_trash2 = new TCanvas("c_trash2", "c_trash2", 100, 100);
    c_trash2->cd();
    std::vector <float> x0;
    std::vector <float> x1;
    std::vector <float> x2;
    std::vector <float> x3;
    std::vector <TString> eq_fit_b;
    std::vector <TString> eq_fit_b2;
    std::vector <TF1*> func_fit_b;
    std::vector <TF1*> func_fit_b2;
    for (int i=0; i<histoVector_MaxNeighbourBarycenter.size(); i++)
    {
        TF1 *fit = new TF1("fit","[0]+[1]*x",0,0.5);
        TF1 *fit2 = new TF1("fit2","[0]+[1]*x",-0.5,0);
        histoVector_MaxNeighbourBarycenter[i]->Fit("fit","REM");
        histoVector_MaxNeighbourBarycenter[i]->Fit("fit2","REM");

        x0.push_back(fit->GetParameter(0));
        x1.push_back(fit->GetParameter(1));
        x2.push_back(fit2->GetParameter(0));
        x3.push_back(fit2->GetParameter(1));
        eq_fit_b.push_back(Form("%f + %f*x", x0[i], x1[i]));
        eq_fit_b2.push_back(Form("%f + %f*x", x2[i], x3[i]));

        TF1 *func = new TF1("func", eq_fit_b[i], -3, 3);
        TF1 *func2 = new TF1("func2", eq_fit_b2[i], -3, 3);
        func->SetLineColor(i);
        func2->SetLineColor(i);
        if (i>=5 && i<=10) func->SetLineColor(i-4);
        if (i>=11 && i<=13) func->SetLineColor(i-10);
        if (i>=14 && i<=22) func->SetLineColor(i-13);
        if (i>=23 && i<=29) func->SetLineColor(i-22);
        if (i>=5 && i<=10) func2->SetLineColor(i-4);
        if (i>=11 && i<=13) func2->SetLineColor(i-10);
        if (i>=14 && i<=22) func2->SetLineColor(i-13);
        if (i>=23 && i<=29) func2->SetLineColor(i-22);
        func->SetLineWidth(2);
        func2->SetLineWidth(2);
        func_fit_b.push_back(func);
        func_fit_b2.push_back(func2);
    }
    for (int i=0; i<histoVector_MaxNeighbourBarycenter.size(); i++) cout<< "[0,+3] Fit of histoVector_MaxNeighbourBarycenter["<<i<<"]  : y= "<<eq_fit_b[i]<<endl;
    for (int i=0; i<histoVector_MaxNeighbourBarycenter.size(); i++) cout<< "[-3,0] Fit of histoVector_MaxNeighbourBarycenter["<<i<<"]  : y= "<<eq_fit_b2[i]<<endl;

    func_fit_b[0]->SetTitle("No cleaning");
    func_fit_b[0]->GetXaxis()->SetTitle("i_{max} - barycenter");
    func_fit_b[0]->GetYaxis()->SetTitle("1^{st} Max neighbour / 1^{st} Min neighbour");
    func_fit_b[0]->GetXaxis()->SetRangeUser(-2,3);
    func_fit_b[0]->GetYaxis()->SetRangeUser(0.8,10);
    
    TCanvas *c_fit2_TIB = new TCanvas("c_fit2_TIB", "c_fit2_TIB", 1100, 700);
    c_fit2_TIB->cd();
    gStyle->SetOptStat(0);
    func_fit_b[0]->Draw();
    for (int i=1; i<=4; i++)
    {
        func_fit_b[i]->Draw("same");
        func_fit_b2[i]->Draw("same");
    }
    auto l_TIB2= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=1; i<=4; i++) l_TIB2->AddEntry(func_fit_b[i],Form("TIB layer %d",i),"l");
    l_TIB2->SetBorderSize(0);
    l_TIB2->Draw();

    TCanvas *c_fit2_TOB = new TCanvas("c_fit2_TOB", "c_fit2_TOB", 1100, 700);
    c_fit2_TOB->cd();
    gStyle->SetOptStat(0);
    func_fit_b[0]->Draw();
    for (int i=5; i<=10; i++)
    {
        func_fit_b[i]->Draw("same");
        func_fit_b2[i]->Draw("same");
    }
    auto l_TOB2= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=5; i<=10; i++) l_TOB2->AddEntry(func_fit_b[i],Form("TOB layer %d",i-4),"l");
    l_TOB2->SetBorderSize(0);
    l_TOB2->Draw();

    TCanvas *c_fit2_TID = new TCanvas("c_fit2_TID", "c_fit2_TID", 1100, 700);
    c_fit2_TID->cd();
    gStyle->SetOptStat(0);
    func_fit_b[0]->Draw();
    for (int i=11; i<=13; i++)
    {
        func_fit_b[i]->Draw("same");
        func_fit_b2[i]->Draw("same");
    }
    auto l_TID2= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=11; i<=13; i++) l_TID2->AddEntry(func_fit_b[i],Form("TID ring %d",i-10),"l");
    l_TID2->SetBorderSize(0);
    l_TID2->Draw();

    TCanvas *c_fit2_TEC_wheel = new TCanvas("c_fit2_TEC_wheel", "c_fit2_TEC_wheel", 1100, 700);
    c_fit2_TEC_wheel->cd();
    gStyle->SetOptStat(0);
    func_fit_b[0]->Draw();
    for (int i=14; i<=22; i++)
    {
        func_fit_b[i]->Draw("same");
        func_fit_b2[i]->Draw("same");
    }
    auto l_TEC2= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=14; i<=22; i++) l_TEC2->AddEntry(func_fit_b[i],Form("TEC wheel %d",i-13),"l");
    l_TEC2->SetBorderSize(0);
    l_TEC2->Draw();

    TCanvas *c_fit2_TEC_ring = new TCanvas("c_fit2_TEC_ring", "c_fit2_TEC_ring", 1100, 700);
    c_fit2_TEC_ring->cd();
    gStyle->SetOptStat(0);
    func_fit_b[0]->Draw();
    for (int i=23; i<=29; i++)
    {
        func_fit_b[i]->Draw("same");
        func_fit_b2[i]->Draw("same");
    }
    auto l_TEC_ring2= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=23; i<=29; i++) l_TEC_ring2->AddEntry(func_fit[i],Form("TEC ring %d",i-22),"l");
    l_TEC_ring2->SetBorderSize(0);
    l_TEC_ring2->Draw();



    TCanvas *c_Qloss=new TCanvas("c_Qloss","c_Qloss",1100,700);
    c_Qloss->cd();
    gStyle->SetOptStat(1111111);
    Compare_eloss_over_sigdigi->Draw("hist");
    Compare_eloss_over_sigdigi_nosat->Draw("hist same");
    Compare_eloss_over_sigdigi->SetAxisRange(0,2,"X");
    auto l1= new TLegend(0.6,0.35,0.85,0.55);
    l1->AddEntry(Compare_eloss_over_sigdigi,Form("All clusters, #mu=%3.2f  #sigma=%4.3f",Compare_eloss_over_sigdigi->GetMean(),Compare_eloss_over_sigdigi->GetStdDev()),"l");
    l1->AddEntry(Compare_eloss_over_sigdigi_nosat,Form("Only ones not saturated, #mu=%3.2f  #sigma=%4.3f",Compare_eloss_over_sigdigi_nosat->GetMean(),Compare_eloss_over_sigdigi_nosat->GetStdDev()),"l");
    l1->SetBorderSize(0);
    l1->Draw();*/
    

    TCanvas *c_QII_Over=new TCanvas("c_QII_Over","c_QII_Over",1100,700);
    c_QII_Over->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Over_SigDigi_FL->SetTitle("No cleaning - Crosstalk correction");
    Compare_QII_Over_SigDigi_FL->GetXaxis()->SetTitle("Q_{corr}/Q_{sum}");
    Compare_QII_Over_SigDigi_FL->GetYaxis()->SetTitle("Nbr of entries (normalized)");
    //Compare_QII_Over_SigDigi->Draw("hist");
    Compare_QII_Over_SigDigi_FL->Draw("hist");
    //Compare_QII_Over_SigDigi_L->Draw("hist same");
    //Compare_QII_Over_SigDigi_R->Draw("hist same");
    //Compare_QII_Over_SigDigi_C->Draw("hist same");
    Compare_QII_Over_SigDigi_FR->Draw("hist same");
    auto lb= new TLegend(0.6,0.35,0.85,0.55);
    //lb->AddEntry(Compare_QII_Over_SigDigi_L,Form("Left: #sigma= %4.2f",Compare_QII_Over_SigDigi_L->GetStdDev()),"l");
    //lb->AddEntry(Compare_QII_Over_SigDigi_R,Form("Right: #sigma= %4.2f",Compare_QII_Over_SigDigi_R->GetStdDev()),"l");
    //lb->AddEntry(Compare_QII_Over_SigDigi_C,Form("Center: #sigma= %4.2f",Compare_QII_Over_SigDigi_C->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FL,Form("FullLeft: #sigma= %4.2f",Compare_QII_Over_SigDigi_FL->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FR,Form("FullRight: #sigma= %4.2f",Compare_QII_Over_SigDigi_FR->GetStdDev()),"l");
    lb->SetBorderSize(0);
    lb->Draw();
    
    /*TCanvas *c_QII_Minus=new TCanvas("c_QII_Minus","c_QII_Minus",1100,700);
    c_QII_Minus->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Minus_SigDigi_FR->SetTitle("No cleaning");
    Compare_QII_Minus_SigDigi_FR->GetXaxis()->SetTitle("Q_{corr} - Q_{sum}");
    Compare_QII_Minus_SigDigi_FR->GetYaxis()->SetTitle("Nbr of entries (normalized)");
    //Compare_QII_Minus_SigDigi->Draw("hist");
    Compare_QII_Minus_SigDigi_FR->Draw("hist");
    Compare_QII_Minus_SigDigi_L->Draw("hist same");
    Compare_QII_Minus_SigDigi_R->Draw("hist same");
    Compare_QII_Minus_SigDigi_C->Draw("hist same");
    Compare_QII_Minus_SigDigi_FL->Draw("hist same");
    Compare_QII_Minus_SigDigi_FR->SetAxisRange(0,0.5,"Y");
    auto la= new TLegend(0.6,0.35,0.85,0.55);
    la->AddEntry(Compare_QII_Minus_SigDigi_L,Form("Left: #sigma= %4.2f",Compare_QII_Minus_SigDigi_L->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_R,Form("Right: #sigma= %4.2f",Compare_QII_Minus_SigDigi_R->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_C,Form("Center: #sigma= %4.2f",Compare_QII_Minus_SigDigi_C->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FL,Form("FullLeft: #sigma= %4.2f",Compare_QII_Minus_SigDigi_FL->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FR,Form("FullRight: #sigma= %4.2f",Compare_QII_Minus_SigDigi_FR->GetStdDev()),"l");
    la->SetBorderSize(0);
    la->Draw();*/


// -----------------------------------------------------------------

    /*TCanvas *c_MaxOverVmax_Barycenter=new TCanvas("c_MaxOverVmax_Barycenter","c_MaxOverVmax_Barycenter",1100,700);
    c_MaxOverVmax_Barycenter->cd();
    MaxOverMaxNeighbour_VS_Barycenter->Draw("colz");
    MaxOverMaxNeighbour_VS_Barycenter->GetZaxis()->SetRangeUser(0,250);
    gStyle->SetOptStat(0);

    TCanvas *c_trash3 = new TCanvas("c_trash3", "c_trash3", 100, 100);
    c_trash3->cd();
    std::vector <float> a;
    std::vector <float> b;
    std::vector <float> bp;
    std::vector <float> c;
    std::vector <float> d;
    std::vector <float> dp;
    std::vector <TString> eq_fit_1;
    std::vector <TString> eq_fit_2;
    std::vector <TF1*> func_fit_1;
    std::vector <TF1*> func_fit_2;
    for (int i=0; i<histoVector_MaxVmaxBarycenter.size(); i++)
    {
        //TF1 *fit = new TF1("fit","[0]*TMath::Exp([1]*x) + 1",0.03,0.4);
        //TF1 *fit2 = new TF1("fit2","[0]*TMath::Exp([1]*x) + 1",-0.4,-0.03);
        TF1 *fit = new TF1("fit","[0] + [1]*TMath::Abs(x)",-0.4,0.4);
        //TF1 *fit2 = new TF1("fit2","pol1",-0.4,0);
        histoVector_MaxVmaxBarycenter[i]->Fit("fit","REM");
        //histoVector_MaxVmaxBarycenter[i]->Fit("fit2","REM");

        a.push_back(fit->GetParameter(0));
        b.push_back(fit->GetParameter(1));
        //c.push_back(fit2->GetParameter(0));
        //d.push_back(fit2->GetParameter(1));
        //if (i != 1 && i != 7 && i!=10)
        //{
            //eq_fit_1.push_back(Form("%f * TMath::Exp(%f *x) + 1", a[i], b[i]));
            //eq_fit_2.push_back(Form("%f * TMath::Exp(%f *x) + 1", c[i], d[i]));
            eq_fit_1.push_back(Form("%f + %f * TMath::Abs(x)", a[i], b[i]));
            //eq_fit_2.push_back(Form("%f * x + %f", d[i], c[i]));
        //}
        //else
        //{
            //eq_fit_1.push_back(Form("%f * TMath::Exp(-%f *x) + 1", c[i], d[i]));
            //eq_fit_2.push_back(Form("%f * TMath::Exp(%f *x) + 1", c[i], d[i]));
        //}
        if (i == 0)
        {
            TF1 *func = new TF1("func", eq_fit_1[i], -3, 3);
            func->SetLineColor(0);
            func_fit_1.push_back(func);
            func_fit_2.push_back(func);
        }
        else 
        {
            TF1 *func = new TF1("func", eq_fit_1[i], -3, 3);
            //TF1 *func2 = new TF1("func2", eq_fit_2[i], -3, 0);
            func->SetLineColor(i);
            //func2->SetLineColor(i);
            if (i>=5 && i<=10) func->SetLineColor(i-4);
            if (i>=11 && i<=13) func->SetLineColor(i-10);
            if (i>=14 && i<=20) func->SetLineColor(i-13);
            //if (i>=5 && i<=10) func2->SetLineColor(i-4);
            //if (i>=11 && i<=13) func2->SetLineColor(i-10);
            //if (i>=14 && i<=20) func2->SetLineColor(i-13);
            func->SetLineWidth(2);
            //func2->SetLineWidth(2);
            func_fit_1.push_back(func);
            func_fit_2.push_back(func);
        }
    }
    for (int i=0; i<histoVector_MaxVmaxBarycenter.size(); i++) cout<< "[0,+3] Fit of histoVector_MaxVmaxBarycenter["<<i<<"]  : y= "<<eq_fit_1[i]<<endl;
    //for (int i=0; i<histoVector_MaxVmaxBarycenter.size(); i++) cout<< "[-3,0] Fit of histoVector_MaxVmaxBarycenter["<<i<<"]  : y= "<<eq_fit_2[i]<<endl;

    func_fit_1[0]->SetTitle("No cleaning");
    func_fit_1[0]->GetXaxis()->SetTitle("i_{max} - barycenter");
    func_fit_1[0]->GetYaxis()->SetTitle("Max / V_{max}");
    func_fit_1[0]->GetXaxis()->SetRangeUser(-3,3);
    func_fit_1[0]->GetYaxis()->SetRangeUser(0.8,20);
    
    TCanvas *c_fit3_TIB = new TCanvas("c_fit3_TIB", "c_fit3_TIB", 1100, 700);
    c_fit3_TIB->cd();
    gStyle->SetOptStat(0);
    func_fit_1[0]->Draw();
    for (int i=1; i<=4; i++)
    {
        func_fit_1[i]->Draw("same");
        //func_fit_2[i]->Draw("same");
    }
    auto l_TIB3= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=1; i<=4; i++) l_TIB3->AddEntry(func_fit_1[i],Form("TIB layer %d",i),"l");
    l_TIB3->SetBorderSize(0);
    l_TIB3->Draw();

    TCanvas *c_fit3_TOB = new TCanvas("c_fit3_TOB", "c_fit3_TOB", 1100, 700);
    c_fit3_TOB->cd();
    gStyle->SetOptStat(0);
    func_fit_1[0]->Draw();
    for (int i=5; i<=10; i++)
    {
        func_fit_1[i]->Draw("same");
        //func_fit_2[i]->Draw("same");
    }
    auto l_TOB3= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=5; i<=10; i++) l_TOB3->AddEntry(func_fit_1[i],Form("TOB layer %d",i-4),"l");
    l_TOB3->SetBorderSize(0);
    l_TOB3->Draw();

    TCanvas *c_fit3_TID = new TCanvas("c_fit3_TID", "c_fit3_TID", 1100, 700);
    c_fit3_TID->cd();
    gStyle->SetOptStat(0);
    func_fit_1[0]->Draw();
    for (int i=11; i<=13; i++)
    {
        func_fit_1[i]->Draw("same");
        //func_fit_2[i]->Draw("same");
    }
    auto l_TID3= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=11; i<=13; i++) l_TID3->AddEntry(func_fit_1[i],Form("TID ring %d",i-10),"l");
    l_TID3->SetBorderSize(0);
    l_TID3->Draw();

    TCanvas *c_fit3_TEC = new TCanvas("c_fit3_TEC", "c_fit3_TEC", 1100, 700);
    c_fit3_TEC->cd();
    gStyle->SetOptStat(0);
    func_fit_1[0]->Draw();
    for (int i=14; i<=20; i++)
    {
        func_fit_1[i]->Draw("same");
        //func_fit_2[i]->Draw("same");
    }
    auto l_TEC3= new TLegend(0.65,0.4,0.85,0.85);
    for (int i=14; i<=20; i++) l_TEC3->AddEntry(func_fit_1[i],Form("TEC ring %d",i-13),"l");
    l_TEC3->SetBorderSize(0);
    l_TEC3->Draw();
    

    TCanvas *c_check = new TCanvas("c_check", "c_check", 1100, 700);
    c_check->cd();
    gStyle->SetOptStat(0);
    histoVector_MaxVmaxBarycenter[13]->Draw("colz");
    histoVector_MaxVmaxBarycenter[13]->GetXaxis()->SetTitle("i_{max} - barycenter");
    histoVector_MaxVmaxBarycenter[13]->GetYaxis()->SetTitle("Max / V_{max}");
    histoVector_MaxVmaxBarycenter[13]->SetTitle("No cleaning");
    histoVector_MaxVmaxBarycenter[13]->GetZaxis()->SetRangeUser(0,150);
    func_fit_1[13]->Draw("same");
    //func_fit_2[12]->Draw("same");

    TCanvas *c_ClusterUnderCorrection = new TCanvas("c_ClusterUnderCorrection", "c_ClusterUnderCorrection", 1100, 700);
    c_ClusterUnderCorrection->cd();
    gStyle->SetOptStat(0);
    ClusterUnderCorrection_layer_B->Scale(1./610980);
    ClusterUnderCorrection_layer->Scale(1./529529);
    ClusterUnderCorrection_layer_BandN->Scale(1./529529);
    ClusterUnderCorrection_layer_B->Draw();
    ClusterUnderCorrection_layer_B->SetTitle("No cleaning");
    ClusterUnderCorrection_layer->Draw("same");
    ClusterUnderCorrection_layer_BandN->Draw("same");
    auto l_ClusterUnderCorrection = new TLegend(0.45,0.6,0.89,0.85);
    l_ClusterUnderCorrection->AddEntry(ClusterUnderCorrection_layer_B,Form("wBarycenter: 269812 / 610980 (%3.1f %%)",(double)100 * 269812/610980),"l");
    l_ClusterUnderCorrection->AddEntry(ClusterUnderCorrection_layer,Form("wNeighbours: 97435 / 592529 (%3.1f %%)",(double)100 * 97435/592529),"l");
    l_ClusterUnderCorrection->AddEntry(ClusterUnderCorrection_layer_BandN,Form("wBarycenter AND wNeighbours: 96659 / 592529 (%3.1f %%)",(double)100 * 96659/592529),"l");
    l_ClusterUnderCorrection->SetBorderSize(0);
    l_ClusterUnderCorrection->Draw();*/

    TCanvas *c_ClusterUnderCorrection_sizeVSeta= new TCanvas("c_ClusterUnderCorrection_sizeVSeta", "c_ClusterUnderCorrection_sizeVSeta", 1100, 700);
    c_ClusterUnderCorrection_sizeVSeta->cd();
    gStyle->SetOptStat(0);
    ClusterUnderCorrection_sizeVSeta->Draw("colz");

    TCanvas *c_ClusterGoodCorrection_sizeVSeta= new TCanvas("c_ClusterGoodCorrection_sizeVSeta", "c_ClusterGoodCorrection_sizeVSeta", 1100, 700);
    c_ClusterGoodCorrection_sizeVSeta->cd();
    gStyle->SetOptStat(0);
    ClusterGoodCorrection_sizeVSeta->Draw("colz");

    TCanvas *c_ClusterUnderCorrection_Neighbours = new TCanvas("c_ClusterUnderCorrection_Neighbours", "c_ClusterUnderCorrection_Neighbours", 1100, 700);
    c_ClusterUnderCorrection_Neighbours->cd();
    gStyle->SetOptStat(0);
    ClusterUnderCorrection_layer->Divide(ClusterUnderCorrection_layer_all);
    ClusterUnderCorrection_layer->Draw("E1");

    TCanvas *c_shape_TIBL2 = new TCanvas("c_shape_TIBL2", "c_shape_TIBL2", 1100, 700);
    c_shape_TIBL2->cd();
    gStyle->SetOptStat(0);
    MaxcorrOverMax_VS_eta_TIBL2->GetZaxis()->SetRangeUser(0,200);
    MaxcorrOverMax_VS_eta_TIBL2->Draw("colz");

    TCanvas *c_shape_TIDR2 = new TCanvas("c_shape_TIDR2", "c_shape_TIDR2", 1100, 700);
    c_shape_TIDR2->cd();
    gStyle->SetOptStat(0);
    MaxcorrOverMax_VS_shape_TIDR2->GetZaxis()->SetRangeUser(0,200);
    MaxcorrOverMax_VS_shape_TIDR2->Draw("colz");
    
    

// -----------------------------------------------------------------

/*  TCanvas *c_Over5 = new TCanvas("c_Over5", "c_Over5", 1100, 700);
    c_Over5->cd();
    gStyle->SetOptStat(0);
    StripSat_Over5_layer_eta->Draw("colz");

    TCanvas *c_Qloss=new TCanvas("c_Qloss","c_Qloss",1100,700);
    c_Qloss->cd();
    gStyle->SetOptStat(1111111);
    Compare_eloss_over_sigdigi->Draw("hist");
    Compare_eloss_over_sigdigi_nosat->Draw("hist same");
    Compare_eloss_over_sigdigi->SetAxisRange(0,2,"X");
    auto l1= new TLegend(0.6,0.35,0.85,0.55);
    l1->AddEntry(Compare_eloss_over_sigdigi,Form("All clusters, #mu=%3.2f  #sigma=%4.3f",Compare_eloss_over_sigdigi->GetMean(),Compare_eloss_over_sigdigi->GetStdDev()),"l");
    l1->AddEntry(Compare_eloss_over_sigdigi_nosat,Form("Only ones not saturated, #mu=%3.2f  #sigma=%4.3f",Compare_eloss_over_sigdigi_nosat->GetMean(),Compare_eloss_over_sigdigi_nosat->GetStdDev()),"l");
    l1->SetBorderSize(0);
    l1->Draw();
    

    TCanvas *c_QII_Over=new TCanvas("c_QII_Over","c_QII_Over",1100,700);
    c_QII_Over->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Over_SigDigi_FL->SetTitle("No cleaning");
    Compare_QII_Over_SigDigi_FL->GetXaxis()->SetTitle("Q_{corr}/Q_{sum}");
    Compare_QII_Over_SigDigi_FL->GetYaxis()->SetTitle("Nbr of entries (normalized)");
    //Compare_QII_Over_SigDigi->Draw("hist");
    Compare_QII_Over_SigDigi_FL->Draw("hist");
    Compare_QII_Over_SigDigi_L->Draw("hist same");
    Compare_QII_Over_SigDigi_R->Draw("hist same");
    Compare_QII_Over_SigDigi_C->Draw("hist same");
    Compare_QII_Over_SigDigi_FR->Draw("hist same");
    auto lb= new TLegend(0.6,0.35,0.85,0.55);
    lb->AddEntry(Compare_QII_Over_SigDigi_L,Form("Left: #sigma= %4.2f",Compare_QII_Over_SigDigi_L->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_R,Form("Right: #sigma= %4.2f",Compare_QII_Over_SigDigi_R->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_C,Form("Center: #sigma= %4.2f",Compare_QII_Over_SigDigi_C->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FL,Form("FullLeft: #sigma= %4.2f",Compare_QII_Over_SigDigi_FL->GetStdDev()),"l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FR,Form("FullRight: #sigma= %4.2f",Compare_QII_Over_SigDigi_FR->GetStdDev()),"l");
    lb->SetBorderSize(0);
    lb->Draw();
    
    TCanvas *c_QII_Minus=new TCanvas("c_QII_Minus","c_QII_Minus",1100,700);
    c_QII_Minus->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Minus_SigDigi_FR->SetTitle("No cleaning");
    Compare_QII_Minus_SigDigi_FR->GetXaxis()->SetTitle("Q_{corr} - Q_{sum}");
    Compare_QII_Minus_SigDigi_FR->GetYaxis()->SetTitle("Nbr of entries (normalized)");
    //Compare_QII_Minus_SigDigi->Draw("hist");
    Compare_QII_Minus_SigDigi_FR->Draw("hist");
    Compare_QII_Minus_SigDigi_L->Draw("hist same");
    Compare_QII_Minus_SigDigi_R->Draw("hist same");
    Compare_QII_Minus_SigDigi_C->Draw("hist same");
    Compare_QII_Minus_SigDigi_FL->Draw("hist same");
    Compare_QII_Minus_SigDigi_FR->SetAxisRange(0,0.5,"Y");
    auto la= new TLegend(0.6,0.35,0.85,0.55);
    la->AddEntry(Compare_QII_Minus_SigDigi_L,Form("Left: #sigma= %4.2f",Compare_QII_Minus_SigDigi_L->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_R,Form("Right: #sigma= %4.2f",Compare_QII_Minus_SigDigi_R->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_C,Form("Center: #sigma= %4.2f",Compare_QII_Minus_SigDigi_C->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FL,Form("FullLeft: #sigma= %4.2f",Compare_QII_Minus_SigDigi_FL->GetStdDev()),"l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FR,Form("FullRight: #sigma= %4.2f",Compare_QII_Minus_SigDigi_FR->GetStdDev()),"l");
    la->SetBorderSize(0);
    la->Draw();*/


/*
    TCanvas *c_check_L=new TCanvas("c_check_L","c_check_L",1100,700);
    c_check_L->cd();
    gStyle->SetOptStat(0);
    check_C_Original->SetTitle("1st highest neighbor of the max in %");
    check_C_Original->Draw("hist");
    check_FL_Original->Draw("hist same");
    check_L_Original->Draw("hist same");
    check_R_Original->Draw("hist same");
    check_FR_Original->Draw("hist same");
    auto lc= new TLegend(0.6,0.35,0.85,0.55);
    lc->AddEntry(check_L_Original,"Left","l");
    lc->AddEntry(check_R_Original,"Right","l");
    lc->AddEntry(check_C_Original,"Center (only left)","l");
    lc->AddEntry(check_FL_Original,"Full Left","l");
    lc->AddEntry(check_FR_Original,"Full Right","l");
    lc->SetBorderSize(0);
    lc->Draw();
    
    
    TCanvas *c_Ih_VS_p_1=new TCanvas("c_Ih_VS_p_1","c_Ih_VS_p_1",1100,700);
    c_Ih_VS_p_1->cd();
    gStyle->SetOptStat(0);
    Ih_VS_p_1->Draw("colz");


    TCanvas *c_Ih_VS_p_2=new TCanvas("c_Ih_VS_p_2","c_Ih_VS_p_2",1100,700);
    c_Ih_VS_p_2->cd();
    gStyle->SetOptStat(0);
    Ih_VS_p_2->Draw("colz");

    TCanvas *c_FullLeft_Border=new TCanvas("c_FullLeft_Border","c_FullLeft_Border",1100,700);
    c_FullLeft_Border->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullLeft_Border->Draw("hist text");
    line[7]->Draw("same");
    l3->Draw();

    TCanvas *c_FullRight_Border=new TCanvas("c_FullRight_Border","c_FullRight_Border",1100,700);
    c_FullRight_Border->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullRight_Border->Draw("hist text");
    line[8]->Draw("same");
    l3->Draw();

    TCanvas *c_FullLeft_noBorder=new TCanvas("c_FullLeft_noBorder","c_FullLeft_noBorder",1100,700);
    c_FullLeft_noBorder->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullLeft_noBorder->Draw("hist text");
    line[5]->Draw("same");
    l3->Draw();

    TCanvas *c_FullRight_noBorder=new TCanvas("c_FullRight_noBorder","c_FullRight_noBorder",1100,700);
    c_FullRight_noBorder->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullRight_noBorder->Draw("hist text");
    line[6]->Draw("same");
    l3->Draw();


    TCanvas *c_FullRight_layer=new TCanvas("c_FullRight_layer","c_FullRight_layer",1100,700);
    c_FullRight_layer->cd();
    gStyle->SetOptStat(0);
    Sigdigi_layer_maxsat->Draw("hist");
    Sigdigi_FullLeft_layer->Draw("hist SAME");
    Sigdigi_FullRight_layer->Draw("hist SAME");
    auto l2= new TLegend(0.6,0.6,0.85,0.85);
    l2->AddEntry(Sigdigi_layer_maxsat,"All saturated clusters","l");
    l2->AddEntry(Sigdigi_FullLeft_layer,"Ones with max on left","l");
    l2->AddEntry(Sigdigi_FullRight_layer,"Ones with max on right","l");
    l2->SetBorderSize(0);
    l2->Draw();
    */
    
}

void histo_MuonPU_2strips()
{
    TFile *Sat_MuonPU_alltracks = new TFile("ROOT_SVG/Sat_MuonPU_twoStrips_alltracks.root", "READ");

    TH1F *Sigdigi_left_1 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_left_1");
    TH1F *Sigdigi_left_2 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_left_2");
    TH1F *Sigdigi_right_1 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_right_1");
    TH1F *Sigdigi_right_2 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_right_2");
    TH1F *Sigdigi_center_1 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_center_1");
    TH1F *Sigdigi_center_2 = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_center_2");
    TH1F *Sigdigi_FullLeft = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullLeft");
    TH1F *Sigdigi_FullRight = (TH1F*)Sat_MuonPU_alltracks->Get("Sigdigi_FullRight");
    TH1F *cpt_sat = (TH1F*)Sat_MuonPU_alltracks->Get("cpt_sat");

    TH1F *Compare_QII_Over_SigDigi = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi");
    TH1F *Compare_QII_Over_SigDigi_L = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_L");
    TH1F *Compare_QII_Over_SigDigi_R = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_R");
    TH1F *Compare_QII_Over_SigDigi_C = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_C");
    TH1F *Compare_QII_Over_SigDigi_FL = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_FL");
    TH1F *Compare_QII_Over_SigDigi_FR = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Over_SigDigi_FR");

    TH1F *Compare_QII_Minus_SigDigi = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi");
    TH1F *Compare_QII_Minus_SigDigi_L = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_L");
    TH1F *Compare_QII_Minus_SigDigi_R = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_R");
    TH1F *Compare_QII_Minus_SigDigi_C = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_C");
    TH1F *Compare_QII_Minus_SigDigi_FL = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_FL");
    TH1F *Compare_QII_Minus_SigDigi_FR = (TH1F*)Sat_MuonPU_alltracks->Get("Compare_QII_Minus_SigDigi_FR");

    vector <TLine*> line;
    for (int i=1;i<=8;i++)
    {
        TLine *l1 = new TLine(0, 255*cpt_sat->GetBinContent(i), 50, 255*cpt_sat->GetBinContent(i));
        l1->SetLineColor(kRed);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        line.push_back(l1);
    }

    auto l3= new TLegend(0.7,0.7,0.83,0.83);
    l3->AddEntry(line[0],"255 ADC","l");
    l3->SetBorderSize(0);

    Sigdigi_left_1->SetLineColor(kBlack);
    Sigdigi_left_1->SetLineWidth(2);
    Sigdigi_left_1->SetTitle("1st left neighbor/1st right neighbor > 10\% && 2nd max on left");
    Sigdigi_left_2->SetLineColor(kBlack);
    Sigdigi_left_2->SetLineWidth(2);
    Sigdigi_left_2->SetTitle("1st left neighbor/1st right neighbor > 10\% && 2nd max on right");
    Sigdigi_right_1->SetLineColor(kBlack);
    Sigdigi_right_1->SetLineWidth(2);
    Sigdigi_right_1->SetTitle("1st right neighbor/1st left neighbor > 10\% && 2nd max on left");
    Sigdigi_right_2->SetLineColor(kBlack);
    Sigdigi_right_2->SetLineWidth(2);
    Sigdigi_right_2->SetTitle("1st right neighbor/1st left neighbor > 10\% && 2nd max on right");
    Sigdigi_center_1->SetLineColor(kBlack);
    Sigdigi_center_1->SetLineWidth(2);
    Sigdigi_center_1->SetTitle("|1- 1st left neighbor/1st right neighbor|<10\% && 2nd max on left");
    Sigdigi_center_2->SetLineColor(kBlack);
    Sigdigi_center_2->SetLineWidth(2);
    Sigdigi_center_2->SetTitle("|1- 1st left neighbor/1st right neighbor|<10\% && 2nd max on right");
    Sigdigi_FullLeft->SetLineColor(kBlack);
    Sigdigi_FullLeft->SetLineWidth(2);
    Sigdigi_FullLeft->SetTitle("Max on left");
    Sigdigi_FullRight->SetLineColor(kBlack);
    Sigdigi_FullRight->SetLineWidth(2);
    Sigdigi_FullRight->SetTitle("Max on right");
    Compare_QII_Over_SigDigi->SetLineColor(kBlack);
    Compare_QII_Over_SigDigi->SetLineWidth(2);
    Compare_QII_Over_SigDigi->SetTitle("QII/Sum_sigdigi");
    Compare_QII_Over_SigDigi_L->SetLineColor(kBlack);
    Compare_QII_Over_SigDigi_L->SetLineWidth(2);
    Compare_QII_Over_SigDigi_L->Scale(1./Compare_QII_Over_SigDigi_L->Integral());
    Compare_QII_Over_SigDigi_R->SetLineColor(kRed);
    Compare_QII_Over_SigDigi_R->SetLineWidth(2);
    Compare_QII_Over_SigDigi_R->Scale(1./Compare_QII_Over_SigDigi_R->Integral());
    Compare_QII_Over_SigDigi_C->SetLineColor(kOrange);
    Compare_QII_Over_SigDigi_C->SetLineWidth(2);
    Compare_QII_Over_SigDigi_C->Scale(1./Compare_QII_Over_SigDigi_C->Integral());
    Compare_QII_Over_SigDigi_FL->SetLineColor(kGreen);
    Compare_QII_Over_SigDigi_FL->SetLineWidth(2);
    Compare_QII_Over_SigDigi_FL->Scale(1./Compare_QII_Over_SigDigi_FL->Integral());
    Compare_QII_Over_SigDigi_FR->SetLineColor(kBlue);
    Compare_QII_Over_SigDigi_FR->SetLineWidth(2);
    Compare_QII_Over_SigDigi_FR->Scale(1./Compare_QII_Over_SigDigi_FR->Integral());
    
    Compare_QII_Minus_SigDigi->SetLineColor(kBlack);
    Compare_QII_Minus_SigDigi->SetLineWidth(2);
    Compare_QII_Minus_SigDigi->SetTitle("QII - Sum_sigdigi");
    Compare_QII_Minus_SigDigi_L->SetLineColor(kBlack);
    Compare_QII_Minus_SigDigi_L->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_L->Scale(1./Compare_QII_Minus_SigDigi_L->Integral());
    Compare_QII_Minus_SigDigi_R->SetLineColor(kRed);
    Compare_QII_Minus_SigDigi_R->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_R->Scale(1./Compare_QII_Minus_SigDigi_R->Integral());
    Compare_QII_Minus_SigDigi_C->SetLineColor(kOrange);
    Compare_QII_Minus_SigDigi_C->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_C->Scale(1./Compare_QII_Minus_SigDigi_C->Integral());
    Compare_QII_Minus_SigDigi_FL->SetLineColor(kGreen);
    Compare_QII_Minus_SigDigi_FL->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_FL->Scale(1./Compare_QII_Minus_SigDigi_FL->Integral());
    Compare_QII_Minus_SigDigi_FR->SetLineColor(kBlue);
    Compare_QII_Minus_SigDigi_FR->SetLineWidth(2);
    Compare_QII_Minus_SigDigi_FR->Scale(1./Compare_QII_Minus_SigDigi_FR->Integral());

    /*TCanvas *c_left_1=new TCanvas("c_left_1","c_left_1",1100,700);
    c_left_1->cd();
    gStyle->SetOptStat(0);
    Sigdigi_left_1->Draw("hist text");
    line[0]->Draw("same");
    l3->Draw();

    TCanvas *c_left_2=new TCanvas("c_left_2","c_left_2",1100,700);
    c_left_2->cd();
    gStyle->SetOptStat(0);
    Sigdigi_left_2->Draw("hist text");
    line[1]->Draw("same");
    l3->Draw();

    TCanvas *c_right_1=new TCanvas("c_right_1","c_right_1",1100,700);
    c_right_1->cd();
    gStyle->SetOptStat(0);
    Sigdigi_right_1->Draw("hist text");
    line[2]->Draw("same");
    l3->Draw();

    TCanvas *c_right_2=new TCanvas("c_right_2","c_right_2",1100,700);
    c_right_2->cd();
    gStyle->SetOptStat(0);
    Sigdigi_right_2->Draw("hist text");
    line[3]->Draw("same");
    l3->Draw();

    TCanvas *c_center_1=new TCanvas("c_center_1","c_center_1",1100,700);
    c_center_1->cd();
    gStyle->SetOptStat(0);
    Sigdigi_center_1->Draw("hist text");
    line[4]->Draw("same");
    l3->Draw();

    TCanvas *c_center_2=new TCanvas("c_center_2","c_center_2",1100,700);
    c_center_2->cd();
    gStyle->SetOptStat(0);
    Sigdigi_center_2->Draw("hist text");
    line[5]->Draw("same");
    l3->Draw();

    TCanvas *c_FullLeft=new TCanvas("c_FullLeft","c_FullLeft",1100,700);
    c_FullLeft->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullLeft->Draw("hist text");
    line[6]->Draw("same");
    l3->Draw();

    TCanvas *c_FullRight=new TCanvas("c_FullRight","c_FullRight",1100,700);
    c_FullRight->cd();
    gStyle->SetOptStat(0);
    Sigdigi_FullRight->Draw("hist text");
    line[7]->Draw("same");
    l3->Draw();*/

    TCanvas *c_QII_Over=new TCanvas("c_QII_Over","c_QII_Over",1100,700);
    c_QII_Over->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Over_SigDigi->Draw("hist");
    //Compare_QII_Over_SigDigi_FR->Draw("hist");
    //Compare_QII_Over_SigDigi_L->Draw("hist same");
    //Compare_QII_Over_SigDigi_R->Draw("hist same");
    //Compare_QII_Over_SigDigi_C->Draw("hist same");
    //Compare_QII_Over_SigDigi_FL->Draw("hist same");
    auto lb= new TLegend(0.6,0.35,0.85,0.55);
    lb->AddEntry(Compare_QII_Over_SigDigi_L,"Left","l");
    lb->AddEntry(Compare_QII_Over_SigDigi_R,"Right","l");
    lb->AddEntry(Compare_QII_Over_SigDigi_C,"Center","l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FL,"Full Left","l");
    lb->AddEntry(Compare_QII_Over_SigDigi_FR,"Full Right","l");
    lb->SetBorderSize(0);
    //lb->Draw();
    
    TCanvas *c_QII_Minus=new TCanvas("c_QII_Minus","c_QII_Minus",1100,700);
    c_QII_Minus->cd();
    gStyle->SetOptStat(0);
    Compare_QII_Minus_SigDigi->Draw("hist");
    //Compare_QII_Minus_SigDigi_FR->Draw("hist");
    //Compare_QII_Minus_SigDigi_L->Draw("hist same");
    //Compare_QII_Minus_SigDigi_R->Draw("hist same");
    //Compare_QII_Minus_SigDigi_C->Draw("hist same");
    //Compare_QII_Minus_SigDigi_FL->Draw("hist same");
    //Compare_QII_Minus_SigDigi_FR->SetAxisRange(0,0.5,"Y");
    auto la= new TLegend(0.6,0.35,0.85,0.55);
    la->AddEntry(Compare_QII_Minus_SigDigi_L,"Left","l");
    la->AddEntry(Compare_QII_Minus_SigDigi_R,"Right","l");
    la->AddEntry(Compare_QII_Minus_SigDigi_C,"Center","l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FL,"Full Left","l");
    la->AddEntry(Compare_QII_Minus_SigDigi_FR,"Full Right","l");
    la->SetBorderSize(0);
    //la->Draw();

}

void ChargeDistribution()
{
    TFile *MuonNoPU = new TFile("ROOT_SVG/MuonPU_checkCorrection_Reco_MuonNoPU.root", "READ");
    TFile *MuonPU = new TFile("ROOT_SVG/MuonPU_checkCorrection_Reco_MuonPU.root", "READ");

    TH1F *Reco_NoPU = (TH1F*)MuonNoPU->Get("ChargeDistributionReco");
    TH1F *Reco_PU = (TH1F*)MuonPU->Get("ChargeDistributionReco");
    TH1F *SigDigi_PU = (TH1F*)MuonPU->Get("ChargeDistributionSigdigi");
    TH1F *SigDigi_PU_zerosupress = (TH1F*)MuonPU->Get("ChargeDistributionSigdigi_zerosupress");

    Reco_NoPU->SetLineColor(kBlack);
    Reco_PU->SetLineColor(kRed);
    SigDigi_PU->SetLineColor(kBlue);
    SigDigi_PU_zerosupress->SetLineColor(kGreen+1);
    Reco_NoPU->GetYaxis()->SetTitle("Entries");
    Reco_NoPU->Scale(1./Reco_NoPU->Integral());
    Reco_PU->Scale(1./Reco_PU->Integral());
    SigDigi_PU->Scale(1./SigDigi_PU->Integral());
    SigDigi_PU_zerosupress->Scale(1./SigDigi_PU_zerosupress->Integral());

    TCanvas *canvas=new TCanvas("canvas","canvas",1100,700);
    canvas->cd();
    gStyle->SetOptStat(0);
    Reco_NoPU->Draw("hist");
    Reco_PU->Draw("hist same");
    SigDigi_PU->Draw("hist same");
    SigDigi_PU_zerosupress->Draw("hist same");
    auto legend= new TLegend(0.6,0.35,0.85,0.55);   
    legend->AddEntry(Reco_PU,Form("Reco PU, MPV=%3.2f",Reco_PU->GetXaxis()->GetBinCenter(Reco_PU->GetMaximumBin())),"l");
    legend->AddEntry(Reco_NoPU,Form("Reco no PU, MPV=%3.2f",Reco_NoPU->GetXaxis()->GetBinCenter(Reco_NoPU->GetMaximumBin())),"l");
    legend->AddEntry(SigDigi_PU,Form("SigDigi, MPV=%3.2f",SigDigi_PU->GetXaxis()->GetBinCenter(SigDigi_PU->GetMaximumBin())),"l");
    legend->AddEntry(SigDigi_PU_zerosupress,Form("SigDigi with zero suppression, MPV=%3.2f",SigDigi_PU_zerosupress->GetXaxis()->GetBinCenter(SigDigi_PU_zerosupress->GetMaximumBin())),"l");
    legend->SetBorderSize(0);
    legend->Draw();

     TH2F *dEdx_sigdigi_VS_recoPU = (TH2F*)MuonPU->Get("dEdx_sigdigi_VS_recoPU");

    // créer un TH1 de projection X du TH2:
    TH1D *projX = dEdx_sigdigi_VS_recoPU->ProjectionX("projX");
    TH1D *projY = dEdx_sigdigi_VS_recoPU->ProjectionY("projY");
    projX->SetLineColor(kBlack);
    projY->SetLineColor(kRed);

    TCanvas *canvas2=new TCanvas("canvas2","canvas2",1100,700);
    canvas2->cd();
    gStyle->SetOptStat(0);
    projX->Draw("hist");
    projY->Draw("hist same");

}

void play_histo_sat()
{
   //histo_Wjets_alltracks();
   //histo_Wjets_reconstructed();
   //histo_MuonPU_1strip();
   ChargeDistribution();
   //histo_MuonPU_2strips();
}
