#!/bin/bash

root -l <<EOC
.L MuonPU_checkCorrection.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodMarch2023_CMSSW_10_6_30/hadd_SingleMu_pt50_200_eta2p1_noPu/nt_SingleMu_pt50_200_eta2p1_noPu.root")
MuonPU_checkCorrection t(chain)
t.Loop()
EOC
