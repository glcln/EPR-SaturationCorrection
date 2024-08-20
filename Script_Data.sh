#!/bin/bash

root -l <<EOC
.L CorrectionOnData.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodFeb2022_CMSSW_10_6_27_probq/SingleMuon/UL2018_RunD/220201_164157/0000/nt_aod_ul2018_*.root")
CorrectionOnData t(chain)
t.Loop()
EOC