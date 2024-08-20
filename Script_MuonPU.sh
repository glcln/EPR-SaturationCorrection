#!/bin/bash

##### MuonPU one saturated strip

root -l <<EOC
.L sat_observation_MuonPU.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodMarch2023_CMSSW_10_6_30/hadd_SingleMu_pt50_200_eta2p1_wPu/nt_SingleMu_pt50_200_eta2p1_wPu.root")
sat_observation_MuonPU t(chain)
t.Loop()
EOC

##### MuonPU two consecutives saturated strips

#root -l <<EOC
#.L sat_observation_MuonPU_twoStrips.C
#TChain *chain = new TChain("stage/ttree")
#chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodMarch2023_CMSSW_10_6_30/hadd_SingleMu_pt50_200_eta2p1_wPu/nt_SingleMu_pt50_200_eta2p1_wPu.root")
#sat_observation_MuonPU_twoStrips t(chain)
#t.Loop()
#EOC
