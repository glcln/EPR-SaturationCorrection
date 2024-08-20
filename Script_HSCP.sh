#!/bin/bash

root -l <<EOC
.L sat_observation_HSCP.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/safe1/cms/gcoulon/Stage2023/ROOT_Files/nt_HSCPgluino_2TeV.root")
sat_observation_HSCP t(chain)
t.Loop()
EOC
