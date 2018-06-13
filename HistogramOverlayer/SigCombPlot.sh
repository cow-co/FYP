#!/bin/bash
echo Starting
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_DELTA_VERTEX_CHI2__0 0 50 100 Low/Comb_B0_DELTA_VERTEX_CHI2__0
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_MIN_CHI2_ 0 0.5 100 Low/Comb_B0_MIN_CHI2_
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_MIN_DV_CHI2_ 0 5 100 Low/Comb_B0_MIN_DV_CHI2_
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_MUON_VTX_ 0 1 100 Low/Comb_B0_MUON_VTX_
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_R_MX1_TO_MX2 0 10 100 Low/Comb_B0_R_MX1_TO_MX2
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_R_MM_TO_MpX 0 10 100 Low/Comb_B0_R_MM_TO_MpX
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_R_MM_TO_MmX 0 10 100 Low/Comb_B0_R_MM_TO_MmX
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_SUM_R_MuMu_TO_MuX 0 10 100 Low/Comb_B0_SUM_R_MuMu_TO_MuX
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_DEL_R_MuX_MuMu 0 10 100 Low/Comb_B0_DEL_R_MuX_MuMu
./StandaloneOverlayer SignalLowQ2.root CombBGLow.root B0_R_MuX_MuX 0 10 100 Low/Comb_B0_R_MuX_MuX
echo Done