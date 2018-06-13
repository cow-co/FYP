#!bin/bash
rm -f nohup.out
nohup gaudirun.py LowQ2_MC.py FilesKstMuMu_MC.py&
disown %1