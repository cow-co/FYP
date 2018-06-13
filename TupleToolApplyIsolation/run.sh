#!bin/bash
rm -f nohup.out
nohup gaudirun.py InclusiveMuMu_MC.py FilesKstMuMu_MC.py&
disown %1