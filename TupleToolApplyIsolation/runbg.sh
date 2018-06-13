#!bin/bash
rm -f nohup.out
nohup gaudirun.py Background_MC.py BackgroundFiles.py&
disown %1