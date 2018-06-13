#!bin/bash
rm -f nohup.out
nohup gaudirun.py BGLowQ2.py BackgroundFiles.py&
disown %1