#!/bin/bash

#rm -rf ./plots/*.png
#scp -r trcn27@hamilton8.dur.ac.uk:/home/trcn27/mf_emergence/fltrace/plots ./

ffmpeg -y -framerate 20 -i ./fancies/fancy%03d.png -b:v 10M mf.mp4



 
 
