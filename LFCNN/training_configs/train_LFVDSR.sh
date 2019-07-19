#!/bin/bash
set -e

LOG=/PATH/TO/LOGS/log.log
LFSRNet_solver=LFVDSR_solver.prototxt

TOOLS=/userhome/caffe/.build_release/tools
#TOOLS=/opt/caffe/build/tools
$TOOLS/caffe train --solver=$LFSRNet_solver -gpu 0 2>&1 | tee $LOG $@