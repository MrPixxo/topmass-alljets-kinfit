#!/usr/bin/bash
export LD_LIBRARY_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/lib/Linux514_x86_64_gcc485/6.24.06/:$LD_LIBRARY_PATH
#export LD_PRELOAD=/usr/lib64/libtbb.so.2:/usr/lib64/libfreetype.so.6:$LD_PRELOAD
export CPLUS_INCLUDE_PATH=/usr/include/:/usr/include/root/:cvmfs/cms.cern.ch/slc7_amd64_gcc530/lcg/root/6.06.00-ikhhed/include/:/cvmfs/cms.cern.ch/slc7_amd64_gcc530/cms/cmssw/CMSSW_8_0_27/src/:$(dirname "$(readlink -f "${BASH_SOURCE[0]}")"):/usr/include/:/usr/include/root/:/usr/include/c++/11/:$CPLUS_INCLUDE_PATH