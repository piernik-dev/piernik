#!/bin/tcsh

setenv IDL_PATH "/scrh/kowalik/piernik0.1/pro_gdl:<IDL_DEFAULT>"

echo " "
echo "WORKING DIRECTORY:"
pwd

if (! -f "./plot_tsl.var") then
  cp ~/piernik0.1/pro/plot_tsl.var .
endif

gdl ~/piernik0.1/pro/plot_tsl.idl

echo " "
