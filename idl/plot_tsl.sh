#!/bin/tcsh

setenv IDL_PATH "~/piernik0.1/pro:<IDL_DEFAULT>"

echo " "
echo "WORKING DIRECTORY:"
pwd

if (! -f "./plot_tsl.var") then
  cp ~/piernik0.1/pro/plot_tsl.var .
endif

idl ~/piernik0.1/pro/plot_tsl.idl

echo " "
