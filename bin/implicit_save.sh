#!/bin/bash

cp $1 buf
sed -i -e "s/\(.*\)\!\(<\|\!\)\(.*\)$/\1/" buf
foo=$(grep '&$' buf)

while [ -n "$foo" ]; do
   sed -i -e 'N;s/&\n/ /;P;D;' buf
   foo=$(grep '&$' buf)
done
ims=$(sed -n "/\(real\|integer*\).*::.*=.*[0-9]/p" buf | grep -v "save\|parameter")
if [[ -n $ims ]]; then
   echo -e "\033[91mQA: Implicit save detected in\033[0m $1"
   # once again for pretty output
   sed -n "/\(real\|integer*\).*::.*=.*[0-9]/p" buf | grep -v "save\|parameter"
   rm -rf buf
   exit 0
fi
rm -rf buf
exit 1
