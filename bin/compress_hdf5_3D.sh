#!/bin/bash

# AJG: This script may compress/recompress/decompress HDF5 plotfiles produced by Piernik
# Restart files (with four-dimensional arrays) are not supported yet.
# Note that if you use SZIP method, you must have szip library accesible in your environment
# not only at compression time, but also for visualization and other tools.
# Warning: Sometimes it may be problematic to perform a restart from compressed restart file.
# One may expect plotfiles to shrink by 30-50% (compression efficiency 150-200%)

#constants
space="                                                                                     "
CSZ="-f SZIP=32,NN"
CGZ="-f SHUF -f GZIP=9"

#defaults
CMETH=$CGZ
autofind=0
clev=110
decompress=0
verbose=0

#usage
if [ $# -lt 1 ] ; then
  echo "Usage $0 [options] Piernik_file.h5 [ ... ]"
  echo "Options:"
  echo "  -gz : compression method $CGZ (default)"
  echo "  -sh : compression method $CSZ"
  echo "  -d|-decompress"
  echo "  -a|-auto : find files automatically"
  echo "  -clev <number> : try to recompress if files found have compression efficiency less than <number>%"
  echo "  -v|-verbose : verbose"
  exit 1
fi

#required tools
for i in h5ls h5repack h5dump ; do
  which $i > /dev/null || exit 2
done

#command line interpreter
flags=1
while [ $flags == 1 ] ; do
  case $1 in
    ("-gz") CMETH=$CGZ ; shift ;;
    ("-sz") CMETH=$CSZ ; shift ;;
    ("-a"|"-auto") autofind=1; shift ;;
    ("-d"|"-decompress") decompress=1; shift ;;
    ("-clev")
      shift
      if [ $# -ge 1 ] ; then
        [ $1 -ge 100 ] && clev=$1 || echo "$1% threshold makes no sense, leaving default ($clev%)"
        shift
      fi ;;
    ("-v"|"-verbose") verbose=1; shift ;;
    (*) flags=0 ;;
  esac
done
[ $decompress == 1 ] && CMETH="-f NONE"

#shortcuts
hreadform()
{
  awk 'BEGIN {\
    n='$1';\
    if (n<0) n=-n;\
    if (n<1024) print n" B";\
    else if (n<1024**2) print int(10.*n/1024)/10." kB";\
    else if (n<1024**3) print int(10.*n/1024**2)/10." MB";\
    else if (n<1024**4) print int(10.*n/1024**3)/10." GB";\
    else print int(10.*n/1024.**4)/10." TB";\
  }'
}

validhdf() { h5dump -n "$1" 2> /dev/null 1>&2 ; echo $?; }

cratio()
{
  h5ls -v "$1" |\
    grep util |\
    sed 's/ *Storage: *//' |\
    sort -nr |\
    awk '{if (NR==1) sz=$1; if ($1==sz) { sl+=$1; sa+=$4;}} END {print int(100.*sl/sa);}'
}

#list of files to operate on
if [ $autofind -ne 0 ] ; then
  list=""
  for i in $( find . -name "*[0-9][0-9][0-9][0-9].h5" ) ; do
    if [ $( validhdf "$i" ) == 0 -a ! -h "$i" ] ; then
      cr=$( cratio "$i" )
      [ $cr -lt $clev ] && list="$list $i"
      ns=$(( ${#oi} - ${#i} ))
      [ $ns -lt 0 ] && ns=0
      echo -ne "checking: $i, estimated compression efficiency $cr% ${space:0:$ns}\r"
      oi=$i
    fi
  done
  echo
else
  list=$*
fi

#performance counters
TRIED=0
MODED=0
INVAL=0
SAVED=0
ITOT=0
OTOT=0

#main loop
for i in $list ; do
  SKIP=0
  TRIED=$(( $TRIED + 1 ))
  if [ $( validhdf "$i" ) -ne 0 ] ; then
    echo "Cannot read $i as an HDF5 file, skipping"
    INVAL=$(( $INVAL + 1 ))
  elif [ ! -f "$i" -o -h "$i" ] ; then
    echo "$i is a link or not a regular file, skipping"
    INVAL=$(( $INVAL + 1 ))
  else
    [ $autofind == 0 -a  $decompress == 0 ] && cr=$( cratio "$i" ) || cr=$clev
    if [ $cr -gt $clev ] ; then
      echo "$i is already compressed with estimated efficiency $cr%, skipping"
    else
      if [  $decompress == 0 ] ; then
#Keep the chunksize below 1MB to prevent I/O performance degradation
        FDIM=$( h5ls "$i" | sed -n 's/.*{\([0-9 ,]*\)}/\1/;/, .*, .*/s/,//gp' | uniq )
        BSIZE=$( echo $FDIM | awk '{print 4*$2*$3*$1}' )
#BEWARE: for restarts the above should be 8*$2*$3*$1
        if [ $BSIZE -gt 0 ] ; then
          CHUNK=$( echo $FDIM | awk 'function min(x,y){return x<y?x:y} {print min($1,64)"x"min($2,32)"x"min($3,32)}' )
          LAYOUT="-l CHUNK=$CHUNK"
        else
          SKIP=1
          INVAL=$(( $INVAL + 1 ))
          echo "Warning: Cannot determine layout of file $i (FDIM=$FDIM, BSIZE=$BSIZE), skipping"
        fi
      else
        LAYOUT=""
      fi
      if [ $SKIP == 0 ] ; then
        ofile=$( mktemp hdf5_XXXXXX )
        [ $verbose -ne 0 ] && echo h5repack -i "$i" -o $ofile $LAYOUT $CMETH
        h5repack -i "$i" -o $ofile $LAYOUT $CMETH || exit 3
        ISIZE=$( du -B 1 "$i" | awk '{print $1}' )
        OSIZE=$( du -B 1 $ofile | awk '{print $1}' )
        touch -r "$i" $ofile
        chmod --reference="$i" $ofile
        chown --reference="$i" $ofile
        DIFFSIZE=$(( $ISIZE - $OSIZE ))
        DIFFSTR=$( hreadform $DIFFSIZE )
        if [ $DIFFSIZE -gt 0 -o $decompress == 1 ] ; then
          [ $DIFFSIZE -lt 0 ] && grshr="grows" || grshr="shrinks"
          echo "$i ${grshr} by ${DIFFSTR}"
          mv $ofile "$i"
          SAVED=$(( $SAVED + $DIFFSIZE ))
          MODED=$(( $MODED + 1 ))
          OTOT=$(( $OTOT + $OSIZE ))
        else
          if [ $DIFFSIZE == 0 ] ; then
            echo "$i would not change its size, operation cancelled"
          else
            echo "$i would grow by ${DIFFSTR}, operation cancelled"
          fi
          \rm -f $ofile
          OTOT=$(( $OTOT + $ISIZE ))
        fi
        ITOT=$(( $ITOT + $ISIZE ))
      fi
    fi
  fi
done

echo "Summary:"
echo "  $MODED out of $TRIED files were modified ($INVAL found invalid)"
echo "  Initial total size was "$( hreadform $ITOT )", now is "$( hreadform $OTOT )
[ $SAVED -ge 0 ] && echo -n "  Saved " || echo -n "  Lost "
echo $( hreadform $SAVED )", compression ratio is "$( awk 'BEGIN { o='$OTOT'; if (o > 0) printf("%.2f%%",100*'$ITOT'/o); else print "unspecified" }' )
