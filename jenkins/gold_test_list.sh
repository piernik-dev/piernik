#!/bin/bash

CONF_DIR=./jenkins/gold_configs

# Set up the workspace
OUT_DIR=./jenkins/workspace/
[ ! -d $OUT_DIR ] && mkdir -p $OUT_DIR
[ ! -d $OUT_DIR ] && exit 1

# Execute all tests in serial or parallel mode
[ -z ${SERIAL+x} ] && SERIAL=0
which parallel > /dev/null 2>&1 || SERIAL=1

if [ $SERIAL -ne 0 ] ; then
    echo "serial-gold"
    for i in ${CONF_DIR}/*.config; do
	echo "  $i"
	eval "./jenkins/gold_test.sh $i > ${OUT_DIR}$(basename $i)'_gold_stdout'"
    done
else
    echo "parallel-gold"
    parallel --load 70% --delay 10 eval "./jenkins/gold_test.sh {} > ${OUT_DIR}{/}'_gold_stdout'" ::: ${CONF_DIR}/*.config
fi

# Create the .csv file with norms for Jenkins
for j in gold riem ; do
    case $j in
         ("gold") jj="gold" ;;
         ("riem") jj="Riemann" ;;
         (*) jj="___" ;;
    esac
    for i in ${CONF_DIR}/*.config ; do
        eval $( grep PROBLEM_NAME $i )
	PROBLEM_NAME=${PROBLEM_NAME//\//___}

	LOG=${OUT_DIR}${PROBLEM_NAME}/${j}_log
	if [ -e $LOG ] ; then
           echo ${PROBLEM_NAME}_${jj}" "$( tail -n 1 $LOG | awk '{print $NF}' )
        else
	   echo ${PROBLEM_NAME}_${jj}" "-0.1
        fi
    done |\
        awk '{\
		    a[$1]=$2;\
		} END {\
		    for (i in a) printf("%s,",i);\
		    print "";\
		    for (i in a) printf("%s,",a[i]);\
		    print "";\
		}' |\
        sed 's/,$//' > ${OUT_DIR}all_${jj}.csv
done

# Print the results
for i in ${OUT_DIR}*/gold_log ; do
    gdist=$( tail -n 1 $i | awk '{print $NF}' )
    printf "%-50s = %s\n" "[GOLD] Total difference for ${i/_gold_log/}" $gdist
    [ ${gdist} != "0" ] && sed -n '/^Difference of /s/^Diff/    Diff/p' $i
done
for i in ${OUT_DIR}*/riem_log ; do
    printf "%-50s = %s\n" "[Riemann] Total difference for ${i/_riem_log/}" $( tail -n 1 $i | awk '{print $NF}' )
done

# Fail if any gold distance is not 0.
for i in ${OUT_DIR}*/gold_log ; do
    [ $( ( grep "^Total difference between" $i || echo 1 ) | awk '{print $NF}' ) == 0 ] || exit 1
done
