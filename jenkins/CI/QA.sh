date
git log --graph --branches --decorate --oneline -30
echo y | ( bin/qa.py $( find src problems -name "*.F90" ) || echo "Error : QA crashed" ) | sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" ; echo
[ $( git diff | wc -l ) -gt 0 ] && echo "Warning : some QA activity detected"
git diff

