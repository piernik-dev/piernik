#!/bin/bash

mode="l"
if [ $# -ge 1 ] ; then
	[ $1 == "d" ] && mode="d"
fi

DEST="../../trunk"

makelist() {
	[ $# == 1 ] && LC_ALL=C diff --exclude=.svn -qrI ' $Id:' "$1" "${DEST}/$1" | awk '/Files.*differ/ {print $2}' | grep -v "~"
}

for i in $( makelist src ) $( makelist problems) ; do
	case $mode in
		("l") echo $i ;;
		("d") xemacs -geometry 190x65 --eval "(ediff-files \"$i\" \"${DEST}/$i\")" ;;
		(*) echo "Don't know what to do with $i"
	esac
done

