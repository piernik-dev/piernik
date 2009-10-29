#!/bin/bash
# copyright A. Gawryszczak 2009

url=$( svn info | grep URL | awk '{print $NF}' )
root=$( LC_ALL=C LANG=C LANGUAGE=C svn info | grep "Repository Root: " | awk '{print $NF}' )
relpath=$( echo $url | sed "s^$root/^^" )
last=$( basename $relpath )

#echo $url $root $relpath $last

cmdb='svn ci -m '"'"'[source:'"$relpath"' '"$last"']: '
cmde="'"' .'

if [ $# -lt 1 ] ; then
	echo "Usage: $0 comment"
	echo ${cmdb}${cmde}
	exit 1
fi

#eval "${cmdb}$*${cmde}"

echo `${cmdb} ${cmde}`
