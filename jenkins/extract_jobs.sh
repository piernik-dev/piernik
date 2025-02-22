#!/bin/bash
for i in ${JENKINS_HOME}/jobs/* ; do
	awk 'BEGIN {p=0}; {if ($0 ~ "<command>") p=1; if (p) print; if ($0 ~ "</command>") p=0}' "$i"/config.xml |\
	sed 's/.*command>//;s/&amp;/\&/g;s/&quot;/"/g;s/&gt;/>/g;s/&apos;/'\''/g' > "$( basename "$i").sh"
done
