#!/usr/bin/awk -f

BEGIN {
    l=0;
} {
    n=index($0, "real");
    if (n>l) l=n;
    a[NR]=$0;
} END {
    form=sprintf("%%-%ds %%s %%6.1f %%6.1f %%6.1f %%08s\n", l-1);
    for (i=1;i<=NR;i++) {
	$0=a[i];
	printf(form, substr($0, 0, index($0, "real")-1), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $(NF));
    }
}
