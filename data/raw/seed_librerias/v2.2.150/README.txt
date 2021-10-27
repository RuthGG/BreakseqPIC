bplib.fa was created from bplib.300.fa (library v 2.2) with:

cat bplib.300.fa | sed 's/^[atcgATCG]\{75\}//g' | sed 's/[atcgATCG]\{75\}$//g' > bplib.fa
