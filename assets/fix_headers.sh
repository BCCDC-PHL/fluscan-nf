for i in *fasta ; do echo $i ; sed -i -E 's/>([^\|]+).+\|([A-Z0-9]+)$/>\1_\2/g' $i ; done

for i in *fasta ; do echo $i ; sed -i -E 's/_M$/_MP/g' $i ; done