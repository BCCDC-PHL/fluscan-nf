"""BLAST tabular output shared by the BLAST-parsing entry points."""

# Must match the -outfmt the blastn and select_db processes pass
BLAST_COLUMNS = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score'.split()
