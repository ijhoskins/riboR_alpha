





# Create WT_1 with rnaseq and metadata
ribog create --name WT_1 --alignmentfile alignments_2.bed --reference appris_human \
	     --lengths transcript_lengths.tsv --annotation annotation.bed --metageneradius 2 \
	     --leftspan 3 --rightspan 2 --lengthmin 2 --lengthmax 5 \
	     --expmeta wt.yaml --ribometa ribo_meta.yaml \
	     wt_1.ribo

ribog rnaseq set --name WT_1 --counts rnaseq_1.tsv --sep $'\t'  wt_1.ribo

# Create WT_2 without metadata or rnaseq
ribog create --name WT_2 --alignmentfile alignments_2.bed --reference appris_human \
	     --lengths transcript_lengths.tsv --annotation annotation.bed --metageneradius 2 \
	     --leftspan 3 --rightspan 2 --lengthmin 2 --lengthmax 5\
	     --expmeta wt.yaml --ribometa ribo_meta.yaml \
	     wt_2.ribo


# Create WT_3 without coverage data
ribog create --name WT_3 --alignmentfile alignments_2.bed --reference appris_human \
	     --lengths transcript_lengths.tsv --annotation annotation.bed --metageneradius 2 \
	     --leftspan 3 --rightspan 2 --lengthmin 2 --lengthmax 5\
	     --expmeta wt.yaml --ribometa ribo_meta.yaml \
             --nocoverage \
	     wt_3.ribo

# Create Hela sample
# with rnaseq and metadata
ribog create --name Hela_1 --alignmentfile alignments_1.bed --reference appris_human \
	     --lengths transcript_lengths.tsv --annotation annotation.bed --metageneradius 2 \
	     --leftspan 3 --rightspan 2 --lengthmin 2 --lengthmax 5\
	     --expmeta hela_meta.yaml --ribometa ribo_meta.yaml \
	     hela_1.ribo

ribog rnaseq set --name Hela_1 --counts rnaseq_1.tsv --sep $'\t'  hela_1.ribo

# Create Hela Sample witohout metadata with rnaseq
ribog create --name Hela_2 --alignmentfile alignments_1.bed --reference appris_human \
	     --lengths transcript_lengths.tsv --annotation annotation.bed --metageneradius 2 \
	     --leftspan 3 --rightspan 2 --lengthmin 2 --lengthmax 5\
	     --ribometa ribo_meta.yaml \
	     hela_2.ribo

ribog rnaseq set --name Hela_2 --counts rnaseq_2.tsv --sep $'\t'  hela_2.ribo

# Finally merge the ribo files

ribog merge sample.ribo hela_1.ribo hela_2.ribo wt_1.ribo wt_2.ribo wt_3.ribo

