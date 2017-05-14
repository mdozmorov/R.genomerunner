# http://quinlanlab.org/tutorials/bedtools/bedtools.html
# curl -O http://quinlanlab.cs.virginia.edu/cshl2013/make-matrix.py
# ln -s ../bed/*.bed .

#all:	make-matrix
#all:	join-matrix-long
all:	join-matrix-wide

SHELL=/bin/bash

make-matrix:
	parallel "bedtools jaccard -a {1} -b {2} \
	         | awk 'NR>1' \
	         | cut -f 3 \
	         > {1}.{2}.jaccard" \
	         ::: *.bed ::: *.bed

join-matrix-long:
	find . \
    | grep jaccard \
    | xargs grep "" \
    | sed -e s"/\.\///" \
    | perl -pi -e "s/.bed./.bed\t/" \
    | perl -pi -e "s/.jaccard:/\t/" \
    > jaccard.matrix.long.txt

join-matrix-wide:
	awk 'NF==3' jaccard.matrix.long.txt \
	| python make-matrix.py > jaccard.matrix.wide.txt

clean:
	rm *.bed.jaccard
