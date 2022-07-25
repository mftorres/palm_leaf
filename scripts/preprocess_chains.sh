#!/bin/bash -l

FILE=$1
FIELDS=$2

echo $FILE

# if hypotheses group hp1 with FIELDS number of fields plus intercept and index
cut -f 1-$FIELDS $FILE > "${FILE}".out

# remove quotes
sed -i 's/"//g' ${FILE}.out

# remove brackets
sed -i 's/[\(\)]//g' ${FILE}.out

# dd index header
sed -i 's/Int/index\tInt/g' ${FILE}.out
sed -i -r 's/\tanimal.+//g' ${FILE}.out

# could be more elegant but I'm on a rush.
# ran in node interactively
# for i in Shape*hp1*-1.txt; do srun -n 1 prepare_chains.sh $i 7; done;
# for i in Shape*hp1*-2.txt; do srun -n 1 prepare_chains.sh $i 7; done;
# for i in Shape*hp2*-1.txt; do srun -n 1 prepare_chains.sh $i 8; done;
# for i in Shape*hp2*-2.txt; do srun -n 1 prepare_chains.sh $i 8; done;
