FAST5PATH=$1
NORMALIZEDSAM=$2
OUTPATH=$3

if [ ! -d $OUTPATH ]
then
  mkdir $OUTPATH
fi

READIDS=NORMALIZEDSAM'.ids.txt'
grep -v '^@' $NORMALIZEDSAM | cut -f 1 > $READIDS

fast5_subset --input $FAST5PATH --save_path $OUTPATH --read_id_list $READIDS --batch_size 100 --recursive
