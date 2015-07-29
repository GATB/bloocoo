bwa index ref.fa
arg="-k 10 -T 20"
bwa mem $arg  ref.fa read_uncor.fa  > read_uncor.sam
bwa mem $arg  ref.fa read_cor.fa  > read_cor.sam

python ../compare_two_sam_files.py read_cor.sam read_uncor.sam ref.fa
