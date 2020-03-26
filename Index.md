


Firts, create a directory and store all the paired-end sequences for each individual in this folder. Also stored the reference genome there 

#### Aligning Sequences

Is important to note that for hundreds of PE sequences the name of the files is usually like teo1_1.fq and teo1_2.fq. However, this scritp considers eliminate the column of the PE direction 
    
     for i in *.fq; do yolo=$(echo "$i" | rev | cut -c 6- | rev | uniq); bwa mem index_tic_geno ${yolo}_1.fq ${yolo}_2.fq -t 20 | samtools view -Sb -@ 20 - > ${yolo}.bam; done &> radseq_10_14.txt &
     
 



