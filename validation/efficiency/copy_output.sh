#!/bin/sh

for data in {anolis,RSV2}
do  
 for param in {1..2}
 do 
    for sim in {1..100}
    do
    #scp login.mahuika.nesi.org.nz:/nesi/project/nesi00390/rong/operator/others/output/${data}_${param}_${sim}.log ~/Desktop/efficiency/${data}/logs/
    #scp login.mahuika.nesi.org.nz:/nesi/project/nesi00390/rong/operator/others/output/${data}_${param}_${sim}.trees ~/Desktop/efficiency/${data}/trees/
    scp login.mahuika.nesi.org.nz:/nesi/project/nesi00390/rong/operator/others/summarytree/s_${data}_${param}_${sim}.trees ~/Desktop/efficiency/${data}/summaries/
    #scp login.mahuika.nesi.org.nz:/nesi/project/nesi00390/rong/operator/others/output/output_${data}_${param}_${sim}.txt ~/Desktop/efficiency/${data}/outputs/
    #scp login.mahuika.nesi.org.nz:/nesi/project/nesi00390/rong/operator/others/output/error_${data}_${param}_${sim}.txt ~/Desktop/efficiency/${data}/times/
    done
  done
done
