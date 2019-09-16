#!/bin/sh
for data in {anolis,Shankarappa,RSV2}
do  
  for param in {1..2}
  do 
    for sim in {1..3}
    do
    sed "s/FILE/${param}_${sim}/g" ./${data}_${param}.xml > ./xml/${data}_${param}_${sim}.xml
    echo "use ${data}_${param}.xml to write ./${data}_${param}_${sim}.xml"
    done
  done
done
