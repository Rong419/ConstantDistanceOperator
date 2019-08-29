#!/bin/sh
# This shell script is used to extract samples from a tree file (output of MCMC run)
# Start from 101-th state 
# End at 10001-th state
# Sample every 100 states



for file in $(seq 101 100 10001) 
do
   sed -n ${file}p /Users/rzha419/Desktop/ratites.trees >>./tree.txt
done


# delete the previous 20 characters on each line
#sed -i '.bak' 's/^....................//g' tree.txt

# delte the last 5 characters
#sed -i 's/.....$//g' 
