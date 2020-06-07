./PrivacyRun file sssp 10000 65 1 1000 1 5 3 privacy
#./PrivacyRun 1k pagerank 10000 20 1 1000 1 3 3 no 
# calculate avg of every 1000 lines
#awk '{sum+=$1} (NR%1000)==0{print sum/1000; sum=0;}' noiserank
