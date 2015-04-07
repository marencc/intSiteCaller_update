gfServer start localhost $2 -tileSize=11 -repMatch=112312 -maxDnaHits=20 -canStop /home/aubreyba/genomeIndices/hg18.2bit &
sleep 600
gfClient -t=dna -q=dna -minIdentity=85 -minScore=27 -dots=1000 -out=psl -nohead localhost $2 / $1 $1.psl
gfServer stop localhost $2
gzip $1.psl
