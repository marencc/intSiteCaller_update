#!/bin/bash
gfServer start localhost $2 -tileSize=11 -repMatch=112312 -maxDnaHits=20 -canStop $3 &

#give enough time for gfServer to fully load entire index (can take awhile if many gfServer instances are using the same index)
while ! ( gfClient -t=dna -q=dna -minIdentity=85 -minScore=27 -dots=1000 -out=psl -nohead localhost $2 . $1 $1.psl )
do
  sleep 60
done

gfServer stop localhost $2
gzip $1.psl
