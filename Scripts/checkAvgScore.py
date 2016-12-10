import os,sys
import numpy as np

scores = []
with open(sys.argv[1]) as fh:
    for line in fh:
        info = line.rstrip('\n').split('\t')
        chrom,start,end,score = info[0],info[1],info[2],info[3]
        if chrom == '15':
            if int(start) >= 23718546 and int(start) <= 23753269:
                scores.append(float(score))
            else:
                pass

print np.mean(scores)
print len(scores)
