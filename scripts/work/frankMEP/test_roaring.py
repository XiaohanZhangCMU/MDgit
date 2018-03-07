import os
import sys
sys.path.append('./PyRoaringBitMap/')

from pyroaring import BitMap

bm1 = BitMap()
bm1.add(3)
bm1.add(18)
bm2 = BitMap([3, 27, 42])
print("bm1       = %s" % bm1)
print("bm2       = %s" % bm2)
print("bm1 & bm2 = %s" % (bm1&bm2))
print("bm1 | bm2 = %s" % (bm1|bm2)) 
