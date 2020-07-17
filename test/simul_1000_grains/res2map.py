
import sys
from ImageD11.grain import write_grain_file, grain, read_grain_file

ubi = sys.argv[1]

grains = read_grain_file( ubi )

res = sys.argv[2]

for line in open(res).readlines():
    if line.find("pos_grains")==0:
        items=line.split()
        j = int(items[0].split("_")[-1])
        grains[j].translation = [float(x)*1000 for x in items[1:]]

write_grain_file(sys.argv[3], grains)
    
