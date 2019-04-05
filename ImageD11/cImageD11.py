

# check for avx2 and import cImageD11 module
import platform, sys
if platform.system() == "Linux":
    if " avx2 " in open("/proc/cpuinfo").read():
        from cImageD11_avx2 import *
        cpu = 'avx2'
    else:
        from cImageD11_sse2 import *
        cpu = 'sse2'

elif platform.system() == "Windows":
    if sys.version_info[:2] > (3,4):
        from cImageD11_avx2 import *
        cpu = 'avx2'
    else:
        from cImageD11_sse2 import *
        cpu = 'sse2'
else:
    
    pass
