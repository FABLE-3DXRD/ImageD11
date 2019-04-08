
import os, sys, platform

# patch for debugging
import distutils.spawn
old_spawn = distutils.spawn.spawn
def my_spawn(*args, **kwargs):
    print( " ".join(args[0]) ) # <-- this is your command right here 
    old_spawn(*args, **kwargs)
distutils.spawn.spawn = my_spawn

import distutils.ccompiler

sources = ("blobs.c cdiffraction.c check_cpu_auto.c closest.c connectedpixels.c "+\
    "darkflat.c localmaxlabel.c sparse_image.c").split()

plat = platform.system()
bits = platform.architecture()[0]
vers = "%d.%d"%(sys.version_info[:2])
tmpdir = "%s_%s_%s"%(plat, bits, vers)
avx2libname = "cImageD11_"+tmpdir+"_avx2"
sse2libname = "cImageD11_"+tmpdir+"_sse2"

if plat == "Linux":
    arg=["-O2", "-fopenmp", "-fPIC" ]
    sse2arg = arg + ["-msse2"]
    avx2arg = arg + ["-mavx2"]
elif plat == "Windows":
    arg=["/Ox", "/Openmp" ]
    sse2arg = arg + ["/arch:SSE2","/Ox"]
    avx2arg = arg + ["/arch:AVX2","/Ox"]
else:
    avx2arg = sse2arg = arg = [ ]

def run_cc( cc, plat, bits, vers, name, flags, libname ):
    cc.compile( sources , output_dir=tmpdir, extra_preargs = flags )
    objs = [os.path.join(tmpdir, f.replace(".c",cc.obj_extension)) for f in sources]    
    cc.create_static_lib( objs, libname, output_dir="." )
    return libname

if __name__=="__main__":
    cc = distutils.ccompiler.new_compiler( verbose=1  )
    cc.add_include_dir( "." )
    sse2lib = run_cc(cc, plat, bits, vers, "sse2", sse2arg, sse2libname )
    avx2lib = run_cc(cc, plat, bits, vers, "avx2", avx2arg, avx2libname )

