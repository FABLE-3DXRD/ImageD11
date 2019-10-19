

"""
We started using some intrinsics in the C code to make it run faster
when processors support specific instructions.  This makes it less
portable and harder to build. These notes explain how things could 
work. 

If a CPU has some feature (sse2, avx2, fma, etc) then we will want
to use it in specific places that we are aware of. We might as well
also let the compiler use those instructions everywhere else where
we have not written specific things too. That means:

- compile everything as avx2 -> avx2 version
- compile everything as sse2 -> sse2 version
- compile everything as ansi -> ansi version

So we are going to build a series of extensions and use #ifdef inside
them.

Which one we get is determined by the cImageD11.py module.
"""

import os, sys, platform
import setuptools
import distutils.ccompiler

need_build = True

sources = ("blobs.c cdiffraction.c closest.c connectedpixels.c"+
 " darkflat.c localmaxlabel.c sparse_image.c splat.c").split()

plat = platform.system()
bits = platform.architecture()[0]
vers = "%d.%d"%(sys.version_info[:2])
tmpdir = "%s_%s_%s"%(plat, bits, vers)
avx2libname = "cImageD11_"+tmpdir+"_avx"
sse2libname = "cImageD11_"+tmpdir+"_sse2"
# ansi ?

compiler = None
for a in sys.argv:
    if "mingw32" in a:
        compiler = "mingw32"



if plat == "Linux" or compiler == "mingw32":
    arg=["-O2", "-fopenmp", "-fPIC", "-std=c99" ]
    sse2arg = arg + ["-msse4.2"]
    avx2arg = arg + ["-mavx"]
    # link args
    lsse2arg = arg + ["-msse4.2"]
    lavx2arg = arg + ["-mavx"]
elif plat == "Windows":
    arg=["/Ox", "/openmp" ]
    # the /arch switches are ignored by the older MSVC compilers
    sse2arg = arg + ["/arch:SSE2",]
    avx2arg = arg + ["/arch:AVX",]
    lsse2arg = []
    lavx2arg = []
else:
    avx2arg = sse2arg = arg = [ ]

def run_cc( cc, plat, bits, vers, name, flags, libname ):
    objs = cc.compile( sources , 
                       output_dir=libname.replace("cImageD11_",""),
                       extra_preargs = flags )
    ok = cc.create_static_lib( objs, libname, output_dir="." )
    return libname

def write_docs( inp, outf ):
    """ One single block of !DOC per item, first word is key
    """
    with open(inp , "r") as pyf:
        fname = None
        docs = {}
        for line in pyf.readlines():
            if line.startswith("!DOC"):
                if fname is not None:
                    docs[fname] += line[5:] 
                else:
                    words = line.split()
                    fname = words[1]
                    docs[fname] = " ".join(words[2:]) + "\n"
            else:
                fname = None
    with open(outf, "w") as docf:
        keys = list(docs.keys())
        keys.sort()
        for fname in keys:
            docf.write( '%s = """%s"""\n'%(fname, docs[fname]))
        docf.write("__all__ = [" + ",".join(['"%s"'%(k) for k in keys ])+  "]")
        

def make_pyf( inp, name ):
    out = open(name+".pyf", "w")
    out.write("python module %s\n"%(name))
    pyf = open(inp , "r").read() 
    out.write( pyf )
    out.write("end python module %s\n"%(name))
    out.close()

def docs():
    write_docs( "cImageD11_interface.pyf", os.path.join("..","ImageD11src", "cImageD11_docstrings.py"))


def main():
    cc = distutils.ccompiler.new_compiler( verbose=1 , compiler=compiler )
    cc.add_include_dir( "." )
    docs()
    make_pyf( "cImageD11_interface.pyf", "cImageD11_sse2")
    make_pyf( "cImageD11_interface.pyf", "cImageD11_avx")
    sse2lib = run_cc(cc, plat, bits, vers, "sse2", sse2arg, sse2libname )
    avx2lib = run_cc(cc, plat, bits, vers, "avx", avx2arg, avx2libname )
    return 0


if __name__=="__main__":
    main()
