
import os
import distutils.spawn
old_spawn = distutils.spawn.spawn
def my_spawn(*args, **kwargs):
    print( " ".join(args[0]) ) # <-- this is your command right here 
    old_spawn(*args, **kwargs)
distutils.spawn.spawn = my_spawn
import distutils.ccompiler

c = distutils.ccompiler.new_compiler( verbose=1  )

c.add_include_dir( "." )

sources = ("blobs.c cdiffraction.c check_cpu_auto.c closest.c connectedpixels.c darkflat.c " + \
          "deriv.c localmaxlabel.c sparse_image.c write_check.c").split()
          
c.compile( sources , output_dir="avx2",
        extra_preargs=["/arch:AVX2","/Ox"],  )

objs = [os.path.join("avx2", f.replace(".c",".obj")) for f in sources]

c.create_static_lib( objs, "libcImageD11_avx2", output_dir="." )