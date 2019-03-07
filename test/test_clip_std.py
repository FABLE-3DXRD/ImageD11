import numpy, ImageD11.cImageD11
numpy.random.seed(42)
im = numpy.random.random((2048,2048)).astype(numpy.float32)
print ImageD11.cImageD11.clip_std(im.ravel(),1)
print im.mean(),im.std(),im.var()
ImageD11.cImageD11.clip_std(im.ravel(),1)
