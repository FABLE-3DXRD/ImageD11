

def filename_series( stem, first, last, fmt="%s%04d.edf"):
    i = first
    while i <= last:
        yield fmt%(stem,i)
        i += 1

def buffered( names, bufsize ):
    buffer = []
    # init
    while len( buffer ) < bufsize:
        name = next(names)
        buffer.append( name )
#        print("buffer", name )
    for i in range(bufsize//2+1):
#        print("from first block")
        yield buffer[i]
    nread = bufsize
    i = i + 1
    for name in names:
        buffer[nread%bufsize] = name
#        print("in main loop")
        yield buffer[i%bufsize]
        nread += 1
        i += 1
    while i < nread:
#        print("last part")
        yield buffer[i%bufsize]
        i += 1


for name in buffered( filename_series( "data", 0, 10 ), 7 ):
    print (name)

    
