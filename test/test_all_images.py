

import glob, os, time, fabio.openimage, gzip, bz2

times = {}
images = []

for fname in glob.glob(os.path.join("testimages","*")):
    if fname.find("header_only")==-1:
        images.append(fname)

images.sort()


print "I/O 1  : Time to read the image"
print "I/O 2  : Time to read the image (repeat"
print "Fabio  : Time for fabio to read the image"
print "Shell  : Time for shell to do decompression"
print "Python : Time for python to do decompression\n"

print "I/O 1  I/O 2  Fabio  Shell  Python"
for im in images:
    # Network/disk io time first
    start = time.clock()
    the_file = open(im,"rb").read()
    times[im] =  [ time.clock()-start ]
    start = time.clock()
    # Network/disk should be cached
    the_file = open(im,"rb").read()
    times[im].append( time.clock() - start )
    start = time.clock()
    try:
        fim = fabio.openimage.openimage(im)
    except KeyboardInterrupt:
        raise
    except:
        print "Problem with",im
        continue
        # raise
    times[im].append( time.clock() - start )
    nt = 3 ; ns = 2
    # Now check for a fabio slowdown effect    
    if im[-3:] == '.gz':
        start = time.clock()
        os.system("gzip -cd %s | wc -c"%(im))
        times[im].append( time.clock()-start )  
        nt += 1; ns -= 1
        start = time.clock()
        the_file = gzip.GzipFile(im).read()
        times[im].append( time.clock()-start )  
        nt += 1; ns -= 1
    if im[-4:] == '.bz2':
        start = time.clock()
        os.system("bzcat -cd %s | wc -c"%(im))
        times[im].append( time.clock()-start )
        nt += 1 ; ns -= 1
        start = time.clock()
        the_file = bz2.BZ2File(im).read()
        times[im].append( time.clock()-start )  
        nt += 1; ns -= 1
    # Speed ratings in megabytes per second
    len(the_file) / 1024 / 1024

    print ("%.4f "*nt + " "*7*ns)%tuple(times[im]), im
    