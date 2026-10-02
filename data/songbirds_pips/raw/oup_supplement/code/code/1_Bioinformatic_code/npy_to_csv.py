import numpy
import sys
path=sys.argv[1]
output=sys.argv[2]
weights=numpy.load(path)
numpy.savetxt(output, weights, delimiter=",")