#!/usr/bin/python

from helper.helper import *
from helper.objects import *
from math import *

def main():
    # read input file
    filename, imgSize, numShapes, shapeArr, triMesh, sphCoords = readFileData("input.txt")

    # init variables
    Matrix = [[0 for x in range(imgSize)] for y in range(imgSize)]
    basis = Basis(vertex(imgSize/2,0,0), vertex(0,imgSize/2,0), vertex(0,0,imgSize/2))
    origin = vertex(0,0,0)
    viewpoint = createViewPoint(sphCoords, origin, basis)

    # create shape(s)
    for i in range(numShapes):

        # create new polygon 
        numSurfaces = shapeArr[i] + 2
        surfaces = [0] * numSurfaces
        newPolygon = polygon(numSurfaces, surfaces)

        # create polygon surfaces
        newPolygon = createPolygonSurfaces(newPolygon, imgSize, origin)

        # create triangle mesh
        # createTriangleMesh(newPolygon, imgSize, triMesh)

        # update Matrix
        rasterize(Matrix, newPolygon, imgSize, viewpoint, origin)

    # write to ppm file
    writePPM(Matrix, filename, imgSize)

# end main

# 
if __name__ == "__main__":
    main()
# end starting script