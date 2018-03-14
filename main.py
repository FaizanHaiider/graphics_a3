#!/usr/bin/python

from helper.helper import *
from helper.objects import *
from math import *

def main():
    origin = vertex(0,0,0)
    filename, imgSize, numShapes, shapeArr, triMesh = readFileData("input.txt")
    Matrix = [[0 for x in range(imgSize)] for y in range(imgSize)]
    
    for i in range(numShapes):

        # create new polygon 
        numSurfaces = shapeArr[i] + 2
        surfaces = [0] * numSurfaces
        newPolygon = polygon(numSurfaces, surfaces)

        # create polygon surfaces
        createPolygonSurfaces(newPolygon, imgSize)

        # create triangle mesh
        createTriangleMesh(newPolygon, imgSize, triMesh)

        # update Matrix
        rasterize(Matrix, newPolygon, imgSize)

    # write to ppm file
    writePPM(Matrix, filename, imgSize)

# end main

# 
if __name__ == "__main__":
    main()
# end starting script