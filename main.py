#!/usr/bin/python

from math import *
from helper.mainFunctions import *
from helper.objects import *
from helper.auxFunctions import *

def main():
    # read file inputs and create output files
    fileInfo, shapeInfo, sphCoords, frustumDim = readFileData("input.txt")
    filename = fileInfo[0] + ".pbm"
    filename_before = filename + "_before.pbm"

    # set up matrices
    Matrix = [[0 for x in range(fileInfo[1])] for y in range(fileInfo[1])]
    Matrix_before = [[0 for x in range(fileInfo[1])] for y in range(fileInfo[1])]
    
    # set up scene
    origin = Vertex(0,0,0)
    basis = Basis(Vertex(1,0,0), Vertex(0,1,0), Vertex(0,0,1))
    viewPoint = createViewPoint(sphCoords, origin, basis)
    viewVolume = ViewVolume(frustumDim[0], frustumDim[1], frustumDim[2])

    # used to create scene
    scaleBy = [Vertex(0.5, 0.5, 0.5), Vertex(0.5, 0.5, 0.5), Vertex(0.5, 0.5, 0.5)]
    transitions = [Vertex(0.5, 0.0, 0.0), Vertex(-0.5, 0.0, 0.0), Vertex(0,0,0)]
    index = 0
    
    # create shape(s)
    polygonArr = []
    for i in range(len(shapeInfo) - 1):
        
        # initialize new polygon
        numSurfaces = shapeInfo[i] + 2
        surfaces = []
        newPoly = Polygon(numSurfaces, surfaces)

        # create polygon surfaces
        newPoly = createPolySurfaces(newPoly, fileInfo[1], origin)

        # create scene
        poly = createScene(newPoly, transitions[index], scaleBy[index])
        index = index + 1

        # rasterize
        rasterize(Matrix, Matrix_before, poly, fileInfo[1], viewPoint, origin, viewVolume)

    # write to before and after files
    writePPM(Matrix_before, filename_before, fileInfo[1])
    writePPM(Matrix, filename, fileInfo[1])

# end main

if __name__ == '__main__':
    main()
