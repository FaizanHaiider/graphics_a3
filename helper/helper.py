#!/usr/bin/python

from objects import *
from auxFunctions import *
from math import *

# create viewPoint
def createViewPoint(sphCoords, origin, basis):
    # required variables
    r = sphCoords[0]
    elvation = radians(sphCoords[1])
    azi = radians(sphCoords[2])

    # calc C
    # calc x = r * cos(azi) * cos(elvation)
    x = r * cos(elvation) * cos(azi)
    # calc y = r * cos(azi) * sin(elvation)
    y = r * sin(azi) * cos(elvation)
    # calc z = r * cos(azi)
    z = r * sin(elvation)
    c = vertex(x,y,z)

    # calc N
    n = createVector(c, origin)
    n = normalize(n)

    # calc V
    k = basis.getk()
    kn = calcInnerProd(k, n)
    knn = scalarMult(kn, n)    
    vPrime = subtrVector(k, knn)
    v = normalize(vPrime)

    # calc U
    u = calcCrossProduct(n, v)

    viewpoint = viewPoint(c,u,v,n)
    return viewpoint
# end createViewPoint

# create surfaces for polygon
def createPolygonSurfaces(newPolygon, imgSize, origin):
    numSides = newPolygon.getNumSides()
    zVal = imgSize/2
    # create front face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, zVal)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # calculate z val for rear face
    zVal = zVal * -1

    # create rear face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, zVal)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # create intermediate surfaces
    for i in range(numSides):
        newPolygon.addSurface(createInterSurfaces(newPolygon, i))
    return newPolygon
# end createPolygonSurfaces

# create triangle mesh
def createTriangleMesh(polygon, imgSize, triMesh):
    numSurfaces = polygon.getNumSurfaces()
    numSides = polygon.getNumSides()
    startIndex = 0
    
    # simple tringle mesh for front/rear faces
    if(numSides > 4):
        createSimpleTriangleMesh(polygon.getSurface(0), imgSize)
        createSimpleTriangleMesh(polygon.getSurface(1), imgSize)
        startIndex = 2       
# end createTriangleMesh

# rasterize polygon
def rasterize(Matrix, poly, imgSize, viewpoint, origin):
    # modify all points so they are relative to the view point

    M, rM = createAugMatrix(viewpoint)
    det, rM = gauss(M, rM)
    newPoly = updatePolyPOV(rM, poly)

    for surface in newPoly.getSurfaceArr():
        verticesLen = surface.getVerticesLen()
        # regular edges
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            v = surface.getVertex(i)
            vv = surface.getVertex(ii)
            BresenhamAlgo(v, vv, imgSize, Matrix, origin)
        # triangle mesh rasterize
        # for triangle in surface.getTriangles():
        #     v = triangle.getVertex(0)
        #     vv = triangle.getVertex(1)
        #     vvv = triangle.getVertex(2)
        #     BresenhamAlgo(v, vv, imgSize, Matrix)
        #     BresenhamAlgo(vv, vvv, imgSize, Matrix)
        #     BresenhamAlgo(vvv, v, imgSize, Matrix)
# end rasterize

# write Matrix to pbm file
def writePPM(Matrix, filename, imgSize):
    # create output file
    with open(filename, "w+") as fp:
        for i in Matrix:
            fp.write(' '.join(map(str, i)))

    # read, format, write file data
    with open(filename, "r+") as fp:
        lines = fp.read()
    lines = lines.replace("00", "0 0")
    lines = lines.replace("01", "0 1")
    lines = lines.replace("10", "1 0")
    lines = lines.replace("11", "1 1")
    with open(filename, "w+") as fp:
        fp.write("P1\n")
        fp.write('{} {}\n'.format(imgSize, imgSize))
        fp.write(lines)
# end writePPM
