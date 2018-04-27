#!/usr/bin/python

from objects import *
from auxFunctions import *
from math import *

# create viewPoint
def createViewPoint(sphCoords, origin, basis):
    # (azi, elvation, r)
    azi = radians(sphCoords[0])
    elvation = radians(sphCoords[1])
    r = sphCoords[2]

    # calc C
    x = r * cos(elvation) * cos(azi)
    y = r * cos(elvation) * sin(azi)
    z = r * sin(elvation)
    c = Vertex(x,y,z)

    # calc N
    tn = createVector(c, origin)
    n = normalize(tn)

    # calc V
    k = basis.getk()
    kn = calcInnerProd(k, n)
    knn = scalarMult(kn, n)   
    vPrime = subtrVector(k, knn)
    v = normalize(vPrime)

    # calc U
    u = calcCrossProduct(n, v)

    viewpoint = ViewPoint(c,u,v,n)
    return viewpoint
# end createViewPoint

# create surfaces for polygon
def createPolySurfaces(newPolygon, imgSize, origin):
    numSides = newPolygon.getNumSides()
    
    # create initial front surface vertices
    frontVertices = calcSurfaceVertices(imgSize, numSides, 0.5)
    newPolygon.addSurface(Surface(frontVertices, numSides))

    # create initial rear surface vertices
    rearVertices = calcSurfaceVertices(imgSize, numSides, -0.5)
    
    # reverse rear vertices for surface normal
    
    # translate rear surface to origin
    newPolygon.addSurface(Surface(rearVertices, numSides))

    # create intermediate surfaces
    for i in range(numSides):
        newPolygon.addSurface(createInterSurfaces(newPolygon, i))
    return newPolygon
# end createPolygonSurfaces

# create scene
def createScene(poly, v, vv):
    numSurfaces = poly.getNumSurfaces()
    surfaces = []
    newPoly = Polygon(numSurfaces, surfaces)

    for surface in poly.getSurfaceArr():
        newVertices = []
        for vertex in surface.getVertices():
            newVertices.append(translateVertex(v, scaleVertex(vv, vertex)))
        newSurface = Surface(newVertices, surface.getVerticesLen())
        newPoly.addSurface(newSurface)

    return newPoly    
# end createScene

# rasterize polygon
def rasterize(Matrix, Matrix_before, poly, imgSize, viewPoint, origin, viewVolume):
    # is polygon in view volume
    if(isPolyInViewVolume(poly, viewVolume) == 0):
        return

    # final scale for scene
    scaleBy = Vertex(imgSize/2,imgSize/2,imgSize/2)

    # modify all points so they are relative to the view 
    M, rM = createAugMatrix(viewPoint)
    det, rM = gauss(M, rM)
    nPoly = updatePolyPOV(rM, poly)
    nnPoly = ConvertToViewWindow(nPoly, viewVolume)

    createBeforeImg(Matrix_before, nnPoly, imgSize, origin)

    for surface in nnPoly.getSurfaceArr():
        # culling
        if(isSurfaceVisible(surface, viewPoint) == 0):
            continue 
        verticesLen = surface.getVerticesLen()
        # regular edges
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            ov = surface.getVertex(i)
            ovv = surface.getVertex(ii)

            # scale to fit in Matrix
            v = scaleVertex(scaleBy, ov)
            vv = scaleVertex(scaleBy, ovv)

            nv = translateVertex(scaleBy, v)
            nvv = translateVertex(scaleBy, vv)

            BresenhamAlgo(nv, nvv, imgSize, Matrix, origin)
# end rasterize

# create before image with no culling
def createBeforeImg(Matrix_before, poly, imgSize, origin):
    scaleBy = Vertex(imgSize/2, imgSize/2, imgSize/2)

    for surface in poly.getSurfaceArr():
        verticesLen = surface.getVerticesLen()
        # regular edges
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            ov = surface.getVertex(i)
            ovv = surface.getVertex(ii)

            # scale to fit in Matrix
            v = scaleVertex(scaleBy, ov)
            vv = scaleVertex(scaleBy, ovv)

            nv = translateVertex(scaleBy, v)
            nvv = translateVertex(scaleBy, vv)

            BresenhamAlgo(nv, nvv, imgSize, Matrix_before, origin)
# end createBeforeImg

# Bresenhams's Algo
def BresenhamAlgo(v, vv, imgSize, Matrix, origin):
    xs, ys = int(v.getx()), int(v.gety())
    xe, ye = int(vv.getx()), int(vv.gety())

    dx = xe - xs
    dy = ye - ys
    isSteep = abs(dy) > abs(dx)

    if(isSteep):
        xs, ys = ys, xs
        xe, ye = ye, xe
    swapped = 0
    if(xs > xe):
        xs, xe = xe, xs
        ys, ye = ye, ys
        swapped = 1
    
    dx = xe - xs
    dy = ye - ys
    error = int(dx/2.0)
    yincr = 1 if(ys < ye) else -1

    y = ys
    points = []
    for x in range(xs, xe + 1):
        tx = (y if isSteep else x)
        ty = (x if isSteep else y)
        tx = (tx if tx < imgSize else tx-1)
        ty = (ty if ty < imgSize else ty-1)
        coord = (tx, ty)
        points.append(coord)
        error  = error - abs(dy)
        if(error < 0):
            y = y + yincr
            error = error + dx
    updateMatrix(Matrix, points, swapped)    
# end BresenhamAlgo

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