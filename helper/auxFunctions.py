#!/usr/bin/python

from objects import *
from math import *
import copy
from fractions import Fraction

# read from input file
def readFileData(filename):
    input_data = []
    # read all lines
    with open(filename, "r") as f:
        input_data = f.readlines()
    
    # fileInfo = [filename, image size]
    filename = input_data[0].strip("\r\n")
    imgSize = int(input_data[1])
    fileInfo = [filename, imgSize]

    # shapeInfo = [shape1, shape2, ..., triMesh]
    shapeInfo = getShapeInfo(input_data[2].strip("\r\n"))
    triMesh = int(input_data[3].strip("\r\n"))
    shapeInfo.append(triMesh)   

    # spherical coordinates = [azi, elvation, r]   
    sphCoords = input_data[4].strip("\r\n").split(",")
    sphCoords[0], sphCoords[1], sphCoords[2] = float(sphCoords[0]), float(sphCoords[1]), float(sphCoords[2])
   
    # frustum dimensions = [far_z, near_z, height]
    frustumDim = input_data[5].strip("\r\n").split(",")
    frustumDim[0], frustumDim[1], frustumDim[2] = float(frustumDim[0]), float(frustumDim[1]), float(frustumDim[2]) * 2 
     
    return fileInfo, shapeInfo, sphCoords, frustumDim
# end readFileData

# format shapes data from input file
def getShapeInfo(inputStr):
    # get number of shapes
    numShapes = inputStr.count(',') + 1
    shapeArr = []

    if(numShapes) == 1:
        shapeArr.append(int(inputStr))
        return shapeArr
    
    inputStr = inputStr.split(",")
    for i in range(numShapes):
        shapeArr.append(int(inputStr[i]))
    return shapeArr
# end getShapeInfo

# calculate cross product two vectors
def calcCrossProduct(vOne, vTwo):
    x = (vOne.gety() * vTwo.getz()) - (vTwo.gety() * vOne.getz())
    y = (vOne.getz() * vTwo.getx()) - (vTwo.getz() * vOne.getx())
    z = (vOne.getx() * vTwo.gety()) - (vTwo.getx() * vOne.gety())
    newVertex = Vertex(x,y,z)
    return newVertex
# end calcCrossProd

# create vector using two vertices
def createVector(vOne, vTwo):
    x = vTwo.getx() - vOne.getx()
    y = vTwo.gety() - vOne.gety()
    z = vTwo.getz() - vOne.getz()

    return Vertex(x,y,z)
# end createVector

# find norm of vertex
def calcNorm(v):
    return sqrt(calcInnerProd(v,v))
# end calcNorm 

# scalar multiplication
def scalarMult(const, v):
    nv = Vertex(v.getx() * const, v.gety() * const, v.getz() * const)
    return nv
# end scalarMult

# normalize unit vector
def normalize(v):
    # find norm
    norm = calcNorm(v)
    # constant
    const = (1/norm)
    # normalize
    return scalarMult(const, v)
# end normalize

# calc inner product
def calcInnerProd(v, vv):
    x = v.getx() * vv.getx()
    y = v.gety() * vv.gety()
    z = v.getz() * vv.getz()
    return x+y+z
# end calcInnerProd

# subtract two vectors
def subtrVector(v, vv):
    nvv = Vertex(vv.getx() * (-1), vv.gety() * (-1), vv.getz() * (-1))
    fv = Vertex(v.getx() + nvv.getx(), v.gety() + nvv.gety(), v.getz() + nvv.getz())
    return fv
# end subtrVector

# square function
def square(num):
    return num**2
# end square

# calculate surface vertices and scale
def calcSurfaceVertices(imgSize, numSides, zVal):
    surfaceVertices = []
    for i in range(numSides):
        x = ((imgSize/2) * cos((2 * pi * i) / numSides) + (imgSize/2)) / imgSize
        y = ((imgSize/2) * sin((2 * pi * i) / numSides) + (imgSize/2)) / imgSize
        newVertex = translateVertex(Vertex(-0.5, -0.5, 0.0), Vertex(x,y,zVal))
        surfaceVertices.append(newVertex)
    return surfaceVertices 
# end calcSurfaceVertices

# translate a vertex
def translateVertex(tv, v):
    rowx = Vertex(1,0,0,tv.getx())
    rowy = Vertex(0,1,0,tv.gety())
    rowz = Vertex(0,0,1,tv.getz())
    rowf = Vertex(0,0,0,1)

    x = (v.getx() * rowx.getx()) + (v.gety() * rowx.gety()) + (v.getz() * rowx.getz()) + (v.getf() * rowx.getf())
    y = (v.getx() * rowy.getx()) + (v.gety() * rowy.gety()) + (v.getz() * rowy.getz()) + (v.getf() * rowy.getf())
    z = (v.getx() * rowz.getx()) + (v.gety() * rowz.gety()) + (v.getz() * rowz.getz()) + (v.getf() * rowz.getf())
    f = (v.getx() * rowf.getx()) + (v.gety() * rowf.gety()) + (v.getz() * rowf.getz()) + (v.getf() * rowf.getf())

    return Vertex(x,y,z,f)
# end translateVertex

# scale vertex
def scaleVertex(tv, v):
    rowx = Vertex(tv.getx(),0,0,0)
    rowy = Vertex(0,tv.gety(),0,0)
    rowz = Vertex(0,0,tv.getz(),0)
    rowf = Vertex(0,0,0,1)

    x = (v.getx() * rowx.getx()) + (v.gety() * rowx.gety()) + (v.getz() * rowx.getz()) + (v.getf() * rowx.getf())
    y = (v.getx() * rowy.getx()) + (v.gety() * rowy.gety()) + (v.getz() * rowy.getz()) + (v.getf() * rowy.getf())
    z = (v.getx() * rowz.getx()) + (v.gety() * rowz.gety()) + (v.getz() * rowz.getz()) + (v.getf() * rowz.getf())
    f = (v.getx() * rowf.getx()) + (v.gety() * rowf.gety()) + (v.getz() * rowf.getz()) + (v.getf() * rowf.getf())

    return Vertex(x,y,z,f)
# end scaleVertex

# rotate around x-axis
def rotateZ(angle, v):
    newAngle = radians(angle)
    rowx = Vertex(cos(newAngle), (-1) * sin(newAngle), 0, 0)
    rowy = Vertex(sin(newAngle), cos(newAngle), 0, 0)
    rowz = Vertex(0, 0, 1, 0)
    rowf  = Vertex(0, 0, 0, 1)

    x = (v.getx() * rowx.getx()) + (v.gety() * rowx.gety()) + (v.getz() * rowx.getz()) + (v.getf() * rowx.getf())
    y = (v.getx() * rowy.getx()) + (v.gety() * rowy.gety()) + (v.getz() * rowy.getz()) + (v.getf() * rowy.getf())
    z = (v.getx() * rowz.getx()) + (v.gety() * rowz.gety()) + (v.getz() * rowz.getz()) + (v.getf() * rowz.getf())
    f = (v.getx() * rowf.getx()) + (v.gety() * rowf.gety()) + (v.getz() * rowf.getz()) + (v.getf() * rowf.getf())

    return Vertex(x,y,z,f)
# end rotateX

# rotate around x-axis
def rotateY(angle, v):
    newAngle = radians(angle)
    rowx = Vertex(cos(newAngle), 0, sin(newAngle), 0)
    rowy = Vertex(0, 1, 0, 0)
    rowz = Vertex((-1) * sin(newAngle), 0, cos(newAngle), 0)
    rowf  = Vertex(0, 0, 0, 1)

    x = (v.getx() * rowx.getx()) + (v.gety() * rowx.gety()) + (v.getz() * rowx.getz()) + (v.getf() * rowx.getf())
    y = (v.getx() * rowy.getx()) + (v.gety() * rowy.gety()) + (v.getz() * rowy.getz()) + (v.getf() * rowy.getf())
    z = (v.getx() * rowz.getx()) + (v.gety() * rowz.gety()) + (v.getz() * rowz.getz()) + (v.getf() * rowz.getf())
    f = (v.getx() * rowf.getx()) + (v.gety() * rowf.gety()) + (v.getz() * rowf.getz()) + (v.getf() * rowf.getf())

    return Vertex(x,y,z,f)
# end rotateX

# rotate around x-axis
def rotateX(angle, v):
    newAngle = radians(angle)
    rowx = Vertex(1, 0, 0, 0)
    rowy = Vertex(0, cos(newAngle), (-1) * sin(newAngle), 0)
    rowz = Vertex(0, sin(newAngle), cos(newAngle), 0)
    rowf  = Vertex(0, 0, 0, 1)

    x = (v.getx() * rowx.getx()) + (v.gety() * rowx.gety()) + (v.getz() * rowx.getz()) + (v.getf() * rowx.getf())
    y = (v.getx() * rowy.getx()) + (v.gety() * rowy.gety()) + (v.getz() * rowy.getz()) + (v.getf() * rowy.getf())
    z = (v.getx() * rowz.getx()) + (v.gety() * rowz.gety()) + (v.getz() * rowz.getz()) + (v.getf() * rowz.getf())
    f = (v.getx() * rowf.getx()) + (v.gety() * rowf.gety()) + (v.getz() * rowf.getz()) + (v.getf() * rowf.getf())

    return Vertex(x,y,z,f)
# end rotateX

# create surfaces between front/rear faces
def createInterSurfaces(newPolygon, i):
    numSides = newPolygon.getNumSides()
    tmpFront = newPolygon.getSurface(0)
    tmpRear = newPolygon.getSurface(1) 

    ii = i + 1
    if(ii == numSides):
        ii = 0
    vertices = [tmpRear.getVertex(i), tmpRear.getVertex(ii), tmpFront.getVertex(ii), tmpFront.getVertex(i)] 
    return Surface(vertices, 4)
# end createInterSurfaces 

# create aug matrix for gauss
def createAugMatrix(viewpoint):
    M = [[0 for x in range(3)] for y in range(3)]
    M[0][0] = viewpoint.getu().getx()
    M[1][0] = viewpoint.getu().gety()
    M[2][0] = viewpoint.getu().getz()

    M[0][1] = viewpoint.getv().getx()
    M[1][1] = viewpoint.getv().gety()
    M[2][1] = viewpoint.getv().getz()

    M[0][2] = viewpoint.getn().getx()
    M[1][2] = viewpoint.getn().gety()
    M[2][2] = viewpoint.getn().getz()

    resultM = [[0 for x in range(4)] for y in range(3)]

    resultM[0][3] = viewpoint.getc().getx() * -1
    resultM[1][3] = viewpoint.getc().gety() * -1
    resultM[2][3] = viewpoint.getc().getz() * -1 

    resultM[0][0], resultM[1][1], resultM[2][2] = 1, 1, 1
    resultM[0][1], resultM[0][2], resultM[1][0], resultM[1][2], resultM[2][0], resultM[2][1] = 0, 0, 0, 0, 0, 0
    
    return M, resultM
# end createAugMatrix

# update polygon coords with new POV coords
def updatePolyPOV(rM, poly):
    numSurfaces = poly.getNumSurfaces()
    surfaces = []
    newPoly = Polygon(numSurfaces, surfaces)

    for surface in poly.getSurfaceArr():
        numVertices = surface.getVerticesLen()
        vertices = []
        for v in surface.getVertices():
            nx, ny, nz = 0,0,0
            nx = (rM[0][0] * v.getx()) + (rM[0][1] * v.gety()) + (rM[0][2] * v.getz()) + rM[0][3]
            ny = (rM[1][0] * v.getx()) + (rM[1][1] * v.gety()) + (rM[1][2] * v.getz()) + rM[1][3]
            nz = (rM[2][0] * v.getx()) + (rM[2][1] * v.gety()) + (rM[2][2] * v.getz()) + rM[2][3]
            vertices.append(Vertex(nx, ny, nz))
        newSurface = Surface(vertices, numVertices)
        newPoly.addSurface(newSurface)
    return newPoly  
# end updatePolyPOV

# Gauss-Jordan Elimination
def gauss(a, b):
    a = copy.deepcopy(a)
    b = copy.deepcopy(b)
    n = len(a)
    p = len(b[0])
    det = 1
    for i in range(n - 1):
        k = i
        for j in range(i + 1, n):
            if abs(a[j][i]) > abs(a[k][i]):
                k = j
        if k != i:
            a[i], a[k] = a[k], a[i]
            b[i], b[k] = b[k], b[i]
            det = -det
 
        for j in range(i + 1, n):
            t = a[j][i]/a[i][i]
            for k in range(i + 1, n):
                a[j][k] -= t*a[i][k]
            for k in range(p):
                b[j][k] -= t*b[i][k]
 
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            t = a[i][j]
            for k in range(p):
                b[i][k] -= t*b[j][k]
        t = 1/a[i][i]
        det *= a[i][i]
        for j in range(p):
            b[i][j] *= t
    return det, b
 
def zeromat(p, q):
    return [[0]*q for i in range(p)]
 
def matmul(a, b):
    n, p = len(a), len(a[0])
    p1, q = len(b), len(b[0])
    if p != p1:
        raise ValueError("Incompatible dimensions")
    c = zeromat(n, q)
    for i in range(n):
        for j in range(q):
                c[i][j] = sum(a[i][k]*b[k][j] for k in range(p))
    return c
 
def mapmat(f, a):
    return [list(map(f, v)) for v in a]
 
def ratmat(a):
    return mapmat(Fraction, a)
# end Gauss Jordon Elimination

# update Matrix after rasterize()
def updateMatrix(Matrix, points, swapped):
    if(swapped):
        points.reverse()
    for point in points:
        Matrix[point[0]][point[1]] = 1
# end updateMatrix

# check if polygon is in view volume
def isPolyInViewVolume(poly, viewVolume):
    for surface in poly.getSurfaceArr():
        for vertex in surface.getVertices():
            if(viewVolume.inViewVolume(vertex) == 0):
                return 0
    return 1
# end PolyInViewVolume

# project coordinates to view window
def ConvertToViewWindow(poly, viewVolume):
    numSurfaces = poly.getNumSurfaces()
    surfaces = []
    nPoly = Polygon(numSurfaces, surfaces)

    for surface in poly.getSurfaceArr():
        numVertices = surface.getVerticesLen()
        vertices = []
        for v in surface.getVertices():
            nx, ny, nz = 0,0,0
            nx = (v.getx() / v.getz()) * viewVolume.getd()
            ny = (v.gety() / v.getz()) * viewVolume.getd()
            z = viewVolume.getd()
            vertices.append(Vertex(nx, ny, nz))
        surfaceNormal = calcSurfaceNormal(vertices)
        newSurface = Surface(vertices, numVertices, surfaceNormal)
        nPoly.addSurface(newSurface)
    return nPoly
# end ConvertToViewVolume

# calculate surface's normal
def calcSurfaceNormal(vertexArr):
    p, pp, ppp = vertexArr[0], vertexArr[1], vertexArr[2]
    v, vv = createVector(p, pp), createVector(p, ppp)

    return normalize(calcCrossProduct(v, vv))
# end calcSurfaceNormal

# is surface visible 
def isSurfaceVisible(surface, viewPoint):
    p = surface.getVertex(0)
    normal = surface.getSurfaceNormal()

    cp = createVector(viewPoint.getc(), p)
    result = calcInnerProd(cp, normal)

    if(result >= 0.0):
        return 0
    return 1
# end isSurfaceVisible 

