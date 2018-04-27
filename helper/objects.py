# vertex/vector class
class Vertex:
    def __init__(self, x, y, z, f=1):
        self.x = x
        self.y = y
        self.z = z
        self.f = f
    
    def getx(self):
        return self.x
    def gety(self):
        return self.y
    def getz(self):
        return self.z
    def getf(self):
        return self.f
    def printVertex(self):
        print '({}, {}, {}, {})'.format(self.x, self.y, self.z, self.f)

# surface class
class Surface:
    def __init__(self, vertices, numVertices, surfaceNormal=None):
        self.vertices = vertices
        self.numVertices = numVertices
        self.triIndex = 0
        self.triangles = []
        self.normal = surfaceNormal
    
    def getVertices(self):
        return self.vertices
    def getVerticesLen(self):
        return self.numVertices
    def getVertex(self, index):
        return self.vertices[index]
    def getTriangles(self):
        return self.triangles
    def getNumTriangles(self):
        return self.triIndex
    def getSurfaceNormal(self):
        return self.normal
    def addTriangle(self, triangle):
        self.triangles.append(triangle)
        self.triIndex = self.triIndex + 1
    def printSurface(self):
        for vertex in self.vertices:
            vertex.printVertex()

# polygon class
class Polygon:
    def __init__(self, numSurfaces, surfaces):
        self.numSurfaces = numSurfaces
        self.surfaces = surfaces
        self.numSides = numSurfaces - 2
    
    def getNumSurfaces(self):
        return self.numSurfaces
    def getNumSides(self):
        return self.numSides
    def getSurfaceArr(self):
        return self.surfaces
    def getSurface(self, index):
        return self.surfaces[index]
    def addSurface(self, surface):
        self.surfaces.append(surface)
       
# basis class
class Basis:
    def __init__(self, i, j, k):
        self.i = i
        self.j = j
        self.k = k
    
    def geti(self):
        return self.i
    def getj(self):
        return self.j
    def getk(self):
        return self.k
    def printBasis(self):
        self.i.printVertex()
        self.j.printVertex()
        self.k.printVertex()

# view point class
class ViewPoint:
    def __init__(self, c, u, v, n):
        self.c = c
        self.u = u
        self.v = v
        self.n = n
    
    def getc(self):
        return self.c
    def getu(self):
        return self.u
    def getv(self):
        return self.v
    def getn(self):
        return self.n
    def printViewPoint(self):
        self.c.printVertex()
        self.u.printVertex()
        self.v.printVertex()
        self.n.printVertex()

# view volume class
class ViewVolume:
    def __init__(self, nearZ, farZ, height):
        self.d = nearZ
        self.f = farZ
        self.h = height
        self.nh = height * -1
    def getd(self):
        return self.d
    def getf(self):
        return self.f
    def geth(self):
        return self.h
    def printViewVolume(self):
        print self.d,self.f,self.h
    def inViewVolume(self, v):
        lowerBound = (self.nh / self.d) * v.getz()
        upperBound = (self.h / self.d) * v.getz()  

        tmp = 0
        if(upperBound < lowerBound):
            tmp = upperBound
            upperBound = lowerBound
            lowerBound = tmp

        if(v.getx() >= lowerBound) and (v.getx() <= upperBound):
            pass
        else:
            return 0    
        if(v.gety() >= lowerBound) and (v.gety() <= upperBound):
            pass
        else:
            return 0
        if(v.getz() >= self.d) and (v.getz() <= self.f):
            pass
        else:
            return 0
        return 1


