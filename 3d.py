# 3d.py, by Miguel Enzo Aruelo: https://github.com/migearu

# cool computer graphics programming test
# seriously though, computer graphics programming is interesting

# divided it up into different sections, for example, since CMU doesn't provide classes for 2d and 3d points, along with matrices, i had to implement them on my own
# then using those classes, i created the basic variables needed for the rest of the program
# then i made the actual polyhedra, then added stuff to them, starting with a numerical representation, then placing that representation on the screen, then adding faces and culling etc.
# after that, is the main section, which just sets up the program, renders everything, and provides controls.

import math
import random
import time

# NOTES

# this section of code below is used to calculate vertex connections in the regular polyhedra (except for cubes)

# lineSegmentsFound = 0
# lineSegmentsFoundRepeatable = 0
# foundArray = []
# for i in range(len(TDPV)):
#     for j in range(len(TDPV)):
#         if i == j:
#             continue
#         distance = Vector3.FindDistance(self.vertices[i], self.vertices[j])
#         if (distance - 2*self.radius < 0.1 and distance - 2*self.radius > -0.1) or (distance + 2*self.radius < 0.1 and distance + 2*self.radius > -0.1):
#             if not [j, i] in foundArray:
#                 lineSegmentsFound += 1
#                 foundArray.append([i, j])
#             lineSegmentsFoundRepeatable += 1
#             if lineSegmentsFoundRepeatable >= 5:
#                 continue
# print(foundArray)

# everything here is basically math

appdata = {
    'steps': 0,
    'drag': {
        'pX': 0,
        'pY': 0,
        'X': 0,
        'Y': 0,
        'dX': 0,
        'dY': 0
    },
    'globalWireframe': True,
    'globalFaceRender': True,
    'globalDebug': False
}

#------------------------- MATHEMATICAL OBJECT CLASSES -------------------------#

class Matrix:
    def __init__(self, matrix):
        self.matrix = matrix
    
    def Multiply(m1, m2):
        result = []
        for i in range(len(m1.matrix)):
            row = []
            for j in range(len(m2.matrix[0])):
                sum = 0
                for k in range(len(m2.matrix)):
                    sum += m1.matrix[i][k] * m2.matrix[k][j]
                row.append(sum)
            result.append(row)
        return Matrix(result)
    
    def ColumnVectorToVector3(self):
        if len(self.matrix[0]) > 1:
            raise AttributeError('A column vector must only contain 1 entry per row.')
        if len(self.matrix) != 3:
            raise AttributeError('A column vector must contain 3 rows for conversion to Vector3 to be possible')
        return Vector3(self.matrix[0][0], self.matrix[1][0], self.matrix[2][0])
    
    def RotationMatrix(rotVect):
        Rx = math.radians(rotVect.z)
        Ry = math.radians(rotVect.y)
        Rz = math.radians(rotVect.x)
        return Matrix([
            [
                math.cos(Ry) * math.cos(Rz),
                math.sin(Rx) * math.sin(Ry) * math.cos(Rz) - math.cos(Rx) * math.sin(Rz),
                math.cos(Rx) * math.sin(Ry) * math.cos(Rz) + math.sin(Rx) * math.sin(Rz)
            ],
            [
                math.cos(Ry) * math.sin(Rz),
                math.sin(Rx) * math.sin(Ry) * math.sin(Rz) + math.cos(Rx) * math.cos(Rz),
                math.cos(Rx) * math.sin(Ry) * math.sin(Rz) - math.sin(Rx) * math.cos(Rz)
            ],
            [
                -math.sin(Ry),
                math.sin(Rx) * math.cos(Ry),
                math.cos(Rx) * math.cos(Ry)
            ]
        ])

class Vector2:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def add(self, vect):
        self.x += vect.x
        self.y += vect.y
        return self
        
    def mul(self, num):
        self.x *= num
        self.y *= num
        return self

class Vector3:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def add(self, vect):
        self.x += vect.x
        self.y += vect.y
        self.z += vect.z
        return self
        
    def sub(self, vect):
        self.x -= vect.x
        self.y -= vect.y
        self.z -= vect.z
        return self
    
    def mul(self, num):
        self.x *= num
        self.y *= num
        self.z *= num
        return self
        
    def dotProduct(self, vect):
        return (self.x * vect.x) + (self.y * vect.y) + (self.z * vect.z)
    
    def __add__(self, vect):
        return self.add(vect)
    
    def To2DProjection(self, camPos, camAngle):
        E = Vector3(0, 0, 1/math.tan(math.radians(fov)/2)) # display surface
        x = self.x - camPos.x
        y = self.y - camPos.y
        z = self.z - camPos.z
        Sx = math.sin(math.radians(camAngle.x))
        Sy = math.sin(math.radians(camAngle.y))
        Sz = math.sin(math.radians(camAngle.z))
        Cx = math.cos(math.radians(camAngle.x))
        Cy = math.cos(math.radians(camAngle.y))
        Cz = math.cos(math.radians(camAngle.z))
        
        Dx = Cy*(Sz*y + Cz*x) - Sy*z
        Dy = Sx*(Cy*z + Sy*(Sz*y + Cz*x)) + Cx*(Cz*y - Sz*x)
        Dz = Cx*(Cy*z + Sy*(Sz*y + Cz*x)) - Sx*(Cz*y - Sz*x)
        
        return Vector2((E.z/Dz)*Dx + E.x, (E.z/Dz)*Dy + E.y)
        
        # s = self.z*math.tan(math.radians(fov)/2)
        # return Vector2(self.x/s, self.y/s)
        
    def ToColumnVector(self):
        return Matrix([
            [self.x],
            [self.y],
            [self.z]
        ])
    
    def subtract(self, vect):
        self.x -= vect.x
        self.y -= vect.y
        self.z -= vect.z
        return self
    
    def FindDistance(vect1, vect2):
        return math.sqrt((vect2.x - vect1.x)**2 + (vect2.y - vect1.y)**2 + (vect2.z - vect1.z)**2)
    
    def clone(self):
        return Vector3(self.x, self.y, self.z)
        
    def invert(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z
        return self

#------------------------- MATHEMATICAL OBJECT CLASSES -------------------------#

#--------------------------- CONSTANTS ---------------------------#
# it is highly recommended you do not touch anything here.

goldenRatio = (1 + 5 ** 0.5) / 2

# each entry represents a line [i, j], which connects vertex i, to vertex j
dodecahedronEdges = [
    [0, 8], [0, 12], [0, 16], [1, 8], [1, 14],
    [1, 17], [2, 9], [2, 12], [2, 18], [3, 9],
    [3, 14], [3, 19], [4, 10], [4, 13], [4, 16],
    [5, 10], [5, 15], [5, 17], [6, 11], [6, 13],
    [6, 18], [7, 11], [7, 15], [7, 19], [8, 10],
    [9, 11], [12, 14], [13, 15], [16, 18], [17, 19],
]
dodecahedronFaces = [
    [16, 4, 10, 8, 0],
    [0, 12, 2, 18, 16],
    [0, 8, 1, 14, 12],
    [14, 1, 17, 19, 3],
    [14, 3, 9, 2, 12],
    [1, 8, 10, 5, 17],
    [19, 17, 5, 15, 7],
    [10, 4, 13, 15, 5],
    [3, 19, 7, 11, 9],
    [11, 7, 15, 13, 6],
    [11, 6, 18, 2, 9],
    [6, 13, 4, 16, 18]
]
icosahedronEdges = [
    [0, 1], [0, 4], [0, 5], [0, 8], [0, 10],
    [1, 6], [1, 7], [1, 8], [1, 10], [2, 3],
    [2, 4], [2, 5], [2, 9], [2, 11], [3, 6],
    [3, 7], [3, 9], [3, 11], [4, 5], [4, 8],
    [4, 9], [5, 10], [5, 11], [6, 7], [6, 8],
    [6, 9], [7, 10], [7, 11], [8, 9], [10, 11],
]
icosahedronFaces = [
    [11, 5, 2],
    [11, 2, 3],
    [11, 3, 7],
    [11, 7, 10],
    [11, 10, 5],
    [10, 0, 5],
    [5, 0, 4],
    [5, 4, 2],
    [2, 4, 9],
    [2, 9, 3],
    [3, 9, 6],
    [3, 6, 7],
    [9, 8, 6],
    [9, 4, 8],
    [4, 0, 8],
    [1, 8, 0],
    [1, 0, 10],
    [1, 10, 7],
    [1, 7, 6],
    [1, 6, 8]
]

#--------------------------- CONSTANTS ---------------------------#

#------------------------- RENDER OPTIONS -------------------------#

app.stepsPerSecond = 60 # fps
fov = 90
camPosition = Vector3(0, 0, 0)
camOrientation = Vector3(0, 0, 0)
center = Vector2(200, 200)
lightPos = Vector3(400, 0, 400)

#------------------------- RENDER OPTIONS -------------------------#

#------------------------- POLYHEDRON CLASSES -------------------------#

class RectPrism:
    #   4-------5
    #  /|      /|
    # 0-------1 |
    # | |     | |
    # | 6-----|-7
    # |/      |/
    # 2-------3
    def __init__(self, vtx1, vtx2, rotVect):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.vertices = [
            Vector3(vtx1.x, vtx1.y, vtx1.z),
            Vector3(vtx2.x, vtx1.y, vtx1.z),
            Vector3(vtx1.x, vtx2.y, vtx1.z),
            Vector3(vtx2.x, vtx2.y, vtx1.z),
            Vector3(vtx1.x, vtx1.y, vtx2.z),
            Vector3(vtx2.x, vtx1.y, vtx2.z),
            Vector3(vtx1.x, vtx2.y, vtx2.z),
            Vector3(vtx2.x, vtx2.y, vtx2.z),
        ]
        # take the average of all x, y and z pounts
        # subtract that to all vertices
        # rotate around that point
        # add back the average
        # boom, rotation.
        shapeCenter = Vector3((vtx1.x+vtx2.x)/2, (vtx1.y+vtx2.y)/2, (vtx1.z+vtx2.z)/2)
        for i in range(len(self.vertices)):
            vertex = self.vertices[i].subtract(shapeCenter)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(shapeCenter)
        self.vertex1 = vtx1
        self.vertex2 = vtx2
        self.rotationVector = rotVect
        self.renderedLineSegments = [Line(0, 0, 0, 0) for i in range(12)]
    
    def changeProperties(self, vtx1, vtx2, rotVect):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.vertices = [
            Vector3(vtx1.x, vtx1.y, vtx1.z),
            Vector3(vtx2.x, vtx1.y, vtx1.z),
            Vector3(vtx1.x, vtx2.y, vtx1.z),
            Vector3(vtx2.x, vtx2.y, vtx1.z),
            Vector3(vtx1.x, vtx1.y, vtx2.z),
            Vector3(vtx2.x, vtx1.y, vtx2.z),
            Vector3(vtx1.x, vtx2.y, vtx2.z),
            Vector3(vtx2.x, vtx2.y, vtx2.z),
        ]
        # take the average of all x, y and z pounts
        # subtract that to all vertices
        # rotate around that point
        # add back the average
        # boom, rotation.
        shapeCenter = Vector3((vtx1.x+vtx2.x)/2, (vtx1.y+vtx2.y)/2, (vtx1.z+vtx2.z)/2)
        for i in range(len(self.vertices)):
            vertex = self.vertices[i].subtract(shapeCenter)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(shapeCenter)
        self.vertex1 = vtx1
        self.vertex2 = vtx2
        self.rotationVector = rotVect
        
    def rotate(self, rotVect):
        shapeCenter = Vector3((self.vertex1.x+self.vertex2.x)/2, (self.vertex1.y+self.vertex2.y)/2, (self.vertex1.z+self.vertex2.z)/2)
        self.rotationVector += rotVect
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        for i in range(len(self.vertices)):
            vertex = self.vertices[i].subtract(shapeCenter)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(shapeCenter)
    
    def render(self):
        TDPV = [] # short for 2D Projected Vertices
        for vertex in self.vertices:
            TDPV.append(vertex.To2DProjection(camPosition, camOrientation).mul(400).add(center))
        # render section
        # Z lines
        for i in range(4):
            self.renderedLineSegments[i].x1 = TDPV[i].x
            self.renderedLineSegments[i].y1 = TDPV[i].y
            self.renderedLineSegments[i].x2 = TDPV[i+4].x
            self.renderedLineSegments[i].y2 = TDPV[i+4].y
        # X lines
        for i in range(4):
            self.renderedLineSegments[i+4].x1 = TDPV[i*2].x
            self.renderedLineSegments[i+4].y1 = TDPV[i*2].y
            self.renderedLineSegments[i+4].x2 = TDPV[i*2+1].x
            self.renderedLineSegments[i+4].y2 = TDPV[i*2+1].y
        # Y lines
        for i in range(2):
            self.renderedLineSegments[i+8].x1 = TDPV[i].x
            self.renderedLineSegments[i+8].y1 = TDPV[i].y
            self.renderedLineSegments[i+8].x2 = TDPV[i+2].x
            self.renderedLineSegments[i+8].y2 = TDPV[i+2].y
        for i in range(2):
            self.renderedLineSegments[i+10].x1 = TDPV[i+4].x
            self.renderedLineSegments[i+10].y1 = TDPV[i+4].y
            self.renderedLineSegments[i+10].x2 = TDPV[i+6].x
            self.renderedLineSegments[i+10].y2 = TDPV[i+6].y

class Dodecahedron:
    # not drawing it,,,
    def __init__(self, position, radius, rotVect, color=Vector3(135, 206, 235)):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.radius = radius
        self.position = position
        self.vertices = [
            Vector3(1, 1, 1),
            Vector3(-1, 1, 1),
            Vector3(1, -1, 1),
            Vector3(-1, -1, 1),
            Vector3(1, 1, -1),
            Vector3(-1, 1, -1),
            Vector3(1, -1, -1),
            Vector3(-1, -1, -1),
            Vector3(0, goldenRatio, 1/goldenRatio),
            Vector3(0, -goldenRatio, 1/goldenRatio),
            Vector3(0, goldenRatio, -1/goldenRatio),
            Vector3(0, -goldenRatio, -1/goldenRatio),
            Vector3(1/goldenRatio, 0, goldenRatio),
            Vector3(1/goldenRatio, 0, -goldenRatio),
            Vector3(-1/goldenRatio, 0, goldenRatio),
            Vector3(-1/goldenRatio, 0, -goldenRatio),
            Vector3(goldenRatio, 1/goldenRatio, 0),
            Vector3(-goldenRatio, 1/goldenRatio, 0),
            Vector3(goldenRatio, -1/goldenRatio, 0),
            Vector3(-goldenRatio, -1/goldenRatio, 0)
        ]
        for i in range(len(self.vertices)):
            vertex = self.vertices[i]
            vertex.mul(radius)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(position)
        self.renderedLineSegments = [Line(0, 0, 0, 0) for i in range(30)]
        self.renderedPolygons = [Polygon() for i in range(12)]
        self.debugLabels = [Label('undefined', 0, 0, size=20, opacity=0) for i in range(20)]
        self.rotationVector = rotVect
        for i in range(12):
            self.renderedPolygons[i].fill = rgb(randrange(0, 255), randrange(0, 255), randrange(0, 255))
        self.color = color
    
    def changeProperties(self, position, radius, rotVect, color=Vector3(135, 206, 235)):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.radius = radius
        self.position = position
        self.vertices = [
            Vector3(1, 1, 1),
            Vector3(-1, 1, 1),
            Vector3(1, -1, 1),
            Vector3(-1, -1, 1),
            Vector3(1, 1, -1),
            Vector3(-1, 1, -1),
            Vector3(1, -1, -1),
            Vector3(-1, -1, -1),
            Vector3(0, goldenRatio, 1/goldenRatio),
            Vector3(0, -goldenRatio, 1/goldenRatio),
            Vector3(0, goldenRatio, -1/goldenRatio),
            Vector3(0, -goldenRatio, -1/goldenRatio),
            Vector3(1/goldenRatio, 0, goldenRatio),
            Vector3(1/goldenRatio, 0, -goldenRatio),
            Vector3(-1/goldenRatio, 0, goldenRatio),
            Vector3(-1/goldenRatio, 0, -goldenRatio),
            Vector3(goldenRatio, 1/goldenRatio, 0),
            Vector3(-goldenRatio, 1/goldenRatio, 0),
            Vector3(goldenRatio, -1/goldenRatio, 0),
            Vector3(-goldenRatio, -1/goldenRatio, 0)
        ]
        for i in range(len(self.vertices)):
            vertex = self.vertices[i]
            vertex.mul(radius)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(position)
        self.rotationVector = rotVect
        self.color = color
    
    def rotate(self, rotVect):
        self.rotationVector += rotVect
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        for i in range(len(self.vertices)):
            vertex = self.vertices[i].subtract(self.position)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(self.position)
        
    def render(self, wireframe=True, faces=True, debug=False):
        TDPV = [] # short for 2D Projected Vertices
        for vertex in self.vertices:
            TDPV.append(vertex.To2DProjection(camPosition, camOrientation).mul(400).add(center))
        if debug:
            for i in range(len(self.vertices)):
                self.debugLabels[i].opacity = 100
                self.debugLabels[i].value = i
                self.debugLabels[i].centerX = TDPV[i].x
                self.debugLabels[i].centerY = TDPV[i].y
        else:
            for i in range(len(self.vertices)):
                self.debugLabels[i].opacity = 0
                self.debugLabels[i].centerX = 0
                self.debugLabels[i].centerY = 0
        if wireframe:
            for i in range(len(dodecahedronEdges)):
                self.renderedLineSegments[i].x1 = TDPV[dodecahedronEdges[i][0]].x
                self.renderedLineSegments[i].y1 = TDPV[dodecahedronEdges[i][0]].y
                self.renderedLineSegments[i].x2 = TDPV[dodecahedronEdges[i][1]].x
                self.renderedLineSegments[i].y2 = TDPV[dodecahedronEdges[i][1]].y
        else:
            for i in range(len(dodecahedronEdges)):
                self.renderedLineSegments[i].x1 = 0
                self.renderedLineSegments[i].y1 = 0
                self.renderedLineSegments[i].x2 = 0
                self.renderedLineSegments[i].y2 = 0
        if faces:
            for i in range(len(dodecahedronFaces)):
                pointList = []
                faceNormal = Vector3(0, 0, 0)
                faceCenter = Vector3(0, 0, 0)
                closestVertex = self.vertices[dodecahedronFaces[i][0]].clone()
                for k in range(5): # this section calculates the normal of the face (used for back-face culling), and the center (used for "lighting")
                    currentVtx = dodecahedronFaces[i][k]
                    nextVtx = dodecahedronFaces[i][(k+1)%5]
                    cvtxPos = self.vertices[currentVtx]
                    nvtxPos = self.vertices[nextVtx]
                    faceNormal.x += (cvtxPos.y - nvtxPos.y) * (cvtxPos.z + nvtxPos.z)
                    faceNormal.y += (cvtxPos.z - nvtxPos.z) * (cvtxPos.x + nvtxPos.x)
                    faceNormal.z += (cvtxPos.x - nvtxPos.x) * (cvtxPos.y + nvtxPos.y)
                    faceCenter.x += cvtxPos.x
                    faceCenter.y += cvtxPos.y
                    faceCenter.z += cvtxPos.z
                faceCenter.mul(1/5)
                # lots of painful math; basically had to teach myself quick linear algebra for a lot of this stuff
                # taking the dot product of the camera to vertex vector and the face normal vector can determine if a face's normal is facing away
                # if it is, that means that the face must be facing away from the camera, and then we can cull it
                # this saves processing power, and also assists in properly giving the full 3d effect.
                camToVtxVect = closestVertex.sub(camPosition)
                facingAway = (camToVtxVect.dotProduct(faceNormal) >= 0)
                if not facingAway:
                    for j in range(5):
                        pointList.append([TDPV[dodecahedronFaces[i][j]].x, TDPV[dodecahedronFaces[i][j]].y])
                    self.renderedPolygons[i].pointList = pointList
                    lightAmnt = Vector3.FindDistance(faceCenter, lightPos)/400
                    self.renderedPolygons[i].fill = rgb(
                        min(self.color.x/lightAmnt, 255),
                        min(self.color.y/lightAmnt, 255),
                        min(self.color.z/lightAmnt, 255)
                    )
                    #self.renderedPolygons[i].toBack()
                    #self.renderedPolygons[i].fill = rgb(randrange(0, 255), randrange(0, 255), randrange(0, 255))
                else:
                    for j in range(5):
                        pointList.append([0, 0])
                    self.renderedPolygons[i].pointList = pointList
        else:
            for i in range(len(dodecahedronFaces)):
                pointList = []
                for j in range(5):
                    pointList.append([0, 0])
                    self.renderedPolygons[i].pointList = pointList

class Icosahedron:
    # not drawing it,,,
    def __init__(self, position, radius, rotVect, color=Vector3(135, 206, 235)):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.radius = radius
        self.position = position
        self.vertices = [
            Vector3(0, 1, goldenRatio),
            Vector3(0, -1, goldenRatio),
            Vector3(0, 1, -goldenRatio),
            Vector3(0, -1, -goldenRatio),
            Vector3(1, goldenRatio, 0),
            Vector3(-1, goldenRatio, 0),
            Vector3(1, -goldenRatio, 0),
            Vector3(-1, -goldenRatio, 0),
            Vector3(goldenRatio, 0, 1),
            Vector3(goldenRatio, 0, -1),
            Vector3(-goldenRatio, 0, 1),
            Vector3(-goldenRatio, 0, -1)
        ]
        for i in range(len(self.vertices)):
            vertex = self.vertices[i]
            vertex.mul(radius)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3() + position
        self.renderedLineSegments = [Line(0, 0, 0, 0) for i in range(30)]
        self.renderedPolygons = [Polygon() for i in range(20)]
        self.debugLabels = [Label('undefined', 0, 0, size=20, opacity=0) for i in range(12)]
        self.rotationVector = rotVect
        for i in range(20):
            self.renderedPolygons[i].fill = rgb(randrange(0, 255), randrange(0, 255), randrange(0, 255))
        self.color = color
    
    def changeProperties(self, position, radius, rotVect, color=Vector3(135, 206, 235)):
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        self.radius = radius
        self.position = position
        self.vertices = [
            Vector3(0, 1, goldenRatio),
            Vector3(0, -1, goldenRatio),
            Vector3(0, 1, -goldenRatio),
            Vector3(0, -1, -goldenRatio),
            Vector3(1, goldenRatio, 0),
            Vector3(-1, goldenRatio, 0),
            Vector3(1, -goldenRatio, 0),
            Vector3(-1, -goldenRatio, 0),
            Vector3(goldenRatio, 0, 1),
            Vector3(goldenRatio, 0, -1),
            Vector3(-goldenRatio, 0, 1),
            Vector3(-goldenRatio, 0, -1)
        ]
        for i in range(len(self.vertices)):
            vertex = self.vertices[i]
            vertex.mul(radius/goldenRatio)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(position)
        self.rotationVector = rotVect
        self.color = color
    
    def rotate(self, rotVect):
        self.rotationVector += rotVect
        rotationMatrix = Matrix.RotationMatrix(rotVect)
        for i in range(len(self.vertices)):
            vertex = self.vertices[i].subtract(self.position)
            vertexVector = vertex.ToColumnVector()
            rotatedVertex = Matrix.Multiply(rotationMatrix, vertexVector)
            self.vertices[i] = rotatedVertex.ColumnVectorToVector3().add(self.position)
    
    def render(self, wireframe=True, faces=True, debug=False):
        TDPV = [] # short for 2D Projected Vertices
        for vertex in self.vertices:
            TDPV.append(vertex.To2DProjection(camPosition, camOrientation).mul(400).add(center))
        if debug:
            for i in range(len(self.vertices)):
                self.debugLabels[i].opacity = 100
                self.debugLabels[i].value = i
                self.debugLabels[i].centerX = TDPV[i].x
                self.debugLabels[i].centerY = TDPV[i].y
        else:
            for i in range(len(self.vertices)):
                self.debugLabels[i].opacity = 0
                self.debugLabels[i].centerX = 0
                self.debugLabels[i].centerY = 0
        if wireframe:
            for i in range(len(icosahedronEdges)):
                self.renderedLineSegments[i].x1 = TDPV[icosahedronEdges[i][0]].x
                self.renderedLineSegments[i].y1 = TDPV[icosahedronEdges[i][0]].y
                self.renderedLineSegments[i].x2 = TDPV[icosahedronEdges[i][1]].x
                self.renderedLineSegments[i].y2 = TDPV[icosahedronEdges[i][1]].y
        else:
            for i in range(len(icosahedronEdges)):
                self.renderedLineSegments[i].x1 = 0
                self.renderedLineSegments[i].y1 = 0
                self.renderedLineSegments[i].x2 = 0
                self.renderedLineSegments[i].y2 = 0
        if faces:
            for i in range(len(icosahedronFaces)):
                pointList = []
                faceNormal = Vector3(0, 0, 0)
                faceCenter = Vector3(0, 0, 0)
                closestVertex = self.vertices[icosahedronFaces[i][0]].clone()
                for k in range(3): # this section calculates the normal of the face (used for back-face culling), and the center (used for "lighting")
                    # research newell algorithm for the face normal
                    currentVtx = icosahedronFaces[i][k]
                    nextVtx = icosahedronFaces[i][(k+1)%3]
                    cvtxPos = self.vertices[currentVtx]
                    nvtxPos = self.vertices[nextVtx]
                    faceNormal.x += (cvtxPos.y - nvtxPos.y) * (cvtxPos.z + nvtxPos.z)
                    faceNormal.y += (cvtxPos.z - nvtxPos.z) * (cvtxPos.x + nvtxPos.x)
                    faceNormal.z += (cvtxPos.x - nvtxPos.x) * (cvtxPos.y + nvtxPos.y)
                    faceCenter.x += cvtxPos.x
                    faceCenter.y += cvtxPos.y
                    faceCenter.z += cvtxPos.z
                faceCenter.mul(1/3)
                # lots of math; basically had to teach myself quick linear algebra for a lot of this stuff
                # taking the dot product of the camera to vertex vector and the face normal vector can determine if a face's normal is facing away
                # if it is, that means that the face must be facing away from the camera, and then we can cull it
                # this saves processing power, and also assists in properly giving the full 3d effect.
                camToVtxVect = closestVertex.sub(camPosition)
                facingAway = (camToVtxVect.dotProduct(faceNormal) >= 0)
                if not facingAway:
                    for j in range(3):
                        pointList.append([TDPV[icosahedronFaces[i][j]].x, TDPV[icosahedronFaces[i][j]].y])
                    self.renderedPolygons[i].pointList = pointList
                    lightAmnt = Vector3.FindDistance(faceCenter, lightPos)/400
                    self.renderedPolygons[i].fill = rgb(
                        min(self.color.x/lightAmnt, 255),
                        min(self.color.y/lightAmnt, 255),
                        min(self.color.z/lightAmnt, 255)
                    )
                    #self.renderedPolygons[i].toBack()
                    #self.renderedPolygons[i].fill = rgb(randrange(0, 255), randrange(0, 255), randrange(0, 255))
                else:
                    for j in range(3):
                        pointList.append([0, 0])
                    self.renderedPolygons[i].pointList = pointList
        else:
            for i in range(len(icosahedronFaces)):
                pointList = []
                for j in range(3):
                    pointList.append([0, 0])
                self.renderedPolygons[i].pointList = pointList

polyhedra = [
    #RectPrism(Vector3(-50, -50, 350), Vector3(50, 50, 250), Vector3(0, 0, 0))
    
    # ----- uncomment these groups of lines for some demos ----- #
    #Dodecahedron(Vector3(0, 100, 400), 40, Vector3(0, 0, 0), Vector3(144, 238, 144)),
    #Icosahedron(Vector3(0, -100, 400), 40, Vector3(0, 0, 0)),
    
    Dodecahedron(Vector3(0, 0, 400), 80, Vector3(0, 0, 0), Vector3(144, 238, 144)),
    
    #Icosahedron(Vector3(0, 0, 400), 80, Vector3(0, 0, 0)),
]

def initalize():
    print("Press W to enable/disable wireframe\nPress F to enable/disable face rendering\nPress D to enable/disable debug mode")

initalize()

def onStep():
    for obj in polyhedra:
        obj.render(wireframe=appdata['globalWireframe'], faces=appdata['globalFaceRender'], debug=appdata['globalDebug'])
        #obj.render()
        obj.rotate(Vector3(1, 1, 1))
        obj.rotate(Vector3(0, appdata['drag']['dX'], -appdata['drag']['dY']))
    
    appdata['steps'] += 0.01*math.pi
    #camPosition.x += appdata['drag']['dX']
    #camPosition.y += appdata['drag']['dY']
    appdata['drag']['dX'] *= 0.95
    appdata['drag']['dY'] *= 0.95

def onMousePress(mx, my):
    appdata['drag']['pX'] = mx
    appdata['drag']['pY'] = my
    appdata['drag']['X'] = mx
    appdata['drag']['Y'] = my
    
def onMouseDrag(mx, my):
    appdata['drag']['X'] = mx
    appdata['drag']['Y'] = my
    appdata['drag']['dX'] = appdata['drag']['pX'] - appdata['drag']['X']
    appdata['drag']['dY'] = appdata['drag']['pY'] - appdata['drag']['Y']
    appdata['drag']['pX'] = mx
    appdata['drag']['pY'] = my
    
def onMouseRelease(mx, my):
    appdata['drag']['X'] = 0
    appdata['drag']['Y'] = 0
    #appdata['drag']['dX'] = 0
    #appdata['drag']['dY'] = 0
    appdata['drag']['pX'] = 0
    appdata['drag']['pY'] = 0

def onKeyPress(key):
    if key == 'w':
        appdata['globalWireframe'] = not appdata['globalWireframe']
    if key == 'f':
        appdata['globalFaceRender'] = not appdata['globalFaceRender']
    if key == 'd':
        appdata['globalDebug'] = not appdata['globalDebug']