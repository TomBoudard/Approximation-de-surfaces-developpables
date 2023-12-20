import openmesh as om

def developabilityDetectFunction(mesh):
    ...
    return 0

#Thaïs
def computeMeshNormals(mesh):
    ...
    return 0

#Thaïs
def initHeap(mesh):
    """
    Take a mesh and returns a heap of the vertices according to their developability
    """
    return 0

def updateHeap(mesh, vertex):
    ...
    return 0


def updateVertex(mesh, vertex):
    """
    Take a mesh and a vertex and move this vertex along its normal
    """

    coordinates = mesh.point(vertex)
    normal = mesh.vertex_normal(vertex)

    delta = 1

    newCoordinates = coordinates + delta*normal

    mesh.set_point(vertex, newCoordinates)


def updateNeighbour(mesh, vertex):
    ...     
    return 0


def main():
    maxIter = 20
    nbIter = 0
    filename = "../Objects/sphere.off"
    mesh = om.read_trimesh(filename)

    developabilityDetectFunction(mesh)
    computeMeshNormals(mesh)
    heapVertex = initHeap(mesh)

    while(...):
        updateVertex(...)
        updateNeighbour(...)

        heapVertex = initHeap(mesh)

        nbIter += 1

    return 0

