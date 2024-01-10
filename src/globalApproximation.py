import openmesh as om
import numpy as np

from tools import *


def surfaceDiscrepancy(originalMesh, editedMesh):
    discrepancy = 0
    for originalEdge, editedEdge in zip(originalMesh.edges(), editedMesh.edges()):
        originalHalfEdge = originalMesh.halfedge_handle(originalEdge, 0)
        editedHalfEdge = editedMesh.halfedge_handle(editedEdge, 0)

        originalVertexHandleFrom = originalMesh.from_vertex_handle(originalHalfEdge)
        originalVertexHandleTo = originalMesh.to_vertex_handle(originalHalfEdge)

        editedVertexHandleFrom = editedMesh.from_vertex_handle(editedHalfEdge)
        editedVertexHandleTo = editedMesh.to_vertex_handle(editedHalfEdge)

        editedNorm = np.linalg.norm(vector(editedMesh, editedVertexHandleFrom, editedVertexHandleTo))
        originalNorm = np.linalg.norm(vector(originalMesh, originalVertexHandleFrom, originalVertexHandleTo))

        discrepancy += (editedNorm - originalNorm) * (editedNorm - originalNorm)

    return discrepancy


def developabilityDetectFunction(mesh, maxDevelopabilty=0):
    developability = 0
    listSum = []
    for vertex in mesh.vertices():
        if not mesh.is_boundary(vertex):
            sum = 2 * np.pi
            neighbours = [vh for vh in mesh.vv(vertex)]
            for i in range(len(neighbours)):
                vh1 = neighbours[i]
                vh2 = neighbours[(i + 1) % len(neighbours)]
                vector1 = vector(mesh, vertex, vh1)
                vector2 = vector(mesh, vertex, vh2)
                angle = np.arccos(
                    np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
                sum -= angle
            listSum.append(sum)
            developability += sum
        else:
            listSum.append(False)
    if maxDevelopabilty == 0:
        maxSum = max(listSum)
    else:
        maxSum = maxDevelopabilty
    minSum = min(listSum)
    for vertex in mesh.vertices():
        if listSum[vertex.idx()]:
            normalizedSum = (listSum[vertex.idx()] - minSum) / (maxSum - minSum)
            r, g, b = gradientColor(normalizedSum)
            color = [r, g, b, 1.]
        else:
            color = [0., 0., 0., 1.]
        mesh.set_color(vertex, color)
    return maxSum

def rhoInitalization(mesh):
    nbEdge = 0
    sumEdgeLength = 0
    developability = developabilityDetectFunction(mesh)

    for edge in mesh.edges():
        nbEdge += 1

        halfEdge = mesh.halfedge_handle(edge, 0)

        vertexHandleFrom = mesh.from_vertex_handle(halfEdge)
        vertexHandleTo = mesh.to_vertex_handle(halfEdge)

        edgeNorm = np.linalg.norm(vector(mesh, vertexHandleFrom, vertexHandleTo))

        sumEdgeLength += edgeNorm * edgeNorm

    return sumEdgeLength/(nbEdge*developability)



def main():
    filename = "../Objects/sphere.off"
    nbMaxIteration = 10

    originalMesh = om.read_trimesh(filename)
    editedMesh = om.read_trimesh(filename)

    rho = rhoInitalization(originalMesh)

    developability = developabilityDetectFunction(editedMesh)

    surfaceDiscrepancy(originalMesh, editedMesh)

    om.write_mesh("../output/test.off", editedMesh, vertex_color=True)

    return 0


if __name__ == "__main__":
    main()
