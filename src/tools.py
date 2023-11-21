import openmesh as om
import numpy as np

def vector(mesh, vertexStart, vertexEnd):
    return (mesh.point(vertexEnd)[0] - mesh.point(vertexStart)[0],
            mesh.point(vertexEnd)[1] - mesh.point(vertexStart)[1],
            mesh.point(vertexEnd)[2] - mesh.point(vertexStart)[2])