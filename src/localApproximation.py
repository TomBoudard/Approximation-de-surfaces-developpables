import openmesh as om
import numpy as np
from tools import *
import heapq
from gaussian_map import mean_curvature

def vector(mesh, vertexStart, vertexEnd):
    return (mesh.point(vertexEnd)[0] - mesh.point(vertexStart)[0],
            mesh.point(vertexEnd)[1] - mesh.point(vertexStart)[1],
            mesh.point(vertexEnd)[2] - mesh.point(vertexStart)[2])

def vectorNew(mesh, vertexStart, vertexEnd, h=0):
    position_vertexStart = getVertexNewPosition(mesh, vertexStart, h)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd, h)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])

def vectorNewMain(mesh, vertexStart, vertexEnd, h=0):
    position_vertexStart = getVertexNewPosition(mesh, vertexStart, h)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd, 0)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])

def vectorNewNeighbour(mesh, vertexStart, vertexEnd, h=0):
    position_vertexStart = getVertexNewPosition(mesh, vertexStart, 0)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd, h)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])


def developabilityDetectFunction(mesh, minDevelopability=0, maxDevelopabilty=0):
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

    if minDevelopability == 0:
        minSum = min(listSum)
    else:
        minSum = minDevelopability

    for vertex in mesh.vertices():
        if listSum[vertex.idx()]:
            normalizedSum = (listSum[vertex.idx()] - minSum) / (maxSum - minSum)
            r, g, b = gradientColor(normalizedSum)
            color = [r, g, b, 1.]
        else:
            color = [0., 0., 0., 1.]
        mesh.set_color(vertex, color)
    return minSum, maxSum

def getDevelopability(mesh, vertex, vectorFunction, h=0):
    ''' Return developability of a vertex from the mesh, function named g in the paper '''
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            vector1 = vectorFunction(mesh, vertex, vh1, h=h)
            vector2 = vectorFunction(mesh, vertex, vh2, h=h)
            vector1 /= np.linalg.norm(vector1)
            vector2 /= np.linalg.norm(vector2)
            angle = np.arccos(np.dot(vector1, vector2))
            developability -= angle
        # developability = developability**2 #vraiment utile?
        return developability
    else:
        return 0

def getDevelopabilityNeighbour(mesh, vertex, mainVertex, h=0):
    ''' Return developability of a vertex from the mesh, function named g in the paper '''
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            if vh1 == mainVertex:
                vector1 = vectorNewNeighbour(mesh, vertex, vh1, h=h)
            else:
                vector1 = vectorNewNeighbour(mesh, vertex, vh1, h=0)
            if vh2 == mainVertex:
                vector2 = vectorNewNeighbour(mesh, vertex, vh2, h=h)
            else:
                vector2 = vectorNewNeighbour(mesh, vertex, vh2, h=0)
            vector1 /= np.linalg.norm(vector1)
            vector2 /= np.linalg.norm(vector2)
            angle = np.arccos(np.dot(vector1, vector2))
            developability -= angle
        # developability = developability**2 #vraiment utile?
        return developability
    else:
        return 0

def initMesh(mesh):
    ''' set two new properties to the vertices of our mesh: the movement scale along each vertex's unit normal n, initialized to 0
    and its developability
    '''
    for vertex in mesh.vertices():
        mesh.set_vertex_property('movement_scale', vertex, 0)
        developability = getDevelopability(mesh, vertex, vectorNew)
        mesh.set_vertex_property('developability', vertex, developability)
        # print("Movement scale ", mesh.vertex_property('movement_scale')[vertex.idx()])
        # normal = getNormal(mesh, vertex) #FIXME
        # mesh.set_vertex_property('normal', vertex, normal)
        # print("Normal ", mesh.vertex_property('normal')[vertex.idx()])
    mesh.update_vertex_normals()
    # for vertex in mesh.vertices():
    #     pass
    # print("Developability ", mesh.vertex_property('developability')[vertex.idx()])
    return mesh


def initDic(mesh):
    """
    Take a mesh and returns a dictionary keyed on the vertices' id, with the [g(· · · )]2 measure as values
    """
    vertices = {vertex.idx() : [mesh.vertex_property('developability')[vertex.idx()]**2, vertex] for vertex in mesh.vertices()}
    # print(vertices) #debug
    return vertices

def getVertex(dictionary):
    ''' return the vertex with the maximum [g(· · · )]2'''
    max_value_vertex_id = None
    max_value = -1  # Initialize max_value as negative infinity

    # print("DICTIONNAIRE: ", dictionary.items())
    for vertex_id, value in dictionary.items():
        # Check if the value is a list and the first element is a numeric type (e.g., int or float)
        if value[0] > max_value:
            max_value = value[0]
            max_value_vertex_id = vertex_id

    return dictionary[max_value_vertex_id][1]

def updateDevelopabilityNeighbour(mesh, moved_vertex, dictionary):
    ''' Update the developability of the vertex and its neighbours in the dictionary according to the movement scale'''
    vertices_to_update = [vh for vh in mesh.vv(moved_vertex)]
    for vertex in vertices_to_update:
        developability = getDevelopability(mesh, vertex, vectorNew)
        mesh.set_vertex_property('developability', vertex, developability)
        if vertex.idx() in dictionary.keys():
            dictionary[vertex.idx()][0] = developability


def localOptimizationConstraint(mesh, vertex, h=0):
    '''return the constraint T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2 of local optimization'''
    new_vertex_developability = getDevelopability(mesh, vertex, vectorNewMain, h=h) # g(q + δn q )
    neighbours = [vh for vh in mesh.vv(vertex)] #les qj
    sum = 0
    for v in neighbours:
        sum += getDevelopabilityNeighbour(mesh, v, vertex, h=h)**2 #FIXME Use h for Vertex
    constraint = new_vertex_developability**2 + sum

    return constraint

def updateVertex(mesh, vertex, dictionary, maxMovementScale):
    '''
    Compute and update the movement scale of the given vertex in the mesh
    δ = δ0 - T(δ0)/dT(δ0)
    '''
    # Compute the movement scale
    prev_movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()] # delta_0
    constraint = localOptimizationConstraint(mesh, vertex, prev_movement_scale) # T(delta_0)
    derivative_constraint = central_difference(mesh, localOptimizationConstraint, vertex)# dT(delta_0)
    # prev_developability = getDevelopability(mesh, vertex, vectorNew)
    if (derivative_constraint != 0):
        new_movement_scale = prev_movement_scale - constraint/derivative_constraint
        if new_movement_scale > maxMovementScale:
            new_movement_scale = maxMovementScale
        elif new_movement_scale < -maxMovementScale:
            new_movement_scale = -maxMovementScale
        mesh.set_vertex_property('movement_scale', vertex, new_movement_scale) #update the movement scale
        new_developability = getDevelopability(mesh, vertex, vectorNew)
        # if (abs(new_developability) > abs(prev_developability)):
        #     mesh.set_vertex_property('movement_scale', vertex, -new_movement_scale)  # update the movement scale
        #     new_developability = getDevelopability(mesh, vertex, True)
            # if (abs(new_developability) > abs(prev_developability)):
            #     mesh.set_vertex_property('movement_scale', vertex, prev_movement_scale)  # do not update the movement scale
            #     # new_developability = 0 #Not true but the vertex can't be more optimised
        mesh.set_vertex_property('developability', vertex, new_developability)
        dictionary[vertex.idx()][0] = new_developability
        return 1

    else:
        mesh.set_vertex_property('movement_scale', vertex, 0) #update the movement scale
        return 0


def getNormal(mesh, vertex):
    normal = np.array([.0, .0, .0])
    for face in mesh.vf(vertex):
        faceVertex = []
        for vertex in mesh.fv(face):
            point = mesh.point(vertex)
            faceVertex.append(np.array([point[0], point[1], point[2]]))

        firstVector = faceVertex[1]-faceVertex[0]
        secondVector = faceVertex[2]-faceVertex[0]
        normal = normal + np.cross(firstVector, secondVector)

    return normal/np.linalg.norm(normal)

def updateNormals(mesh):
    # for vertex in mesh.vertices():
    #     normal = getNormal(mesh, vertex)
    #     mesh.set_vertex_property('normal', vertex, normal)
    mesh.update_vertex_normals()

    return 0

def central_difference(mesh, T, vertex, h=1e-5):
    '''
    Use central difference to approximate dT
     f(x + h) - f(x - h)) / (2 * h)
    '''
    a = T(mesh, vertex, h=h)
    b = T(mesh, vertex, h=-h)
    c = a-b
    return c/(2*h)

def updateVerticesPositions(mesh):
    '''
    Update the positions of all the vertices by their movement scales
    '''
    for vertex in mesh.vertices():
        new_position = getVertexNewPosition(mesh, vertex)
        mesh.set_point(vertex, new_position)


def getVertexNewPosition(mesh, vertex, h=0):
    '''
    Return the new position of a vertex according to its movement scale
    '''
    current_position = mesh.point(vertex)
    movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()]
    if movement_scale == None: #A cause de l'initialisation pour pas parcourir deux fois le mesh dans initMesh
        movement_scale = 0
    # normal = mesh.vertex_property('normal')[vertex.idx()]
    # new_position = np.array([current_position[0], current_position[1], current_position[2]])
    # new_position[0] += (movement_scale + h) * normal[0]
    # new_position[1] += (movement_scale + h) * normal[1]
    # new_position[2] += (movement_scale + h) * normal[2]

    normal = mesh.normal(vertex)
    new_position = current_position + (movement_scale + h) * normal

    return new_position

def main():
    print("---------- Début main\n")
    maxIter = 200
    nbIter = 0
    epsilon = 0.1
    maxMovementScale = 0.05

    ###    Read .off file
    filename = "../Objects/eight.off"
    mesh = om.read_trimesh(filename)
    # a, b = add_angles(mesh)
    # initial_object = mean_curvature(mesh, a, b)
    initialMinDevelopability, initialMaxDevelopability = developabilityDetectFunction(mesh)
    om.write_mesh("obj_initial.off", mesh, vertex_color = True)

    #    Compute the vertex developability detect function g(q) at each vertex on the given mesh patch and set the movement scale to each vertex to 0
    initMesh(mesh)

    #    Place all vertices in a dictionary H keyed on the vertices' id and with the [g(· · · )]2 measure as values and get the vertex with max values
    vertices_dic = initDic(mesh)
    vertex = getVertex(vertices_dic)
    max_developability = mesh.vertex_property('developability')[vertex.idx()]
    print("START | Max index: ", vertex.idx(), "Max developability = ", max_developability)

    #while (the developability of the vertex with max developability is greater than ε) and (nbIter < maxIter);
    while ((abs(max_developability) > epsilon) and (nbIter < maxIter) ):
        #    Update worst_vertex's movement scale along its unit normal n according to Eq. 17;
        res = updateVertex(mesh, vertex, vertices_dic, maxMovementScale)
        if (res == 0):
            vertices_dic.pop(vertex.idx())

        else:
            #    Update the cost of q and its adjacent vertices to reflect the movement on q
            #    Update the dictionary
            updateDevelopabilityNeighbour(mesh, vertex, vertices_dic)

        #    Select the new vertex for next iteration
        vertex = getVertex(vertices_dic)
        max_developability = abs(mesh.vertex_property('developability')[vertex.idx()])
        print("Max index: ", vertex.idx(), "Max developability = ", max_developability)
        nbIter += 1
    print("---------- Fin algo\n")
    # Update the positions of all the vertices by their movement scales
    updateVerticesPositions(mesh)
    updateNormals(mesh)
    # c, d = add_angles(mesh)
    # optimize_object = mean_curvature(mesh, c, d)
    developabilityDetectFunction(mesh, initialMinDevelopability, initialMaxDevelopability)
    om.write_mesh("obj_optimized.off", mesh, vertex_color = True)
    # Update the normal vectors of all the vertices on O
    # print("max_developability = ", max_developability)
    print("nbIter = ", nbIter)

if __name__ == "__main__":
    main()