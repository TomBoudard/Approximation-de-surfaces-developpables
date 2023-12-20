import openmesh as om
import numpy as np
from tools import *
import heapq
from gaussian_map import mean_curvature


def vector(mesh, vertexStart, vertexEnd):
    position_vertexStart = getVertexNewPosition(mesh, vertexStart)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])


def getDevelopability(mesh, vertex):
    ''' Return developability of a vertex from the mesh, function named
    g in the paper
    '''
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            vector1 = vector(mesh, vertex, vh1)
            vector2 = vector(mesh, vertex, vh2)
            angle = np.arccos(
                np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
            developability -= angle
        developability = developability**2 #vraiment utile?
        return developability
    else:
        return 0

def initMesh(mesh):
    ''' set two new properties to the vertices of our mesh: the movement scale along each vertex's unit normal n, initialized to 0
    and its developability
    '''
    for vertex in mesh.vertices():
        mesh.set_vertex_property('movement_scale', vertex, 0)
        developability = getDevelopability(mesh, vertex)
        mesh.set_vertex_property('developability', vertex, developability)
        #print(mesh.vertex_property('developability')[vertex.idx()])
        # print(mesh.vertex_property('movement_scale')[vertex.idx()])
    return mesh


def initDic(mesh):
    """
    Take a mesh and returns a dictionary keyed on the vertices' id, with the [g(· · · )]2 measure as values
    """
    vertices = {vertex.idx() : [mesh.vertex_property('developability')[vertex.idx()], vertex] for vertex in mesh.vertices()}
    #print(vertices) #debug
    return vertices

def getVertex(mesh, dictionary):
    ''' return the vertex with the maximum [g(· · · )]2'''
    # vertex_id = max(dictionary, key = lambda x :dictionary.get()[0]) #give the vertex id
    # vertex_id = max(dictionary, key=lambda x: dictionary[x][0])
    # return dictionary[vertex_id][1]
    max_value_vertex_id = None
    max_value = float('-inf')  # Initialize max_value as negative infinity

    for vertex_id, value in dictionary.items():
        if isinstance(value, list) and len(value) >= 1 and isinstance(value[0], (int, float)):
            # Check if the value is a list and the first element is a numeric type (e.g., int or float)
            if value[0] > max_value:
                max_value = value[0]
                max_value_vertex_id = vertex_id

    return dictionary[max_value_vertex_id][1]

def updateDevelopability(mesh, moved_vertex, dictionary):
    ''' Update the developability of the vertex and its neighbours in the dictionary according to the movement scale'''
    vertices_to_update = [vh for vh in mesh.vv(moved_vertex)]
    vertices_to_update.append(moved_vertex)
    for vertex in vertices_to_update:
        developability = getDevelopability(mesh, vertex)
        mesh.set_vertex_property('developability', vertex, developability)
        dictionary[vertex.idx()] = developability

def localOptimizationConstraint(mesh, vertex, delta):
    '''return the constraint T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2 of local optimization'''
    normal = mesh.normal(vertex)
    new_position = mesh.point(vertex) + delta * normal
    mesh.set_point(vertex, new_position)
    print("new_vertex =", new_position)
    new_vertex_developability = getDevelopability(mesh, vertex) # g(q + δn q )
    neighbours = [vh for vh in mesh.vv(vertex)] #les qj
    sum = 0
    for v in neighbours:
        sum += getDevelopability(mesh, v)**2
    constraint = new_vertex_developability**2 + sum
    return constraint

def updateVertex(mesh, vertex):
    '''
    Compute and update the movement scale of the given vertex in the mesh
    δ = δ0 - T(δ0)/dT(δ0)
    '''
    # Compute the movement scale
    prev_movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()] # delta_0
    constraint = localOptimizationConstraint(mesh, vertex, prev_movement_scale) # T(delta_0)
    derivative_constraint = central_difference(mesh, localOptimizationConstraint, vertex, prev_movement_scale)# dT(delta_0)
    new_movement_scale = prev_movement_scale - constraint/derivative_constraint
    mesh.set_vertex_property('movement_scale', vertex, new_movement_scale) #update the movement scale
    return 0

def central_difference(mesh, T, vertex, delta, h=1e-5):
    '''
    Use central difference to approximat dT
     f(x + h) - f(x - h)) / (2 * h)
    '''
    return (T(mesh, vertex, delta + h) - T(mesh, vertex, delta - h)) / (2 * h)



def updateNeighbour(mesh, vertex):
    '''
    Update the developabilty of the neighbourhood of the vertex we moved
    '''
    return 0

def updateVerticesPositions(mesh):
    '''
    Update the positions of all the vertices by their movement scales
    '''
    for vertex in mesh.vertices():
        new_position = getVertexNewPosition(mesh, vertex)
        mesh.set_point(vertex, new_position)

def getVertexNewPosition(mesh, vertex):
    '''
    Return the new position of a vertex according to its movement scale
    '''
    current_position = mesh.point(vertex)
    movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()]
    if movement_scale == None: #A cause de l'initialisation pour pas parcourir deux fois le mesh dans initMesh
        movement_scale = 0
    normal = mesh.normal(vertex)
    # print("current_position = ", current_position)
    # print("movement_scale = ", movement_scale)
    # print("normal = ", normal)
    new_position = current_position + movement_scale * normal
    return new_position

def main():
    print("---------- Début main\n")
    maxIter = 100
    nbIter = 0
    epsilon = 0.001

    ###    Read .off file
    filename = "../Objects/EmpireDress.off"
    mesh = om.read_trimesh(filename)
    a, b = add_angles(mesh)
    initial_object = mean_curvature(mesh, a, b)
    om.write_mesh("obj_initial.off", initial_object, vertex_color = True)

    ###    Compute the vertex developability detect function g(q) at each vertex on the given mesh patch and set the movement scale to each vertex to 0
    initMesh(mesh)

    ###    Compute the unit normal n of each vertex
    #NONNNNN A FAIRE
    mesh.update_vertex_normals()
    # print(mesh.vertex_normals()) #pour debug

    ###    Place all vertices in a dictionary H keyed on the vertices' id and with the [g(· · · )]2 measure as values and get the vertex with max values
    vertices_dic = initDic(mesh)
    vertex = getVertex(mesh, vertices_dic)
    max_developability = mesh.vertex_property('developability')[vertex.idx()]
    #while (the developability of the vertex with max developability is greater than ε) and (nbIter < maxIter);
    while ((max_developability > epsilon) and (nbIter < maxIter) ):
        ###    Update worst_vertex's movement scale along its unit normal n according to Eq. 17;
        updateVertex(mesh, vertex)

        ###    Update the cost of q and its adjacent vertices to reflect the movement on q
        #updateNeighbour(...)

        ###    Update the dictionary
        updateDevelopability(mesh, vertex, vertices_dic)

        #Select the new vertex for nest iteration
        vertex = getVertex(mesh, vertices_dic)
        max_developability = mesh.vertex_property('developability')[vertex.idx()]
        nbIter += 1
    print("---------- Fin algo\n")
    print("max_developability = ", max_developability)
    print("nbIter = ", nbIter)
    # Update the positions of all the vertices by their movement scales
    updateVerticesPositions(mesh)
    c, d = add_angles(mesh)
    optimize_object = mean_curvature(mesh, c, d)
    om.write_mesh("obj_optimized.off", optimize_object, vertex_color = True)
    # Update the normal vectors of all the vertices on O
    #Utile pour nous?

def main2():  # Dummy main with basic mesh to test the functions
    print("---------- Début main\n")
    mesh = om.TriMesh() #créer un mesh
    #initialize a mesh with 4 vertices
    vh0 = mesh.add_vertex([0, 1, 0])
    vh1 = mesh.add_vertex([1, 0, 0])
    vh2 = mesh.add_vertex([2, 1, 0])
    vh3 = mesh.add_vertex([0,-1, 0])
    # vh4 = mesh.add_vertex([2,-1, 0])

    # add a couple of faces to the mesh
    fh0 = mesh.add_face(vh0, vh1, vh2)
    #fh1 = mesh.add_face(vh1, vh3, vh4)
    fh2 = mesh.add_face(vh0, vh3, vh1)
    #vh_list = [vh2, vh1, vh4]
    #fh3 = mesh.add_face(vh_list)
    # TODO implémenter algorithme
    # a, b = add_angles(mesh)
    # mesh.request_face_colors() # pour ajouter des couleurs aux faces?
    # new_mesh = mean_curvature(mesh, a, b)
    # mesh.get_color
    maxIter = 1
    nbIter = 0
    epsilon = 0.001


    initMesh(mesh)
    ###    Compute the unit normal n of each vertex
    mesh.update_vertex_normals()
    #print(mesh.vertex_normals()) #pour debug

    ###    Place all vertices in a dictionary H keyed on the vertices' id and with the [g(· · · )]2 measure as values and get the vertex with max values
    vertices_dic = initDic(mesh)
    vertex = getVertex(mesh, vertices_dic)
    max_developability = mesh.vertex_property('developability')[vertex.idx()]
    #while (the developability of the vertex with max developability is greater than ε) and (nbIter < maxIter);
    while ((max_developability > epsilon) and (nbIter < maxIter) ):
        vertex = getVertex(mesh, vertices_dic)
        ###    Update worst_vertex's movement scale along its unit normal n according to Eq. 17;
        updateVertex(mesh, vertex)
        ###    Update the cost of q and its adjacent vertices to reflect the movement on q
        #updateNeighbour(...)
        ###    Update the dictionary
        updateDevelopability(mesh, vertex, vertices_dic)
        nbIter += 1
    for vertex in mesh.vertices():
        print(mesh.vertex_property('movement_scale')[vertex.idx()])
    om.write_mesh("test_dummy_main.off", mesh) # write .off file
    #mean_curvature(mesh)
    print("---------- Fin main\n")

if __name__ == "__main__":
    main()