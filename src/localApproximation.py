import openmesh as om
import numpy as np
from tools import *
import heapq
from gaussian_map import mean_curvature
from globalApproximation import developabilityDetectFunction


def vector(mesh, vertexStart, vertexEnd):
    # position_vertexStart = mesh.point(vertexStart)
    # position_vertexEnd = mesh.point(vertexEnd)  #Tenter de changer ça
    position_vertexStart = getVertexNewPosition(mesh, vertexStart)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])

def new_vector(mesh, vertexStart, vertexEnd, delta):
    normal = mesh.normal(vertexStart)
    old_position = mesh.point(vertexStart)
    position_vertexStart = old_position + delta * normal
    # position_vertexStart = getVertexNewPosition(mesh, vertexStart)
    position_vertexEnd =  getVertexNewPosition(mesh, vertexEnd)
    # position_vertexEnd = mesh.point(vertexEnd) #############               tenter changer ça

    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])

def new_vector2(mesh, vertexStart, vertexEnd, delta):
    normal = mesh.normal(vertexEnd)
    old_position = mesh.point(vertexEnd)
    position_vertexEnd = old_position + delta * normal
    # position_vertexStart = getVertexNewPosition(mesh, vertexStart)
    position_vertexStart =  getVertexNewPosition(mesh, vertexStart)
    # position_vertexStart = mesh.point(vertexStart) #############               tenter changer ça

    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])



def getDevelopability(mesh, vertex): # Utilise les positions de bases!!
    ''' Return developability of a vertex from the mesh, function named
    g in the paper
    '''
    # print("Entrée dans local getDevelopability")

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
        developability = developability
        return developability
    else:
        return 0


def getNewDevelopability(mesh, vertex, delta): # Utilise les positions de bases!!
    ''' Return developability of a vertex from the mesh, function named
    g in the paper
    '''
    # print("Entrée dans local getewDevelopability")
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            vector1 = new_vector(mesh, vertex, vh1, delta)
            vector2 = new_vector(mesh, vertex, vh2, delta)
            angle = np.arccos(
                np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
            developability -= angle
        developability = developability
        return developability
    else:
        return 0


def getNewDevelopabilityNeighboors(mesh, neighbour , moved_vertex, delta):
    # print("Entrée dans local getNewDevelopabilityNeighboors")
    if not mesh.is_boundary(neighbour): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(neighbour)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            if (vh1 == moved_vertex): #Le voisin est notre vertex qui a bougé
                #On prend en compte le mouvement scale de celui ci dans calcul de developpabilité de notre vertex
                # print("Le voisin 1 est notre vertex qui a bougé")
                vector1 = new_vector2(mesh, neighbour, vh1, delta)
            else:
                vector1 = vector(mesh, neighbour, vh1)
            if (vh2 == moved_vertex): #Le voisin est notre vertex qui a bougé
                #On prend en compte le mouvement scale de celui ci dans calcul de developpabilité de notre vertex
                # print("Le voisin 2 est notre vertex qui a bougé")
                vector2 = new_vector2(mesh, neighbour, vh2, delta)
            else:
                vector2 = vector(mesh, neighbour, vh2)
            # print("vector1 = ", vector1)
            # print("vector2 = ", vector2)

            angle = np.arccos(
                np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
            developability -= angle
        developability = developability
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
    return mesh


def initDic(mesh):
    """
    Take a mesh and returns a dictionary keyed on the vertices' id, with the [g(· · · )]2 measure as values 
    """
    vertices = {vertex.idx() : [(mesh.vertex_property('developability')[vertex.idx()])**2, vertex] for vertex in mesh.vertices()}
    # print(vertices) #debug
    return vertices

def getVertex(mesh, dictionary):
    ''' return the vertex with the maximum [g(· · · )]2'''
    # print("Entrée dans la fonction getVertex")

    max_value_vertex_id = None
    max_value = -1  # Initialize max_value as negative infinity

    for vertex_id, value in dictionary.items():
        # Check if the value is a list and the first element is a numeric type (e.g., int or float)
        if value[0] > max_value:
            max_value = value[0]
            max_value_vertex_id = vertex_id
    return dictionary[max_value_vertex_id][1]

def updateDevelopability(mesh, moved_vertex, dictionary):
    ''' Update the developability of the vertex and its neighbours in the dictionary according to the movement scale'''
    vertices_to_update = [vh for vh in mesh.vv(moved_vertex)]
    vertices_to_update.append(moved_vertex)
    # old_position = mesh.point(moved_vertex)
    # normal = mesh.normal(moved_vertex)
    # movement_scale = mesh.vertex_property('movement_scale')[moved_vertex.idx()]
    # new_position = old_position + movement_scale * normal
    # mesh.set_point(moved_vertex, new_position)
    for vertex in vertices_to_update:
        developability = getDevelopability(mesh, vertex) #Sinon remplacer par getNewDevelopability
        mesh.set_vertex_property('developability', vertex, developability)
        dictionary[vertex.idx()][0] = developability**2
    # mesh.set_point(moved_vertex, old_position)
    return dictionary

def localOptimizationConstraint(mesh, moved_vertex, delta):
    '''return the constraint T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2 of local optimization'''
    # print("Entrée dans local optimization")
    # normal = mesh.normal(vertex)
    # old_position = mesh.point(vertex)
    # new_position = mesh.point(vertex) + delta * normal
    # mesh.set_point(vertex, new_position)
    # prev_movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()]

    #Calcule dévelopabilité de q + δn q
    new_vertex_developability = getNewDevelopability(mesh, moved_vertex, delta) # g(q + δn q )

    # print("new_vertex_developability =", new_vertex_developability)
    neighbours = [vh for vh in mesh.vv(moved_vertex)] #les qj
    sum = 0
    for v in neighbours:
        sum += getNewDevelopabilityNeighboors(mesh, v , moved_vertex, delta)**2 #Calculer la nouvelle developabilité des voisins en prenant en compte movement scale du vertex
        # print("getNewDevelopability(mesh, v)**2 = ", getNewDevelopability(mesh, v)**2)
    # print("sum =", new_vertex_developability)
    # print("sum =", sum)
    constraint = new_vertex_developability**2 + sum
    # mesh.set_point(vertex, old_position)
    # print("Constraint = ", constraint)
    return constraint

def updateVertex(mesh, vertex):
    '''
    Compute and update the movement scale of the given vertex in the mesh
    δ = δ0 - T(δ0)/dT(δ0)
    return 0 if the vertex need to be fixed
    else return 1
    '''
    print("Entrée dans updateVertex")
    # Compute the movement scale
    prev_movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()] # delta_0
    # print("prev_movement_scale = ", prev_movement_scale)
    constraint = localOptimizationConstraint(mesh, vertex, prev_movement_scale) # T(delta_0)
    derivative_constraint = DerivativelocalOptimizationConstraint(mesh, vertex, prev_movement_scale)# dT(delta_0)
    if (derivative_constraint != 0):
        new_movement_scale = prev_movement_scale - constraint/derivative_constraint
        mesh.set_vertex_property('movement_scale', vertex, new_movement_scale) #update the movement scale
    else:
        mesh.set_vertex_property('movement_scale', vertex, 0) #update the movement scale
        print( "!!!!!! derivative_constraint = 0")
        # We remove the vertices from the dict so we don't iterate on it because we fix it
        return 0
    # print("prev_movement_scale =", prev_movement_scale)
    # print("constraint =", constraint)
    # print("derivative_constraint =", derivative_constraint)
    # print("new_movement_scale =", new_movement_scale)

    return 1

def DerivativelocalOptimizationConstraint(mesh, moved_vertex, delta):
    '''return the constraint derivative constrainte (where the constraint is T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2)
    Use central difference to approximat dT
     f(x + h) - f(x - h)) / (2 * h)'''
    h=1e-5

    term1 = localOptimizationConstraint(mesh, moved_vertex, delta + h)
    term2 = localOptimizationConstraint(mesh, moved_vertex, delta - h)
    #  #f(x+h)
    # # normal = mesh.normal(vertex)
    # # old_position = mesh.point(vertex)
    # # new_position = mesh.point(vertex) + (delta +h) * normal
    # # mesh.set_point(vertex, new_position)
    # new_vertex_developability1 = getNewDevelopability(mesh, vertex, delta + h) # g(q + δn q )
    # # print("new_vertex_developability =", new_vertex_developability1)
    # neighbours = [vh for vh in mesh.vv(vertex)] #les qj
    # sum = 0
    # for v in neighbours:
    #     sum += getDevelopability(mesh, v)**2
    #     # print("getDevelopability(mesh, v)**2 = ", getDevelopability(mesh, v)**2)
    # # print("sum =", new_vertex_developability)
    # # print("sum =", sum)
    # term1 = new_vertex_developability1**2 + sum
    # # new_position2 = mesh.point(vertex) + (delta - h) * normal
    # # mesh.set_point(vertex, new_position2)
    # new_vertex_developability2 = getNewDevelopability(mesh, vertex , delta - h) # g(q + δn q )
    # sum = 0
    # for v in neighbours:
    #     sum += getDevelopability(mesh, v)**2
    # # print("sum =", new_vertex_developability)
    # # print("sum =", sum)
    # term2 = new_vertex_developability2**2 + sum
    derivative_constraint = (term1 - term2)/(2*h)
    # mesh.set_point(vertex, old_position)
    return derivative_constraint

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
    new_position = current_position + movement_scale * normal
    return new_position

def main():
    print("---------- Début main\n")
    maxIter = 200
    nbIter = 0
    epsilon = 0.01
    print("Nombre d'itérations max: ", maxIter)
    print("epsilon = ", epsilon)

    ###    Read .off file
    filename = "../Objects/mesh_00040.off"
    print("lecture du fichier: ", filename)
    mesh = om.read_trimesh(filename)
    # a, b = add_angles(mesh)
    # initial_object = mean_curvature(mesh, a, b)
    developabilityDetectFunction(mesh)
    om.write_mesh("obj_initial.off", mesh, vertex_color = True)
    print("génération de la couleur des courbures de Gauss sur l'objet initial")

    ###    Compute the vertex developability detect function g(q) at each vertex on the given mesh patch and set the movement scale to each vertex to 0
    initMesh(mesh)

    ###    Compute the unit normal n of each vertex
    mesh.update_vertex_normals()
    # print("Calcul des normals:\n")
    # print(mesh.vertex_normals()) #pour debug

    ###    Place all vertices in a dictionary H keyed on the vertices' id and with the [g(· · · )]2 measure as values and get the vertex with max values
    vertices_dic = initDic(mesh)
    vertex = getVertex(mesh, vertices_dic)
    max_developability = mesh.vertex_property('developability')[vertex.idx()]
    print("L'id du vertex initiale à optimiser:", vertex.idx())
    print("Sa developabilité: (=max_dev) ", max_developability)
    initial_max_developability = max_developability
    #while (the developability of the vertex with max developability is greater than ε) and (nbIter < maxIter);
    print("---------- Début algo\n")
    while ((abs(max_developability) > epsilon) and (nbIter < maxIter) ):
        print("Iteration numéro: ", nbIter)
        print("L'id du vertex à optimiser:", vertex.idx())
        print("Sa developabilité 1: (=max_dev) ", max_developability)
        # print("Sa position: ", mesh.point(vertex))
        movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()] # delta_0
        # print("son movement scale =", movement_scale)
        ###    Update worst_vertex's movement scale along its unit normal n according to Eq. 17;
        res = updateVertex(mesh, vertex)
        max_developability = mesh.vertex_property('developability')[vertex.idx()]
        # print("Sa developabilité 1bis: (=max_dev) ", max_developability)
        if (res == 0):
            #this vertex need to be fixed so we remove it from the dic:
            vertices_dic.pop(vertex.idx())
        else:  
            ###    Update the cost of q and its adjacent vertices to reflect the movement on q
            #updateNeighbour(...)
            ###    Update the dictionary
            vertices_dic = updateDevelopability(mesh, vertex, vertices_dic)
            max_developability = mesh.vertex_property('developability')[vertex.idx()]
            # print("Sa developabilité 2: (=max_dev) ", max_developability)
        vertex = getVertex(mesh, vertices_dic)

        max_developability = mesh.vertex_property('developability')[vertex.idx()]
        # print("Sa developabilité 3: (=max_dev) ", max_developability)

        max_developability = mesh.vertex_property('developability')[vertex.idx()]
        # print("Sa developabilité 4: (=max_dev) ", max_developability)
        nbIter += 1
    print("---------- Fin algo\n")
    print("max_developability initiale = ", initial_max_developability)
    print("max_developability finale = ", max_developability)
    print("nbIter final = ", nbIter)
    # Update the positions of all the vertices by their movement scales
    updateVerticesPositions(mesh)
    # c, d = add_angles(mesh)
    # optimize_object = mean_curvature(mesh, c, d)
    # developabilityDetectFunction(mesh)
    om.write_mesh("obj_optimized.off", mesh, vertex_color = True)
    # Update the normal vectors of all the vertices on O
    #Utile pour nous?

if __name__ == "__main__":
    main() 