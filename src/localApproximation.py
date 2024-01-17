import openmesh as om
import numpy as np

def developabilityDetectFunction(mesh, minDevelopability=0, maxDevelopabilty=0):
    '''
    Add color to the vertices of a mesh according to their developability
    return the max and min developability values of the mesh 
    '''
    developability = 0
    nbVertices = 0
    listSum = []
    for vertex in mesh.vertices():
        nbVertices += 1
        if not mesh.is_boundary(vertex):
            sum = 2 * np.pi
            neighbours = [vh for vh in mesh.vv(vertex)]
            for i in range(len(neighbours)):
                vh1 = neighbours[i]
                vh2 = neighbours[(i + 1) % len(neighbours)]
                vector1 = vectorFixed(mesh, vertex, vh1)
                vector2 = vectorFixed(mesh, vertex, vh2)
                if ((np.linalg.norm(vector1)) == 0 or (np.linalg.norm(vector2) == 0)):
                    raise ZeroDivisionError("Cannot divide by zero")
                vector1 /= np.linalg.norm(vector1)
                vector2 /= np.linalg.norm(vector2)
                angle = np.arccos(np.dot(vector1, vector2))
                sum -= angle
            listSum.append(abs(sum))
            developability += abs(sum)
        else:
            listSum.append(False)
    if maxDevelopabilty == 0:
        maxSum = max(listSum)
    else:
        maxSum = max(maxDevelopabilty, max(listSum))

    minSum = 0

    for vertex in mesh.vertices():
        if listSum[vertex.idx()]:
            normalizedSum = (listSum[vertex.idx()] - minSum) / (maxSum - minSum)
            r, g, b = (normalizedSum, normalizedSum, normalizedSum)
            color = [r, g, b, 1.]
        else:
            color = [0., 0., 0., 1.]
        mesh.set_color(vertex, color)
    developability = developability/nbVertices
    print("DEVELOPABILITE MOYENNE: ",developability)
    return minSum, maxSum


def vector(mesh, vertexStart, vertexEnd):
    '''
    Return the vector vertexEnd - vertexStart taking into account the movement scale of each of the vertices
    '''
    position_vertexStart = getVertexNewPosition(mesh, vertexStart)
    position_vertexEnd = getVertexNewPosition(mesh, vertexEnd)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])


def vectorFixed(mesh, vertexStart, vertexEnd):
    '''
    Return the vector vertexEnd - vertexStart without taking into account the movement scale of each of the vertices
    '''
    position_vertexStart = mesh.point(vertexStart)
    position_vertexEnd = mesh.point(vertexEnd)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])


def new_vector(mesh, vertexStart, vertexEnd, delta, movement_on_start):
    '''
    Return the vector vertexEnd - vertexStart taking into account the movement scale of each of the vertices 
    and a new delta movement scale on one of the vertices
    if movement_on_start is True the delta movement is applied on the vertexStart
    else the delta movement is applied on the vertexEnd 
    '''
    if (movement_on_start == True):
        normal = mesh.normal(vertexStart)
        old_position = mesh.point(vertexStart)
        position_vertexStart = old_position + delta * normal
        position_vertexEnd =  getVertexNewPosition(mesh, vertexEnd)
    else:
        normal = mesh.normal(vertexEnd)
        old_position = mesh.point(vertexEnd)
        position_vertexEnd = old_position + delta * normal
        position_vertexStart =  getVertexNewPosition(mesh, vertexStart)
    return (position_vertexEnd[0] - position_vertexStart[0],
            position_vertexEnd[1] - position_vertexStart[1],
            position_vertexEnd[2] - position_vertexStart[2])


def getDevelopability(mesh, vertex):
    ''' 
    Return developability of a vertex from the mesh taking into account the saved movement scale of every vertices
    this is the function named g in the paper
    '''
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            vector1 = vector(mesh, vertex, vh1)
            vector2 = vector(mesh, vertex, vh2)
            if ((np.linalg.norm(vector1)) == 0 or (np.linalg.norm(vector2) == 0)):
                raise ZeroDivisionError("Cannot divide by zero")
            vector1 /= np.linalg.norm(vector1)
            vector2 /= np.linalg.norm(vector2)
            angle = np.arccos(np.dot(vector1, vector2))
            developability -= angle
        developability = developability
        return developability
    else:
        return 0


def getNewDevelopability(mesh, vertex, delta): 
    ''' 
    Return developability of a vertex from the mesh (function named g in the paper) with a new movement scale applied on the vertex
    This is g(q + δn q ) in the paper
    '''
    if not mesh.is_boundary(vertex): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(vertex)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            vector1 = new_vector(mesh, vertex, vh1, delta, True)
            vector2 = new_vector(mesh, vertex, vh2, delta, True)
            denom = np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))
            if (denom == 0):
                raise ZeroDivisionError("Cannot divide by zero")
            angle = np.arccos(np.dot(vector1, vector2) / denom)
            developability -= angle
        developability = developability
        return developability
    else:
        return 0


def getNewDevelopabilityNeighboors(mesh, neighbour , moved_vertex, delta):
    ''' 
    Return developability of the adjacent vertex called neighbour of a moved_vertex from the mesh
    delta being a new movement scale applied on the moved_vertex
    '''
    if not mesh.is_boundary(neighbour): #else developability = 0
        developability = 2 * np.pi
        neighbours = [vh for vh in mesh.vv(neighbour)]
        for i in range(len(neighbours)):
            vh1 = neighbours[i]
            vh2 = neighbours[(i + 1) % len(neighbours)]
            if (vh1 == moved_vertex): #Le voisin est notre vertex qui a bougé
                #On prend en compte le mouvement scale de celui ci dans calcul de developpabilité de notre vertex
                vector1 = new_vector(mesh, neighbour, vh1, delta, False)
            else:
                vector1 = vector(mesh, neighbour, vh1)
            if (vh2 == moved_vertex): #Le voisin est notre vertex qui a bougé
                #On prend en compte le mouvement scale de celui ci dans calcul de developpabilité de notre vertex
                vector2 = new_vector(mesh, neighbour, vh2, delta, False)
            else:
                vector2 = vector(mesh, neighbour, vh2)
            denom = np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))
            if (denom == 0):
                raise ZeroDivisionError("Cannot divide by zero")
            angle = np.arccos(np.dot(vector1, vector2) / denom)
            developability -= angle
        developability = developability
        return developability
    else:
        return 0


def initMesh(mesh):
    ''' 
    Set two new properties to the vertices of our mesh: the movement scale along each vertex's unit normal n, initialized to 0 
    and the vertices' developability
    '''
    for vertex in mesh.vertices():
        mesh.set_vertex_property('movement_scale', vertex, 0)
        developability = getDevelopability(mesh, vertex)
        mesh.set_vertex_property('developability', vertex, developability)
    return mesh


def initDic(mesh):
    '''
    Take a mesh and returns a dictionary keyed on the vertices' id, with the [g(· · · )]2 measure as values 
    '''
    vertices = {vertex.idx() : [(mesh.vertex_property('developability')[vertex.idx()])**2, vertex] for vertex in mesh.vertices()}
    # print(vertices) #debug
    return vertices


def getVertex(mesh, dictionary):
    """
    return the vertex with the maximum [g(· · · )]2
    """
    max_value_vertex_id = None
    max_value = -1

    for vertex_id, value in dictionary.items():
        if value[0] > max_value:
            max_value = value[0]
            max_value_vertex_id = vertex_id
    return dictionary[max_value_vertex_id][1]


def updateDevelopability(mesh, moved_vertex, dictionary):
    """
    Update the developability of the moved_vertex and its neighbours in the dictionary 
    """
    vertices_to_update = [vh for vh in mesh.vv(moved_vertex)]
    vertices_to_update.append(moved_vertex)
    for vertex in vertices_to_update:
        developability = getDevelopability(mesh, vertex)
        mesh.set_vertex_property('developability', vertex, developability)
        if (dictionary.get(vertex.idx())) != None: 
            # In case one of the vertices had been removed from the dictionnary
            dictionary[vertex.idx()][0] = developability**2
    return dictionary


def localOptimizationConstraint(mesh, moved_vertex, delta):
    """
    return the constraint T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2 of local optimization
    """
    new_vertex_developability = getNewDevelopability(mesh, moved_vertex, delta) # g(q + δn q )
    neighbours = [vh for vh in mesh.vv(moved_vertex)] #les qj
    sum = 0
    for v in neighbours:
        sum += getNewDevelopabilityNeighboors(mesh, v , moved_vertex, delta)**2 #Calculer la nouvelle developabilité des voisins en prenant en compte movement scale du vertex
    constraint = new_vertex_developability**2 + sum
    return constraint

def updateVertex(mesh, vertex, maxMovementScale):
    """
    Compute and update the movement scale of the given vertex in the mesh
    δ = δ0 - T(δ0)/dT(δ0)
    return 0 if dT(δ0) = 0 and the vertex need to be fixed
    else return 1
    """
    # Compute the movement scale
    prev_movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()] # delta_0
    constraint = localOptimizationConstraint(mesh, vertex, prev_movement_scale) # T(delta_0)
    derivative_constraint = DerivativelocalOptimizationConstraint(mesh, vertex, prev_movement_scale)# dT(delta_0)
    if (derivative_constraint != 0):
        new_movement_scale = prev_movement_scale - constraint / derivative_constraint
        if new_movement_scale > maxMovementScale:
            new_movement_scale = maxMovementScale
        elif new_movement_scale < -maxMovementScale:
            new_movement_scale = -maxMovementScale
        mesh.set_vertex_property('movement_scale', vertex, new_movement_scale) #update the movement scale

    else:
        mesh.set_vertex_property('movement_scale', vertex, 0) #update the movement scale to 0
        print( " /!\ Le vertex numéro ", vertex.idx(), " a été fixé")
        # We remove the vertices from the dict so we don't iterate on it because we fix it
        return 0
    return 1


def DerivativelocalOptimizationConstraint(mesh, moved_vertex, delta):
    """
    Return the constraint derivative constrainte (where the constraint is T(δ) = (g(q + δn q ))2 + sum_j(g(q j ))2)
    Using central difference to approximate dT: f(x + h) - f(x - h)) / (2 * h)
    """
    h = 1e-5
    term1 = localOptimizationConstraint(mesh, moved_vertex, delta + h)
    term2 = localOptimizationConstraint(mesh, moved_vertex, delta - h)
    derivative_constraint = (term1 - term2)/(2*h)
    return derivative_constraint


def updateVerticesPositions(mesh):
    """
    Update the positions of all the vertices by their movement scales
    """
    for vertex in mesh.vertices():
        new_position = getVertexNewPosition(mesh, vertex)
        mesh.set_point(vertex, new_position)


def getVertexNewPosition(mesh, vertex):
    """
    Return the new position of a vertex according to its movement scale
    """
    current_position = mesh.point(vertex)
    if (np.isnan(current_position[0])):
        print("La position du vertex numéro ", vertex.idx(), "vaut nan")
    movement_scale = mesh.vertex_property('movement_scale')[vertex.idx()]
    if movement_scale == None: #A cause de l'initialisation pour pas parcourir deux fois le mesh dans initMesh
        movement_scale = 0
    normal = mesh.normal(vertex)
    new_position = current_position + movement_scale * normal
    return new_position

def main():
    print("---------- Début main\n")
    # Parameters
    maxIter = 500 #40
    nbIter = 0
    epsilon = 0.001
    maxMovementScale = 0.5
    print("Nombre d'itérations max: ", maxIter)
    print("Epsilon = ", epsilon)

    ###    Read .off file
    filename = "../Objects/plane2.off"
    print("lecture du fichier: ", filename)
    mesh = om.read_trimesh(filename)

    #Add colors to see developability
    initialMinDevelopability, initialMaxDevelopability = developabilityDetectFunction(mesh)
    om.write_mesh("obj_initial.off", mesh, vertex_color = True)

    ###    Compute the vertex developability detect function g(q) at each vertex on the given mesh patch and set the movement scale to each vertex to 0
    initMesh(mesh)

    ###    Compute the unit normal n of each vertex
    mesh.update_vertex_normals()
    # print(mesh.vertex_normals()) #pour debug

    ###    Place all vertices in a dictionary vertices_dic keyed on the vertices' id and with the [g(· · · )]2 measure as values and get the vertex with max value
    vertices_dic = initDic(mesh)
    vertex = getVertex(mesh, vertices_dic)
    max_developability = mesh.vertex_property('developability')[vertex.idx()]
    initial_max_developability = max_developability
    liste_max_dev = []

    #while (the developability of the vertex with max developability is greater than ε) and (nbIter < maxIter);
    print("---------- Début algo\n")
    while ((abs(max_developability) > epsilon) and (nbIter < maxIter) ):
        print("Iteration numéro: ", nbIter)
        print("L'id du vertex à optimiser:", vertex.idx())
        print("Sa developabilité : ", max_developability)
        liste_max_dev.append(max_developability)

        ###    Update worst_vertex's movement scale along its unit normal n according to Eq. 17;
        is_fixed = updateVertex(mesh, vertex, maxMovementScale)
        if (is_fixed == 0):
            #this vertex need to be fixed so we remove it from the dic:
            vertices_dic.pop(vertex.idx())
        else:  
            ###    Update the dictionary to reflect the movement on vertex and its adjacent vertices
            max_developability = mesh.vertex_property('developability')[vertex.idx()]
            vertices_dic = updateDevelopability(mesh, vertex, vertices_dic) #ICI
            max_developability = mesh.vertex_property('developability')[vertex.idx()]

        ### Get the next vertex with max developability
        vertex = getVertex(mesh, vertices_dic)
        max_developability = mesh.vertex_property('developability')[vertex.idx()]
        nbIter += 1
    print("---------- Fin algo\n")
    print("Initial maximum developability = ", initial_max_developability)
    print("Final maximum developability = ", max_developability)
    print("Final iterations count = ", nbIter)

    ### Update the positions of all the vertices by their movement scales
    updateVerticesPositions(mesh)

    #Add colors to see the new mesh and its developability
    developabilityDetectFunction(mesh, initialMinDevelopability, initialMaxDevelopability)
    om.write_mesh("obj_optimized.off", mesh, vertex_color = True)

    # Update the normal vectors of all the vertices
    mesh.update_vertex_normals()

    #Plot max developability evolution
    # plt.plot(liste_max_dev, marker='o', linestyle='-', color='b')
    # plt.xlabel('Itération')
    # plt.ylabel('Développabilité max')
    # plt.title('Tracé de la valeur de développabilité maximale à chaque itération')
    # plt.show()

if __name__ == "__main__":
    main() 


    