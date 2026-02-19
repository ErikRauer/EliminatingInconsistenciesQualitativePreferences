"""MinPrefFBAS_Exact.py

Functions to find the minimum number of preferences (or edges) to be removed from a SL-IPG in order to make the preferences consistent
These implement the algorithms described in the "Iterative Exact Method For MP and MD" Section
"""

from enum import Enum
from ortools.linear_solver import pywraplp
import networkx as nx
import time
import heapq as hq
import itertools

"""# Graph Functions"""

#function to create a preference graph from a text file consisting of edge entries
#seperated by new lines. Each line of the file contains one edge entry that consists
#of 3 numbers: from vertex, to vertex, preferences (as a binary string)
#example entry of edge from 1->2 with preferences 0,2,3 is stored as 1 2 1011
# INPUT:
#     filename: name of the PFG input file
#     outputFile: name of the file to print the data to, defaults to "output.txt"
# RETURNS:
#     G: the SL-IPG as a networkX DiGraph
def prefFileToGraph(filename, outputFile = "output.txt"):
  # Create an empty directed graph
  G = nx.DiGraph()
  numEdges = 0
  numPrefs = 0

  # Open the file
  with open(filename, 'r') as file:
      for line in file:
          # Split the line into components: from vertex, to vertex, and preferences
          from_vertex, to_vertex, preferences = line.strip().split()

          # Ensure the vertices are integers and preferences is stored as a string
          from_vertex = int(from_vertex)
          to_vertex = int(to_vertex)

          # Add the edge with the 'preferences' attribute
          G.add_edge(from_vertex, to_vertex, preferences=preferences)
          numEdges += 1

  with open(outputFile, 'w') as file:
    file.write("----------- GRAPH -----------")
    file.write("\nNumber of Preferences: " + str(getNumPrefs(G)))
    file.write("\nNumber of Outcomes: " + str(G.number_of_nodes()))
    file.write("\nNumber of Edges: " + str(numEdges) + "\n")
  return G

#function to create a cycle matrix for graph G as a dictionary with keys for every
#edge in the graph and values corresponding to a list of the cycle entries of that edge.
#ex. CM[i][j] means that edge i partakes in cycle j
# INPUT:
#     G: the SL-IPG
# RETURNS:
#     CM: the cycle matrix as a dictionary from vertex pair edges to (currently empty) lists
def initCycleMatrix(G):  #FIX TO HAVE AN ENTRY PER EDGE NOT PREFERECE
  CM = {}  #create a new dictionary for cycle matrix
  for edge in G.edges():
    CM[edge] = []  # Initialize each key with an empty list or any desired default value
  return CM

#Function to remove a SINGLE preference (pref) from a graph (G) and delete edges that,
#once the preference is removed, have no preferences left
#INPUT: G: graph
#       pref: a single number representing the index of the preference in the preference list
# RETURNS:
#     newG: SL_IPG copy of G, but with the edges induced by pref removed
def removePref(G, pref):
  edges_to_remove = []  # Collect edges to remove
  newG = G.copy()   #create a copy of G so we don't affect original graph
  for u, v in newG.edges():
      if (newG.edges[u, v]['preferences'])[pref] == '1':  # Check for constraint 'pref'
          # Modify the 'preferences' string
          prefsAsList = list(newG.edges[u, v]['preferences'])
          prefsAsList[pref] = '0'
          newG.edges[u, v]['preferences'] = ''.join(prefsAsList)

          # Check if all preferences are 0
          if all(char == '0' for char in newG.edges[u, v]['preferences']):
              edges_to_remove.append((u, v))  # Mark edge for removal

  # Remove edges after iteration
  newG.remove_edges_from(edges_to_remove)
  return newG

#Function to remove MULTIPLE preferences (prefs) from a graph (G) and delete edges that,
#once the preference is removed, have no preferences left
#INPUT: G: graph
#       prefs: a list of integers representing the indexes of preferences in the prference list
# RETURNS:
#     newG: SL_IPG copy of G, but with the edges induced by prefs removed
def removePreferencesFrom(G,prefs):
  edges_to_remove = [] #collect edges to remove
  newG = G.copy() #create a copy of G so we don't affect original
  for u, v in newG.edges():
    for pref in prefs:
      if (newG.edges[u, v]['preferences'])[pref] == '1':  # Check for constraint 'pref'
        # Modify the 'preferences' string
        prefsAsList = list(newG.edges[u, v]['preferences'])
        prefsAsList[pref] = '0'
        newG.edges[u, v]['preferences'] = ''.join(prefsAsList)

        # Check if all preferences are 0
        if all(char == '0' for char in newG.edges[u, v]['preferences']):
            edges_to_remove.append((u, v))  # Mark edge for removal

  # print("Removing edges:")
  for edge in edges_to_remove:
    prefs = G.edges[edge[0],edge[1]]['preferences']
    # print(edge, prefs)
  newG.remove_edges_from(edges_to_remove)
  return newG


#Function to get a list of the edges deleted when MULTIPLE preferences (prefs) from a graph (G) are removed
#INPUT: G: graph
#       prefs: a list of integers representing the indexes of preferences in the prference list
# RETURNS:
#     edges_to_remove: list of the edges induced by prefs that would be removed if prefs was removed
def edgesDeletedWhenPrefsRemoved(G,prefs):
  edges_to_remove = [] #collect edges to remove
  newG = G.copy() #create a copy of G so we don't affect original
  for u, v in newG.edges():
    for pref in prefs:
      if (newG.edges[u, v]['preferences'])[pref] == '1':  # Check for constraint 'pref'
        # Modify the 'preferences' string
        prefsAsList = list(newG.edges[u, v]['preferences'])
        prefsAsList[pref] = '0'
        newG.edges[u, v]['preferences'] = ''.join(prefsAsList)

        # Check if all preferences are 0
        if all(char == '0' for char in newG.edges[u, v]['preferences']):
            edges_to_remove.append((u, v))  # Mark edge for removal

  return edges_to_remove

#function to find a feedback arcset in graph G and return a set of all the preferences
#on those edges (as indicies of the 'preferences' string) as well as the backedges
#Look at section 2 of the research paper for heuristics on how to do so
#Current way: Do a DFS and add all back edges
# INPUT:
#     G: the SL-IPG
# RETURNS:
#     back_edges: list of back edges making up a feedback arc set of G
#     back_edge_prefs: list of preferences inducing the edges in back_edges
def getFeedbackPrefset(G):
    visited = set()
    in_stack = set()          # To track nodes in the recursion stack
    back_edges = []
    back_edge_prefs = set()   #set to store the (indicies of) preferences on backedges
    numPrefs = getNumPrefs(G)

    def dfs(node):
        visited.add(node)
        in_stack.add(node)

        for neighbor in G.successors(node):  # Use successors() for directed graphs
            if neighbor not in visited:
                dfs(neighbor)
            elif neighbor in in_stack:  # Back edge found
                back_edges.append((node, neighbor))
                curBackedgePrefStr = ((G.get_edge_data(node, neighbor))['preferences']) #get the binary string of the preferences on the current backedge
                for curPref in range(0,numPrefs):
                  if (curBackedgePrefStr[curPref] == '1'):
                    back_edge_prefs.add(curPref) #add the prefs found on the backedge to the Backedge prefset

        in_stack.remove(node)  # Done processing this node

    # Run DFS from every unvisited node to cover disconnected components
    for node in G.nodes():
        if node not in visited:
            dfs(node)

    return back_edges, back_edge_prefs

#adds one more row to graph G's cycle matrix CM by using the feedback arcset F
#G: networkX graph
#F: array of tuples representing edges in the feedback arcset
#CM: dictionary whose keys are edges and values are arrays representing the cycle matrix
#This is "algorithm 2" on the Ali et al. paper https://doi.org/10.1145/3446429
def extendCycleMatrix(G, F, CM):
  for e in F:
    u, v = e  # u is head of e, v is tail of e

    # Find the shortest path from u to v in G
    try:
      # Use BFS to find the shortest path (will throw exception if none exists)
      p = nx.shortest_path(G, source=v, target=u, weight=None)

      # Create a cycle by adding edge e to the path p
      cycle_verticies = p + [v]  # This completes the cycle by connecting back to v

      #Create an array of the edges in the cycle
      cycle_edges = []
      for i in range(1,len(cycle_verticies)):
        cycle_edges.append((cycle_verticies[i-1], cycle_verticies[i]))

      #Add a new row to the cycle matrix CM for cycle "cycle"
      for edgeEntry in CM:
        if (edgeEntry in cycle_edges):
          CM[edgeEntry].append(1) #if the key edge is in the cycle, make it 1
        else:
          CM[edgeEntry].append(0)


    except nx.NetworkXNoPath:
      # If no path exists from v to u, skip this edge e
      print("error: no path found from " + str(v) + " to " + str(u))
      continue

  return CM



#function to print graph data of graph G
# INPUT:
#     G: the SL-IPG
def printGraph(G):
  # Display nodes and edges
  print("Nodes:", G.nodes())
  print("Edges with weights:")
  for u, v, weight in G.edges(data=True):
    print(f"({u}, {v}, {weight['preferences']})")

#function to get the number of preferences in the pref graph (length of pref string that labels each edge)
# INPUT:
#     G: the SL-IPG
def getNumPrefs(G):
  first_edge = list(G.edges(data=True))[0]  # (u, v, attr_dict)
  attr_dict = first_edge[2]  # Get the attr_dict (third element of the tuple)

  # Access the 'preferences' attribute
  numPrefs = len(attr_dict['preferences'])  # Get the length of the 'preferences' string
  # print(attr_dict['preferences'], numPrefs)
  return numPrefs

# Converts a dictionary mapping edge indexes to 0 or 1 to a the set of edges whose value is 1 in the dictionary
# INPUT:
#     cycleMatrix: the cycle matrix thus far
#     indexDict: dictionary mapping edge indexes to 0 or 1
# RETURNS:
#     edgeSet: set of edges as vertex pairs whose values where 1 in the indexDict
def getEdgesFromIndex(cycleMatrix, indexDict):
  curEdgeIndex = 0
  edgeSet = []
  for key in cycleMatrix:
    if indexDict[curEdgeIndex] == 1:
      edgeSet.append(key)
    curEdgeIndex += 1
  return edgeSet


"""# Min Feedback Pref Set Methods"""

# Enum used to tell function which algorithm to run
class Algorithm(Enum):
  MFAS = 1      # Normal minimum feedback arc set
  PREFS = 2     # Minimize number of preferences removed
  EDGES = 3     # Minimize number of edges removed by removing preferences

#Function to find a minimum feedback pref set based on integer programing and lazy constraint generation
#this is "algorithm 1" from the research paper
#input: G, a directed graph with m arcs and nonnegative arc weights w_j (j=1,2,...,m)
#       objMinPrefs, a boolean that says whether to minimize the number of preferences or the number of edges of the IPG being removed
#       alg, Algorithm Enum which determines which algorithm to use
#output: A minimum weight feedback pref set
#Variables:
  #allPrefLabels = a list of names of all the preferences in this graph
  #feedback_prefset = feedback prefset at any given time [in research paper this is F^i]
  #zLower = lower bound on the objective function (begins as 0)
  #zUpper = upper bound on the objective function (begins as sum(w_j*y_j))
  #curBestSol = the dictionary that contains the current FBPS with the fewest preferences (y^(i) in research paper)
  #curSolPrefs = the list of preferences in curBestSol
  #curZ = the minimized sum(w_j*y_j) up to this point gotten by summing the entries of curBestSol (z^(i) in research paper)
def findMinFBPS(G, alg):
  startTime = time.time()  #get the time I started the algo to calculate total runtime
  backEdges, feedback_prefset = getFeedbackPrefset(G) #compute F0 of G using heuristic
  numPrefs = getNumPrefs(G)
  numCyclesHist = []
  numConstraintsHist = []

  #set solution associated with F(0) as the incumbent y.
  curBestSolPrefs = dict()
  for pref in range(0,numPrefs):
    if pref in feedback_prefset:
      curBestSolPrefs[pref] = 1
    else:
      curBestSolPrefs[pref] = 0

  #set zLower to 0 and zUpper to sum(w_j*y_j)
  zLower = 0
  zUpper = 0
  if alg == Algorithm.MFAS:
    zUpper = len(backEdges)
  elif alg == Algorithm.PREFS or alg == Algorithm.GREEDY:
    zUpper = len(feedback_prefset) #since w_j=1 and y=1 if the pref is in the FBPS and 0 otherwise, sum w_iy_i = len(FBPS)
  elif alg == Algorithm.EDGES:
    zUpper = len(edgesDeletedWhenPrefsRemoved(G, feedback_prefset))
  #call algorithm 2 with G, F0, and empty cycle matrix to get first cycle matrix A1
  cycle_matrix = initCycleMatrix(G)   #initialize cycle matrix as empty
  cycle_matrix = extendCycleMatrix(G, backEdges, cycle_matrix) #get the first cycle matrix (A^1 in research paper)

  curEdgeIndex = 0
  curBestSolEdges = dict()
  for key in cycle_matrix:
    if key in backEdges:
      curBestSolEdges[curEdgeIndex] = 1
    else:
      curBestSolEdges[curEdgeIndex] = 0
    curEdgeIndex += 1

  iterNum = 1 #variable to keep track of your iteration in the while loop
  while iterNum >0:    #on research paper this is "for i=1,2,..., do"
    # print(cycle_matrix)
    print("Iteration Number: ", iterNum)
    curNumConstraints = 0
    curSolPrefs = {}
    curSolEdges = {}
    curZ_i = 0
    numCyclesHist.append(len(list(cycle_matrix.values())[0]))
    if alg == Algorithm.GREEDY:
      curSolPrefs, curZ_i, curSolEdges, curNumPrefSets = greedyApprox(cycle_matrix, G)
      numConstraintsHist.append(curNumPrefSets)
    else:
      curSolPrefs, curZ_i, curSolEdges, curNumConstraints = CMsolverILP(cycle_matrix, G, alg) #solve the relaxed problem P(i). results: solution yi, the associated feedback arc set S and objective value z(i)
      numConstraintsHist.append(curNumConstraints)
    curPrefSet = yToPrefs(curSolPrefs)
    #TODO-OPTIONAL: when integer programming solver is invoked on the line just above, y can be used as a starting point
    zLower = max(zLower, curZ_i) #set the lower bound zLower to max(zLower, zi)
    # print(zLower, zUpper)
    if(zLower == zUpper): #if zLower == zUpper
      totalTime = time.time() - startTime   #stop clock and determine runtime
      print("exited on condition zLower == zUpper")
      if not alg == Algorithm.MFAS:
        removedEdges = edgesDeletedWhenPrefsRemoved(G, yToPrefs(curBestSolPrefs))
      else:
        removedEdges = getEdgesFromIndex(cycle_matrix, curBestSolEdges)
      displayStats(iterNum, cycle_matrix, yToPrefs(curBestSolPrefs), removedEdges, totalTime, alg, numCyclesHist, numConstraintsHist)
      return yToPrefs(curBestSolPrefs), removedEdges #stop (yHat is optimal)  TODO: is y_i different from yHat? if so curSolPrefs=y_i, not yHat

    #Let G(i) denote the graph obtained by removing all preferences of S from original graph G
    if alg == Algorithm.MFAS:
      G_i = G.copy()
      G_i.remove_edges_from(getEdgesFromIndex(cycle_matrix, curSolEdges))
    else:
      G_i = removePreferencesFrom(G, curPrefSet) #set G_i to be G-{curSolPrefs} without affecting G

    #if Gi can be topologically sorted
    if nx.is_directed_acyclic_graph(G_i):
      totalTime = time.time()-startTime   #stop clock and determine runtime
      print("exited on condition nx.is_directed_acyclic_graph(G_i)")
      if not alg == Algorithm.MFAS:
        removedEdges = edgesDeletedWhenPrefsRemoved(G, curPrefSet)
      else:
        removedEdges = getEdgesFromIndex(cycle_matrix, curSolEdges)
      displayStats(iterNum, cycle_matrix, curPrefSet, removedEdges, totalTime, alg, numCyclesHist, numConstraintsHist)
      return curPrefSet, removedEdges   #stop (yi is the optimal solution to P as well)

    #Otherwise, there are still cycles so compute FBPS F(i) of G(i) using heuristic (my function)
    backEdges, feedback_prefset = getFeedbackPrefset(G_i)

    #Set those components of yi to 1 that correspond to a pref in Fi
    for pref in feedback_prefset:
      curSolPrefs[pref] = 1
    curEdgeIndex = 0
    for key in cycle_matrix:
      if key in backEdges:
        curSolEdges[curEdgeIndex] = 1
      curEdgeIndex += 1

    #yi is now a feasible solution to P
    if alg == Algorithm.PREFS or alg == Algorithm.GREEDY:
      curZ = sumDictionary(curSolPrefs)   #Let zHat be the new objective value at yi
    elif alg == Algorithm.EDGES:
      curZ = len(edgesDeletedWhenPrefsRemoved(G, yToPrefs(curSolPrefs)))
    else:
      curZ = sumDictionary(curSolEdges)
    if curZ < zUpper:
      zUpper = curZ
      curBestSolPrefs = curSolPrefs      #set y to yi
      curBestSolEdges = curSolEdges
    lenCMBeforeExtension = len(next(iter(cycle_matrix.values())))  #get the length of the cycle matrix before extending this iteration
    cycle_matrix = extendCycleMatrix(G_i, backEdges, cycle_matrix) #Call algorithm 2 with Gi, Fi, and Ai to get the extended cycle matrix A(i+1)
    #A(i+1) is guarenteed to have at least one additional row compared to A(i) so throw error if it doesn't
    if(lenCMBeforeExtension >= len(next(iter(cycle_matrix.values())))):
      print("Logic error: ExtendCycleMatrix() didn't add any rows to the Cycle Matrix")
    iterNum += 1
  print("exited on return at end of function (unexpected)")
  return #This should never be executed

#Improved version of findMinFBPS to allow more control of the problem through parameters
#same as findMinFBPS if all parameters are left at default values
#Function to find a minimum feedback pref set based on integer programing and lazy constraint generation
#this is "algorithm 1" from the Ali et al. paper
# INPUT: 
#     G: a directed graph with m arcs and nonnegative arc weights w_j (j=1,2,...,m)
#     alg: Algorithm Enum which determines which ILP to use
#     outputFile: name of output file to write data to
#     useMinPrefSize: whether to use a previous min feedback preference set solution as a constraint
#     minPrefSize: this is the number of preferences removed by a previous solution, only used if useMinPrefSize is true
#     useMinEdgeSize: whether to use a previous min feedback edge set solution as a constraint
#     minEdgeSize: this is the number of edges removed by a previous solution, only used if useMinEdgeSize is true
#     minimizeObjective: whether to minimize or maximize the objective function
# RETURNS: 
#     curPrefSet: A minimum weight feedback pref set
#     removedEdges: Edges removed when preferences of curPrefSet are removed
#Variables:
  #allPrefLabels = a list of names of all the preferences in this graph
  #feedback_prefset = feedback prefset at any given time [in research paper this is F^i]
  #zLower = lower bound on the objective function (begins as 0)
  #zUpper = upper bound on the objective function (begins as sum(w_j*y_j))
  #curBestSol = the dictionary that contains the current FBPS with the fewest preferences (y^(i) in research paper)
  #curSolPrefs = the list of preferences in curBestSol
  #curZ = the minimized sum(w_j*y_j) up to this point gotten by summing the entries of curBestSol (z^(i) in research paper)
def improved_findMinFBPS(G, alg, outputFile = "output.txt", useMinPrefSize = False, minPrefSize = 0, useMinEdgeSize=False, minEdgeSize=0, minimizeObjective=True):
  startTime = time.time()  #get the time I started the algo to calculate total runtime
  backEdges, feedback_prefset = getFeedbackPrefset(G) #compute F0 of G using heuristic
  numPrefs = getNumPrefs(G)
  numCyclesHist = []
  numConstraintsHist = []

  #set solution associated with F(0) as the incumbent y.
  curBestSolPrefs = dict()
  for pref in range(0,numPrefs):
    if pref in feedback_prefset:
      curBestSolPrefs[pref] = 1
    else:
      curBestSolPrefs[pref] = 0

  #set zLower to 0 and zUpper to sum(w_j*y_j)
  zLower = 0
  zUpper = 0
  if alg == Algorithm.MFAS:
    zUpper = len(backEdges)
  elif alg == Algorithm.PREFS or alg == Algorithm.GREEDY:
    zUpper = len(feedback_prefset) #since w_j=1 and y=1 if the pref is in the FBPS and 0 otherwise, sum w_iy_i = len(FBPS)
  elif alg == Algorithm.EDGES:
    zUpper = len(edgesDeletedWhenPrefsRemoved(G, feedback_prefset))
  #call algorithm 2 with G, F0, and empty cycle matrix to get first cycle matrix A1
  cycle_matrix = initCycleMatrix(G)   #initialize cycle matrix as empty
  cycle_matrix = extendCycleMatrix(G, backEdges, cycle_matrix) #get the first cycle matrix (A^1 in research paper)

  solver, curNumConstraints, y_p, y_e = setupSolver(cycle_matrix, G, alg, useMinPrefs=useMinPrefSize, minPrefSize=minPrefSize, useMinEdges=useMinEdgeSize, minEdgeSize=minEdgeSize, minimizeObj=minimizeObjective)
  prevNumCycs = 0

  curEdgeIndex = 0
  curBestSolEdges = dict()
  for key in cycle_matrix:
    if key in backEdges:
      curBestSolEdges[curEdgeIndex] = 1
    else:
      curBestSolEdges[curEdgeIndex] = 0
    curEdgeIndex += 1

  iterNum = 1 #variable to keep track of your iteration in the while loop
  while iterNum >0:    #on research paper this is "for i=1,2,..., do"
    # print(cycle_matrix)
    print("Iteration Number: ", iterNum)
    # curNumConstraints = 0
    curSolPrefs = {}
    curSolEdges = {}
    curZ_i = 0
    numCyclesHist.append(len(list(cycle_matrix.values())[0]))
    if alg == Algorithm.GREEDY:
      curSolPrefs, curZ_i, curSolEdges, curNumPrefSets = greedyApprox(cycle_matrix, G)
      numConstraintsHist.append(curNumPrefSets)
    else:
      curSolPrefs, curZ_i, curSolEdges, curNumConstraints = CMsolverILP(cycle_matrix, G, alg, solver, prevNumCycs, curNumConstraints, y_p, y_e) #solve the relaxed problem P(i). results: solution yi, the associated feedback arc set S and objective value z(i)
      numConstraintsHist.append(curNumConstraints)
    curPrefSet = yToPrefs(curSolPrefs)
    #TODO-OPTIONAL: when integer programming solver is invoked on the line just above, y can be used as a starting point
    zLower = max(zLower, curZ_i) #set the lower bound zLower to max(zLower, zi)
    # print(zLower, zUpper)
    if(zLower == zUpper and not useMinPrefSize and not useMinEdgeSize): #if zLower == zUpper
      totalTime = time.time() - startTime   #stop clock and determine runtime
      print("exited on condition zLower == zUpper")
      if not alg == Algorithm.MFAS:
        removedEdges = edgesDeletedWhenPrefsRemoved(G, yToPrefs(curBestSolPrefs))
      else:
        removedEdges = getEdgesFromIndex(cycle_matrix, curBestSolEdges)
      displayStats(iterNum, cycle_matrix, yToPrefs(curBestSolPrefs), removedEdges, totalTime, alg, numCyclesHist, numConstraintsHist, outputFile)
      return yToPrefs(curBestSolPrefs), removedEdges #stop (yHat is optimal)  TODO: is y_i different from yHat? if so curSolPrefs=y_i, not yHat

    #Let G(i) denote the graph obtained by removing all preferences of S from original graph G
    if alg == Algorithm.MFAS:
      G_i = G.copy()
      G_i.remove_edges_from(getEdgesFromIndex(cycle_matrix, curSolEdges))
    else:
      G_i = removePreferencesFrom(G, curPrefSet) #set G_i to be G-{curSolPrefs} without affecting G

    #if Gi can be topologically sorted
    if nx.is_directed_acyclic_graph(G_i):
      totalTime = time.time()-startTime   #stop clock and determine runtime
      print("exited on condition nx.is_directed_acyclic_graph(G_i)")
      if not alg == Algorithm.MFAS:
        removedEdges = edgesDeletedWhenPrefsRemoved(G, curPrefSet)
      else:
        removedEdges = getEdgesFromIndex(cycle_matrix, curSolEdges)
      displayStats(iterNum, cycle_matrix, curPrefSet, removedEdges, totalTime, alg, numCyclesHist, numConstraintsHist, outputFile)
      return curPrefSet, removedEdges   #stop (yi is the optimal solution to P as well)

    prevNumCycs = len(list(cycle_matrix.values())[0])

    #Otherwise, there are still cycles so compute FBPS F(i) of G(i) using heuristic (my function)
    backEdges, feedback_prefset = getFeedbackPrefset(G_i)

    #Set those components of yi to 1 that correspond to a pref in Fi
    for pref in feedback_prefset:
      curSolPrefs[pref] = 1
    curEdgeIndex = 0
    for key in cycle_matrix:
      if key in backEdges:
        curSolEdges[curEdgeIndex] = 1
      curEdgeIndex += 1

    #yi is now a feasible solution to P
    if alg == Algorithm.PREFS or alg == Algorithm.GREEDY:
      curZ = sumDictionary(curSolPrefs)   #Let zHat be the new objective value at yi
    elif alg == Algorithm.EDGES:
      curZ = len(edgesDeletedWhenPrefsRemoved(G, yToPrefs(curSolPrefs)))
    else:
      curZ = sumDictionary(curSolEdges)
    if curZ < zUpper:
      zUpper = curZ
      curBestSolPrefs = curSolPrefs      #set y to yi
      curBestSolEdges = curSolEdges
    lenCMBeforeExtension = len(next(iter(cycle_matrix.values())))  #get the length of the cycle matrix before extending this iteration
    cycle_matrix = extendCycleMatrix(G_i, backEdges, cycle_matrix) #Call algorithm 2 with Gi, Fi, and Ai to get the extended cycle matrix A(i+1)
    #A(i+1) is guarenteed to have at least one additional row compared to A(i) so throw error if it doesn't
    if(lenCMBeforeExtension >= len(next(iter(cycle_matrix.values())))):
      print("Logic error: ExtendCycleMatrix() didn't add any rows to the Cycle Matrix")
    iterNum += 1
  print("exited on return at end of function (unexpected)")
  return #This should never be executed


"""# Integer Linear Programming"""

#Sets up the ILP solver
# INPUT: 
#     G: a directed graph with m arcs and nonnegative arc weights w_j (j=1,2,...,m)
#     CM: cycle matrix found thus far
#     alg: Algorithm Enum which determines which ILP to use
#     useMinPrefSize: whether to use a previous min feedback preference set solution as a constraint
#     minPrefSize: this is the number of preferences removed by a previous solution, only used if useMinPrefSize is true
#     useMinEdgeSize: whether to use a previous min feedback edge set solution as a constraint
#     minEdgeSize: this is the number of edges removed by a previous solution, only used if useMinEdgeSize is true
#     minimizeObjective: whether to minimize or maximize the objective function
# RETURNS: 
#     solver: the setup solver ready to be used
#     num_constraints: number of constraints added to the solver
#     y_p: array of ILP variables corresponding with the preferences
#     y_e: array of ILP variables corresponding with the edges
def setupSolver(CM, G, alg, useMinPrefs=False, minPrefSize=0, useMinEdges=False, minEdgeSize=0, minimizeObj=True):
  solver = pywraplp.Solver.CreateSolver('CBC')

  # Keep track of the number of constraints
  num_constraints = 0

  # Define the variables with integer constraints
  num_edges = len(CM)         #counts the number of keys in the dictionary
  num_prefs = getNumPrefs(G)  #counts number of prefs in the graph
  y_p = [solver.IntVar(0, 1, f'y_p{i}') for i in range(num_prefs)] #creates a ILP variable for all preferences in the graph
  y_e = [solver.IntVar(0, 1, f'y_e{i}') for i in range(num_edges)] #creates an ILP variable for all edges in the graph

  #Define the constraints that ensure edges arent removed unless all preferences labeling that edge are removed
  if not alg == Algorithm.MFAS:
    edgeInd = 0 #variable to iterate through edge ILP variables array (since its a list, not a dict), so y_e[edgeInd] corresponds to the 'edgeInd'th edge (key) in CM when they are "ordered"
    for edgeKey in CM:
      u, v = edgeKey #split the edge edgeKey into its 2 verticies
      curEdgePrefVars = [] # list of extra variables used when trying to minimize the number of edges of IPG being removed
      numPrefs = 0 # number of preferences labeling this edge
      curPrefStr = ((G.get_edge_data(u,v))['preferences']) #get the string stored in the 'preferences' attribute
      for prefInd in range(0,num_prefs):
        if (curPrefStr[prefInd] == '1'):
          solver.Add(y_e[edgeInd] <= y_p[prefInd])
          num_constraints += 1

          # If we want to minimize the number of edges removed, need to ensure edges are removed iff all preferences are removed from it
          if alg == Algorithm.EDGES or useMinEdges:
            # Create a new set of variables for each preference labelling the current edge
            curEdgePrefVars.append(solver.IntVar(0, 1, f'r{edgeInd}_{prefInd}'))
            # First set of constraints ensuring edge variable is same as pref variable iff corresponding new variable is 1
            solver.Add(curEdgePrefVars[numPrefs] - 1 <= y_e[edgeInd] - y_p[prefInd])
            solver.Add(y_e[edgeInd] - y_p[prefInd] <= 1 - curEdgePrefVars[numPrefs])
            num_constraints += 2
          numPrefs += 1
      # If we're minimizing the number of edges removed, need to add another constraint for the current edge
      if alg == Algorithm.EDGES or useMinEdges:
        curEdgeConstraintSum = 0
        for edgePrefVar in curEdgePrefVars:
          curEdgeConstraintSum += edgePrefVar
        # Ensures only one of the new variables is equal to 1
        solver.Add(curEdgeConstraintSum == 1)
        num_constraints += 1
      edgeInd += 1

  if useMinPrefs:
    prefConstraint = 0
    for prefVar in y_p:
      prefConstraint += prefVar
    solver.Add(prefConstraint <= minPrefSize)
    num_constraints += 1

  if useMinEdges:
    edgeConstraint = 0
    for edgeVar in y_e:
      edgeConstraint += edgeVar
    solver.Add(edgeConstraint <= minEdgeSize)
    num_constraints += 1

  equation_to_minimize = 0  #variable to store the objective
  if alg == Algorithm.PREFS:
    for i in range(0,num_prefs):
      equation_to_minimize += y_p[i]
  else:
    for i in range(0, num_edges):
      equation_to_minimize += y_e[i]
  if minimizeObj:
    solver.Minimize(equation_to_minimize) #we want to minimize y_1 + y_2 + ... + y_n
  else:
    solver.Maximize(equation_to_minimize)

  return solver, num_constraints, y_p, y_e



#Function to use Integer Linear Programming to solve the relaxed problem P^(i)
# INPUT: 
#     CM: cycle matrix to be solved (considered P^(i))
#     num_prefs: number of prefs in the graph
#     alg: Algorithm Enum used to determine which objective function and what conditions to use
#     solver: ILP solver
#     num_constraints: number of constraints in the ILP solver
#     y_p: list of preferences variables in the ILP solver
#     y_e: list of edge variables in the ILP solver
#output:
#  solPrefs: list of tuples representing the prefs in the FBPS
#  zVal: objective zVal z^(i) found by doing sum(y_iw_i) from the solY (value of ILP equation we are trying to minimize)
#  solEdges: edges removed via the solution
#  num_constraints: updated number of constraints in the solver
def CMsolverILP(CM, G, alg, solver, prevNumCycs, num_constraints, y_p, y_e):
  # Initialize the solver with integer programming mode
  # solver = pywraplp.Solver.CreateSolver('CBC')

  # # Keep track of the number of constraints
  # num_constraints = 0

  # # Define the variables with integer constraints
  num_edges = len(CM)         #counts the number of keys in the dictionary
  num_prefs = getNumPrefs(G)  #counts number of prefs in the graph
  # y_p = [solver.IntVar(0, 1, f'y_p{i}') for i in range(num_prefs)] #creates a ILP variable for all preferences in the graph
  # y_e = [solver.IntVar(0, 1, f'y_e{i}') for i in range(num_edges)] #creates an ILP variable for all edges in the graph

  # Define the constraints that ensure acyclic graph
  num_CM_rows = len(list(CM.values())[0]) #get the length of array at any key in dict (all should be the same length) which corresponds to num rows in Cycle matrix
  for i in range(prevNumCycs, num_CM_rows): #go through CM & add a constraint of sum(a_i*y_i)>= 1 for each row in the cycle matrix
    constraint_sum = 0 #variable to hold (a_11*y_1 + a_12*y_2+...+a_1n*y_n) to make sure its greater than 1
    j = 0 #variable to iterate through edge ILP variables array (since its a list, not a dictionary), so y_e[i] corresponds to the ith edge (key) in CM when they are "ordered"
    for key in CM:
      constraint_sum += CM[key][i]*y_e[j]
      j+=1 #increment your index in y each time
    solver.Add(constraint_sum >= 1)
    num_constraints += 1

  # #Define the constraints that ensure edges arent removed unless all preferences labeling that edge are removed
  # if not alg == Algorithm.MFAS:
  #   edgeInd = 0 #variable to iterate through edge ILP variables array (since its a list, not a dict), so y_e[edgeInd] corresponds to the 'edgeInd'th edge (key) in CM when they are "ordered"
  #   for edgeKey in CM:
  #     u, v = edgeKey #split the edge edgeKey into its 2 verticies
  #     curEdgePrefVars = [] # list of extra variables used when trying to minimize the number of edges of IPG being removed
  #     numPrefs = 0 # number of preferences labeling this edge
  #     curPrefStr = ((G.get_edge_data(u,v))['preferences']) #get the string stored in the 'preferences' attribute
  #     for prefInd in range(0,num_prefs):
  #       if (curPrefStr[prefInd] == '1'):
  #         solver.Add(y_e[edgeInd] <= y_p[prefInd])
  #         num_constraints += 1

  #         # If we want to minimize the number of edges removed, need to ensure edges are removed iff all preferences are removed from it
  #         if alg == Algorithm.EDGES:
  #           # Create a new set of variables for each preference labelling the current edge
  #           curEdgePrefVars.append(solver.IntVar(0, 1, f'r{edgeInd}_{prefInd}'))
  #           # First set of constraints ensuring edge variable is same as pref variable iff corresponding new variable is 1
  #           solver.Add(curEdgePrefVars[numPrefs] - 1 <= y_e[edgeInd] - y_p[prefInd])
  #           solver.Add(y_e[edgeInd] - y_p[prefInd] <= 1 - curEdgePrefVars[numPrefs])
  #           num_constraints += 2
  #         numPrefs += 1
  #     # If we're minimizing the number of edges removed, need to add another constraint for the current edge
  #     if alg == Algorithm.EDGES:
  #       curEdgeConstraintSum = 0
  #       for edgePrefVar in curEdgePrefVars:
  #         curEdgeConstraintSum += edgePrefVar
  #       # Ensures only one of the new variables is equal to 1
  #       solver.Add(curEdgeConstraintSum == 1)
  #       num_constraints += 1
  #     edgeInd += 1

  # # Define the objective function
  # equation_to_minimize = 0  #variable to store the objective
  # if alg == Algorithm.PREFS:
  #   for i in range(0,num_prefs):
  #     equation_to_minimize += y_p[i]
  # else:
  #   for i in range(0, num_edges):
  #     equation_to_minimize += y_e[i]
  # solver.Minimize(equation_to_minimize) #we want to minimize y_1 + y_2 + ... + y_n

  # Solve the problem
  print(solver.NumVariables(),solver.NumConstraints())
  status = solver.Solve()
  print("ILP Done")

  if status == pywraplp.Solver.OPTIMAL:
    #If successful, Create return variables (solution dictionary Y, list of prefs in solution, etc.)
    solPrefs = {} #create the solution Y (as a dictionary) simmilar to CM
    curYindex = 0  #variable to traverse the y_p ILP variable list
    solEdges = {}

    # for var in solver.variables():
    #   print(var, var.solution_value())

    # print("Prefs being removed:")
    for keyPref in range(0, num_prefs): #go through all the keys, ordered as they were when the y_i's were created
      solPrefs[keyPref] = y_p[curYindex].solution_value()  #assign the pref value in the solution dictionary to its corresponding y_i
      # if y_p[curYindex].solution_value() == 1:     #if the pref is in the FBPS, add it to the solPrefs list
        # print(curYindex)                     #key here is a preference
      curYindex += 1                      #as you go through keys in dictionary, also go through y array
    # print("Solution Pref Vars ", solPrefs)
    # print("Edges being removed:")
    curYindex = 0
    for keyEdge in range(0, num_edges): #go through all the keys, ordered as they were when the y_i's were created
      solEdges[keyEdge] = y_e[curYindex].solution_value()  #assign the pref value in the solution dictionary to its corresponding y_i
      # if y_e[curYindex].solution_value() == 1:
      #   print(curYindex)
      curYindex += 1                      #as you go through keys in dictionary, also go through y array
    # print("Solution Edge Vars", solEdges)
    zVal = solver.Objective().Value() #get the value of the minimized y_1 + y_2 + ... equation
    return solPrefs, zVal, solEdges, num_constraints

  elif status == pywraplp.Solver.FEASIBLE:
      print('A feasible solution was found.')
      return None, None, None, None #TODO, handle this better
  else:
      print('No feasible solution found.')
      return None, None, None, None #TODO, handle this better


#Function to use Integer Linear Programming to solve the relaxed problem P^(i)
#input: CM - cycle matrix to be solved (considered P^(i))
#       num_prefs - number of prefs in the graph
#       alg - Algorithm Enum used to determine which objective function and what conditions to use
#output:
#  solY: dictionary representing the FBAS where prefss have entry 1 if they are in the FBPS and 0 otherwise
#  solPrefs: list of tuples representing the prefs in the FBPS
#  zVal: objective zVal z^(i) found by doing sum(y_iw_i) from the solY (value of ILP equation we are trying to minimize)
def CMsolverILP(CM, G, alg):
  # Initialize the solver with integer programming mode
  solver = pywraplp.Solver.CreateSolver('CBC')

  # Keep track of the number of constraints
  num_constraints = 0

  # Define the variables with integer constraints
  num_edges = len(CM)         #counts the number of keys in the dictionary
  num_prefs = getNumPrefs(G)  #counts number of prefs in the graph
  y_p = [solver.IntVar(0, 1, f'y_p{i}') for i in range(num_prefs)] #creates a ILP variable for all preferences in the graph
  y_e = [solver.IntVar(0, 1, f'y_e{i}') for i in range(num_edges)] #creates an ILP variable for all edges in the graph

  # Define the constraints that ensure acyclic graph
  num_CM_rows = len(list(CM.values())[0]) #get the length of array at any key in dict (all should be the same length) which corresponds to num rows in Cycle matrix
  for i in range(0, num_CM_rows): #go through CM & add a constraint of sum(a_i*y_i)>= 1 for each row in the cycle matrix
    constraint_sum = 0 #variable to hold (a_11*y_1 + a_12*y_2+...+a_1n*y_n) to make sure its greater than 1
    j = 0 #variable to iterate through edge ILP variables array (since its a list, not a dictionary), so y_e[i] corresponds to the ith edge (key) in CM when they are "ordered"
    for key in CM:
      constraint_sum += CM[key][i]*y_e[j]
      j+=1 #increment your index in y each time
    solver.Add(constraint_sum >= 1)
    num_constraints += 1

  #Define the constraints that ensure edges arent removed unless all preferences labeling that edge are removed
  if not alg == Algorithm.MFAS:
    edgeInd = 0 #variable to iterate through edge ILP variables array (since its a list, not a dict), so y_e[edgeInd] corresponds to the 'edgeInd'th edge (key) in CM when they are "ordered"
    for edgeKey in CM:
      u, v = edgeKey #split the edge edgeKey into its 2 verticies
      curEdgePrefVars = [] # list of extra variables used when trying to minimize the number of edges of IPG being removed
      numPrefs = 0 # number of preferences labeling this edge
      curPrefStr = ((G.get_edge_data(u,v))['preferences']) #get the string stored in the 'preferences' attribute
      for prefInd in range(0,num_prefs):
        if (curPrefStr[prefInd] == '1'):
          solver.Add(y_e[edgeInd] <= y_p[prefInd])
          num_constraints += 1

          # If we want to minimize the number of edges removed, need to ensure edges are removed iff all preferences are removed from it
          if alg == Algorithm.EDGES:
            # Create a new set of variables for each preference labelling the current edge
            curEdgePrefVars.append(solver.IntVar(0, 1, f'r{edgeInd}_{prefInd}'))
            # First set of constraints ensuring edge variable is same as pref variable iff corresponding new variable is 1
            solver.Add(curEdgePrefVars[numPrefs] - 1 <= y_e[edgeInd] - y_p[prefInd])
            solver.Add(y_e[edgeInd] - y_p[prefInd] <= 1 - curEdgePrefVars[numPrefs])
            num_constraints += 2
          numPrefs += 1
      # If we're minimizing the number of edges removed, need to add another constraint for the current edge
      if alg == Algorithm.EDGES:
        curEdgeConstraintSum = 0
        for edgePrefVar in curEdgePrefVars:
          curEdgeConstraintSum += edgePrefVar
        # Ensures only one of the new variables is equal to 1
        solver.Add(curEdgeConstraintSum == 1)
        num_constraints += 1
      edgeInd += 1

  # Define the objective function
  equation_to_minimize = 0  #variable to store the objective
  if alg == Algorithm.PREFS:
    for i in range(0,num_prefs):
      equation_to_minimize += y_p[i]
  else:
    for i in range(0, num_edges):
      equation_to_minimize += y_e[i]
  solver.Minimize(equation_to_minimize) #we want to minimize y_1 + y_2 + ... + y_n

  # Solve the problem
  # print(solver.NumVariables(),solver.NumConstraints())
  status = solver.Solve()

  if status == pywraplp.Solver.OPTIMAL:
    #If successful, Create return variables (solution dictionary Y, list of prefs in solution, etc.)
    solPrefs = {} #create the solution Y (as a dictionary) simmilar to CM
    curYindex = 0  #variable to traverse the y_p ILP variable list
    solEdges = {}

    # for var in solver.variables():
    #   print(var, var.solution_value())

    # print("Prefs being removed:")
    for keyPref in range(0, num_prefs): #go through all the keys, ordered as they were when the y_i's were created
      solPrefs[keyPref] = y_p[curYindex].solution_value()  #assign the pref value in the solution dictionary to its corresponding y_i
      # if y_p[curYindex].solution_value() == 1:     #if the pref is in the FBPS, add it to the solPrefs list
        # print(curYindex)                     #key here is a preference
      curYindex += 1                      #as you go through keys in dictionary, also go through y array
    # print("Solution Pref Vars ", solPrefs)
    # print("Edges being removed:")
    curYindex = 0
    for keyEdge in range(0, num_edges): #go through all the keys, ordered as they were when the y_i's were created
      solEdges[keyEdge] = y_e[curYindex].solution_value()  #assign the pref value in the solution dictionary to its corresponding y_i
      # if y_e[curYindex].solution_value() == 1:
      #   print(curYindex)
      curYindex += 1                      #as you go through keys in dictionary, also go through y array
    # print("Solution Edge Vars", solEdges)
    zVal = solver.Objective().Value() #get the value of the minimized y_1 + y_2 + ... equation
    return solPrefs, zVal, solEdges, num_constraints

  elif status == pywraplp.Solver.FEASIBLE:
      print('A feasible solution was found.')
      return None, None, None, None #TODO, handle this better
  else:
      print('No feasible solution found.')
      return None, None, None, None #TODO, handle this better



"""# Additional Helper Functions"""

#function to sum all values in a dictionary.
def sumDictionary(dict):
  total = 0
  for val in dict.values():
    total += val
  return total

#function to turn a solution dictionary (y) wherin preference entries are 1 if
#the preference is included and 0 otherwise, into the list of preferences it represents.
def yToPrefs(solutionDict):
  solPrefs = []
  for key in solutionDict:
    if solutionDict[key] == 1:
      solPrefs.append(key)
  return solPrefs


#function to print off an array that is stored as a dictionary where the keys are indicies
def printDictArr(dictArr):
  for key in dictArr.keys():
    print(key, end = "|")
  print()
  for i in range(0, len(list(dictArr.values())[0])):
    for key in dictArr.keys():
      print(dictArr[key][i], end=" | ")
    print()
  return

def displayStats(numIterations, cycleMatrix, minFBPS, removedEdges, totalTime, alg, numCycles, numConstraints, outputFile = "output.txt"):
  print("----------- STATS -----------")
  print("Algorithm used: " + str(alg))
  print("runtime: " + str(totalTime) + "seconds")
  print("Number of Iterations: " + str(numIterations))
  print("Number of Rows in final Cycle Matrix: " + str(len(list(cycleMatrix.values())[0])))
  print("History of number of rows in the cycle matrix: " + str(numCycles))
  if alg == Algorithm.GREEDY:
    print("Number of preference sets in final Iteration: " + str(numConstraints[numIterations - 1]))
    print("History of number of preference sets: " + str(numConstraints))
  else:
    print("Number of Constraints in final Iteration: " + str(numConstraints[numIterations - 1]))
    print("History of Constraint Nums: " + str(numConstraints))
  # print("Nonzero rows of final Cycle Matrix")
  if not alg == Algorithm.MFAS:
    print("Prefs Removed (size of minFBPS): " + str(len(minFBPS)))
    print("minimum Feedback Prefsetset: " + str(minFBPS))
  print("Number of Edges Removed: " + str(len(removedEdges)))
  print("Feedback Edge Set: " + str(removedEdges))
  print("-----------------------------")
  with open(outputFile,'a') as file:
    file.write("----------- STATS -----------\n")
    file.write("Algorithm used: " + str(alg))
    file.write("\nRuntime: " + str(totalTime) + "seconds\n")
    file.write("Number of Iterations: " + str(numIterations))
    file.write("\nNumber of Rows in Final Cycle Matrix: " + str(len(list(cycleMatrix.values())[0])))
    file.write("\nHistory of Number of Rows in the Cycle Matrix: " + str(numCycles))
    if alg == Algorithm.GREEDY:
      file.write("\nNumber of Preference Sets in Final Iteration: " + str(numConstraints[numIterations - 1]))
      file.write("\nHistory of Number of Preference Sets: " + str(numConstraints))
    else:
      file.write("\nNumber of Constraints in Final Iteration: " + str(numConstraints[numIterations - 1]))
      file.write("\nHistory of Constraint Nums: " + str(numConstraints))
    # print("Nonzero rows of final Cycle Matrix")
    if not alg == Algorithm.MFAS:
      file.write("\nNumber of Prefs Removed (size of minFBPS): " + str(len(minFBPS)))
      file.write("\nMinimum Feedback Pref Set: " + str(minFBPS))
    file.write("\nNumber of Edges Removed: " + str(len(removedEdges)))
    file.write("\nFeedback Edge Set: " + str(removedEdges) + "\n")

"""# Main"""

# G = prefFileToGraph("prefGraph_n5_e8_ID1.txt")
# G = prefFileToGraph("prefGraph_n6_e9_ID2.txt")
# G = prefFileToGraph("prefGraph_n4_e8_ID3.txt")
# G = prefFileToGraph("prefGraph_n9_e8_ID5.txt")
# G = prefFileToGraph("prefGraph_00056-00004001.txt")
# G = prefFileToGraph("n6d2g3-PFG3.txt")
# G = prefFileToGraph("cptheory-n8d2g3p50-PFG1.txt")
# printGraph(G)
# drawPrefGraph(G)

# prefix = "n5d2g3"
for i in range(3, 10):
  for j in range(1,4):
    prefix = "CP-Nets/n" + str(i) + "d2g3/n" + str(i) + "d2g3"
    # prefix = "CP-Theories/n" + str(i) + "d2g3/cptheory-n" + str(i) + "d2g3p50"
    inputFile = prefix + "-PFG" + str(j) + ".txt"
    output_file = prefix + "-" + str(j) + "-max-output.txt"
    print(output_file)

    G = prefFileToGraph(inputFile, output_file)

    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.MFAS, outputFile=output_file)
    newG = G.copy()
    newG.remove_edges_from(removedEdges)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.EDGES, outputFile=output_file)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")
      

    print(len(edgesDeletedWhenPrefsRemoved(G,minFPS)))

    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.PREFS, outputFile=output_file, useMinEdgeSize=True, minEdgeSize=len(edgesDeletedWhenPrefsRemoved(G,minFPS)), minimizeObjective=False)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.PREFS, outputFile=output_file)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.EDGES, outputFile=output_file, useMinPrefSize=True, minPrefSize=len(minFPS), minimizeObjective=False)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")



for i in range(3, 8):
  for j in range(1,4):
    # prefix = "CP-Nets/n" + str(i) + "d2g3/n" + str(i) + "d2g3"
    prefix = "CP-Theories/n" + str(i) + "d2g3/cptheory-n" + str(i) + "d2g3p50"
    inputFile = prefix + "-PFG" + str(j) + ".txt"
    output_file = prefix + "-" + str(j) + "-max-output.txt"
    print(output_file)

    G = prefFileToGraph(inputFile, output_file)

    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.MFAS, outputFile=output_file)
    newG = G.copy()
    newG.remove_edges_from(removedEdges)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.EDGES, outputFile=output_file)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")
      

    print(len(edgesDeletedWhenPrefsRemoved(G,minFPS)))

    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.PREFS, outputFile=output_file, useMinEdgeSize=True, minEdgeSize=len(edgesDeletedWhenPrefsRemoved(G,minFPS)), minimizeObjective=False)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.PREFS, outputFile=output_file)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


    minFPS, removedEdges = improved_findMinFBPS(G.copy(), Algorithm.EDGES, outputFile=output_file, useMinPrefSize=True, minPrefSize=len(minFPS), minimizeObjective=False)
    newG = removePreferencesFrom(G,minFPS)
    if(nx.is_directed_acyclic_graph(newG)):
        print("minFBAS is correct")
    else:
      print("PROBLEM!! minFBAS is NOT correct. cycles still exist")

# minFPS, removedEdges = findMinFBPS(G.copy(), Algorithm.GREEDY)
# newG = removePreferencesFrom(G,minFPS)
# if(nx.is_directed_acyclic_graph(newG)):
#     print("minFBAS is correct")
# else:
#   print("PROBLEM!! minFBAS is NOT correct. cycles still exist")

# drawPrefGraph(newG)
# allPrefs = list(range(0,11))

# printGraph(newG)

# for i in allPrefs:
#   leftoverPrefs = allPrefs.copy()
#   leftoverPrefs.remove(i)

#   for j in leftoverPrefs:
#     testPrefs = leftoverPrefs.copy()
#     testPrefs.remove(j)

#     testG = removePreferencesFrom(G, testPrefs)
#     print("Trying with only " + str(i) + " and " + str(j))
#     if(nx.is_directed_acyclic_graph(testG)):
#       print("-----------------Is acyclic-----------------")
#     else:
#       print("Is not acyclic")


# for removedPref in minFPS:
#   # print(removedPref)
#   testFPS = minFPS.copy()
#   # print(str(testFPS))
#   testFPS.remove(removedPref)
#   # print(str(testFPS))
#   testG = removePreferencesFrom(G,testFPS)
#   print("Trying with " + str(removedPref) + " added back:")
#   if(nx.is_directed_acyclic_graph(testG)):
#     print("Is acyclic")
#   else:
#     print("Is not acyclic")


# G = prefFileToGraph("n9d2g3-PFG1.txt")

# minFPS, removedEdges = findMinFBPS(G.copy(), Algorithm.EDGES)
# newG = removePreferencesFrom(G,minFPS)
# if(nx.is_directed_acyclic_graph(newG)):
#     print("minFBAS is correct")
# else:
#   print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


# G = prefFileToGraph("n9d2g3-PFG2.txt")

# minFPS, removedEdges = findMinFBPS(G.copy(), Algorithm.EDGES)
# newG = removePreferencesFrom(G,minFPS)
# if(nx.is_directed_acyclic_graph(newG)):
#     print("minFBAS is correct")
# else:
#   print("PROBLEM!! minFBAS is NOT correct. cycles still exist")


# G = prefFileToGraph("n9d2g3-PFG3.txt")

# minFPS, removedEdges = findMinFBPS(G.copy(), Algorithm.EDGES)
# newG = removePreferencesFrom(G,minFPS)
# if(nx.is_directed_acyclic_graph(newG)):
#     print("minFBAS is correct")
# else:
#   print("PROBLEM!! minFBAS is NOT correct. cycles still exist") 9903520314283042199192993792