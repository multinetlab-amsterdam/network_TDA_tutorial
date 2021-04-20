#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Topological Data Analysis 

This script is a compilation of some topological data analysis tools.

"""

__author__ = 'Fernando Nobrega & Eduarda Centeno & Giulia Moreni'
__contact__ = 'f.nobregasantos@amsterdamumc.nl or e.centeno@amsterdamumc.nl'
__date__ = '2020/05/10'   ### Date it was created
__status__ = 'Production'


####################
# Review History   #
####################


####################
# Libraries        #
####################

# Standard imports
import itertools

# Third party imports 
import numpy as np # version 1.18.5
import matplotlib.pyplot as plt # version 3.3.2
import scipy.io # version 1.5.0
from sklearn import preprocessing # version 0.24.1
import networkx as nx # version 2.4
import scipy.special

#############################
# Pre-defined settings      #
#############################
# Notice that some TDA scripts are quite heavy, if you run in a server 
# consider using nice command
# niceValue = os.nice(10)

# You can set up the maximum clique size for your analysis bellow.
kmax = 30 # dimensions - max size for the optimized clique algorithm

# Define Functions ----------------------------------------------------

def normalize(matrix):
    """Matrix normalization
    
    Parameters
    ----------
    matrix: numpy matrix
        
    Returns
    -------
    X_scale_maxabs: numpy matrix
        rescaled matrix
    
    """
    # For details in this normalization, see: 
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MaxAbsScaler.html
    # Scale each feature by its maximum absolute value.
    # This estimator scales and translates each feature individually such
    # that the maximal absolute value of each feature in the training set will be 1.0. 
    # It does not shift/center the data, and thus does not destroy any sparsity.


    # This is the data I want to scale
    X_scale = np.copy(matrix)
    # This is the one I can use for the HCP
    
    max_abs_scaler = preprocessing.MaxAbsScaler()
    X_scale_maxabs = max_abs_scaler.fit_transform(X_scale)
   
    return X_scale_maxabs #X_train_minmax


def max_cliques(N, k):
    """
    
    Parameters
    ----------
    N: number of nodes of your network
    
    k: maximum size of the cliques
       
    Returns
    -------
    mclq: total number of possible cliques with size from 0 to k
    
    OBS:
    ---
    The clique algorithm is time consuming (NP) large and dense matrices, 
    this function is an attempt to deal with it
    
    """
    
    mclq = 0
    for i in range(0, k+1):
        # Notice that we sum up to k+1, since Python does not counts 
        # the last values in the range.
        mclq += scipy.special.binom(N, i)
    
    mclq = int(mclq)
    
    return mclq


def Kmaxcliques(G, kmax=kmax):
    """
    
    Parameters
    ----------
    G: networkx graph
    
    kmax: int
        number of dimensions
    
    Returns
    -------
    C: list with all maximal cliques of the graph G with size up to kmax
    
    """
    # Depending on the analysis, we can use a timer for the computation    
    # start_time = time.time()
    # main()
    Nodes = len(G)
    
    Cliques = nx.find_cliques(G)
    Limit = max_cliques(Nodes, kmax)
    Cl=[]
    while True:
        try:
            for i in range(0,Limit):
                clq = next(Cliques)
                if len(clq) <= kmax: # IF YOU DON'T WANNA USE KMAX
                                     # JUST COMMENT THIS
                    Cl.append(clq)
        except StopIteration:
            break
    # If not provided, compute maximal cliques
    # if (C is None) : C = nx.find_cliques(G)
    
    # Sort each clique, make sure it's a tuple
    C = [list(sorted(c)) for c in Cl]
    
    return C

def Kmax_all_cliques(G, kmax=kmax):
    """"
    Enumerate all cliques to a maximum (fixed) size
    """
    C = Kmaxcliques(G)
    Sk = set()
    for k in range(0, max(len(s) for s in C)) :
        # Get all (k+1)-cliques, i.e. k-simplices, from all max cliques mc in C
        # Notice that we are usning set(c) so that we count each clique only once
        [Sk.add(c) for mc in C for c in (itertools.combinations(mc, k+1))]
        # Check that each simplex is in increasing order
        # assert(all((list(s) == sorted(s)) for s in Sk))
        # Assign an ID to each simplex, in lexicographic order
        # S.append(dict(zip(sorted(Sk), range(0, len(Sk)))))
        # Appending the number of cliques of size k+1
        Cliques = [list(i) for i in Sk]
    return Cliques



def Euler_charac(G, kmax=kmax): 
    """
    
    Parameters
    ----------
    G: networkx graph 
    
    kmax: int
        number of dimensions
    
    Returns
    -------
    summary:
        A list with a topological summary for the graph G with Euler 
        characteristics, tau, and number of cliques for each size
    OBS:
    ---
    This function limits the number of cliques to a maximum kmax
    
    """
    
    #start_time = time.time()
    #main()
    Nodes = len(G)
    
    Cliques = nx.find_cliques(G)
    Limit = max_cliques(Nodes, kmax)
    Cl = []
    while True:
        try:
            for i in range(0,Limit):
                clq = next(Cliques)
                if len(clq) <= kmax: # IF YOU DON'T WANNA USE KMAX, 
                                     # JUST COMMENT THIS
                    Cl.append(clq)
        except StopIteration:
            break
    # If not provided, compute maximal cliques
    # if (C is None) : C = nx.find_cliques(G)
    
    # Sort each clique, make sure it's a tuple
    C = [tuple(sorted(c)) for c in Cl]
    
    
    summary = []
    for k in range(0, max(len(s) for s in C)) :
        # Get all (k+1)-cliques, i.e. k-simplices, from all max cliques mc in C
        # Notice that we are usning set(c) so that we count each clique only once
        Sk = set(c for mc in C for c in itertools.combinations(mc, k+1))
        # Check that each simplex is in increasing order
        # assert(all((list(s) == sorted(s)) for s in Sk))
        # Assign an ID to each simplex, in lexicographic order
        # S.append(dict(zip(sorted(Sk), range(0, len(Sk)))))
        # Appending the number of cliques of size k+1
        summary.append(len(Sk))
    tau = sum(summary) # Tau gives the total number of cliques
    kmax = len(summary) # Kmax is the maximum clique size one can find
    ec = 0 # ec is the Euler characteristics
    for i in range(0,len(summary)):
        if i%2 == 0:
                ec += summary[i]
        if i%2 == 1:
                ec += -summary[i]
        #ec+=(-1)**(k % 2)*k
    #print((k))
    summary.insert(0, kmax)
    summary.insert(0, tau)
    summary.insert(0, ec)
    # I want to include new elements after kmax with zero, to say that 
    # there are no simplicies with this size, but all the outputs will 
    # be lists with the same size
    for i in range(kmax, 30): 
        # The first guy is chi, the second is tau, the third is kmax
        summary.insert(kmax+3, 0) 
    
    # The output will be summary starting with EC, tau, kmax, clique_0,
    # Clique_1,Clique_2, Clique_3, and so on...
        
    return summary


def Graph_thresh(e, i): 
    """Creating a binarized graph with a specific threshold 
    
    Parameters
    ----------
    e: int
        threshold value
        
    i: numpy matrix
        connectivity matrix
    
    Returns
    -------
    temp: networkx graph
        
    Notes
    -------
    Returns a graph that maps all elements greater than e to the zero element
    
    """
    
    # Notice that we did not normalize the data. If you want to normalize 
    # just uncomment here
    # ScdatanGA=np.array(normalize(Aut[i]))
    
    data = i
    # be careful to always pass copy of data, othewise will change the data as well
    cpdata = (np.copy(np.abs(data))) 
    cpdata[(np.copy(np.abs(data))) <= (1-e)] = 0.0
   
    thresh_graph= nx.from_numpy_matrix(cpdata[:,:])
    
    return thresh_graph


def densthr(d, i,DIAGNOSTIC=False):
    """Creating a binarized graph with a specific density
    
    Parameters
    ---------   
    d: float
        density value
        
    i: numpy matrix
        connectivity matrix
        
    Returns
    -------
    finaldensity: float
        final density value 
    
    G1: networkx graph
        graph with the specified density
        
    """
    
    np.fill_diagonal(i, 0)
    temp = sorted(i.ravel(), reverse=True) # Will flatten it and rank corr values.
    size = len(i)
    cutoff = np.ceil(d*(size*(size-1)))
    tre = temp[int(cutoff)]
    G0 = nx.from_numpy_matrix(i)
    G0.remove_edges_from(list(nx.selfloop_edges(G0)))
    G1 = nx.from_numpy_matrix(i)
    for u,v,a in G0.edges(data=True):
        if (a.get('weight')) <= tre:
            G1.remove_edge(u, v)
    finaldensity = nx.density(G1)
    if DIAGNOSTIC == True:
        print(finaldensity)
    
    return G1


def Eulerange_thr(i, maxvalue):
    """
    computes the Euler Characteristic and the respective summary metrics 
    for a range of thresholds.
    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
    
    Returns
    -------
    Ec: List with Euler characteristic for a list of thresholds
        
    Notes
    -------
    Filtration process based on thresh
    Notice that we sliced the network in 1/100 steps. We can always do higher or 
    lower steps if required.
    
    """

    Ec = []
    for j in range(0, maxvalue):
        Ec.append(Euler_charac(Graph_thresh(j/100, i)))
        
    return Ec


def Eulerange_dens(i, maxvalue):
    """
    Computes the Euler Characteristic and the respective summary metrics 
    for a range of densities.

    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
    
    Returns
    -------
    Ec: List with the Euler characteristic for a list of densities
        
    Notes
    -------
    Filtration process based on density
    Notice that we sliced the network in 1/100 steps. 

    
    """

    Ec = []
    for j in range(0, maxvalue):
        Ec.append(Euler_charac(densthr(j/100, i)))
        
    return Ec


def Eulerchoice_thr(i, maxvalue, k):
    """
    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
    
    k: int
        euler characteristic=0, total=1,max=2, 3=vertices, 4=edges, 
        5=triangles, etc.
    
    Returns
    -------
    output: Returns a list with  an specific summary metric k for a 
    range of thresholds.
        
    Notes
    -------
    

    """

    temp = Eulerange_thr(i, maxvalue)
    output = [temp[i][k] for i in range(0, maxvalue)]
    
    return output


def Eulerchoice_dens(i, maxvalue, k):
    """
    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
    
    k: int
        euler characteristic=0, total=1, max=2, 3=vertices, 4=edges, 
        5=triangles, etc.
    
    Returns
    -------
    output: Returns a list with  an specific summary metric k for a 
    range of densities.
        
    Notes
    -------
    

    """

    temp = Eulerange_dens(i, maxvalue)
    output = [temp[i][k] for i in range(0, maxvalue)]
    
    return output


def plotEuler_thr(i, maxvalue):
    """Plotting the Euler entropy, i.e. the logarithm of the Euler 
    characteristics for a given threshold interval.
    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
        from 0 to 100
    
    Returns
    -------
    A plot of the Euler entropy based on thr
        
    """
    
    # Change to eulerchoice_dens if intended
    plt.plot(np.log(np.abs(Eulerchoice_thr(i, maxvalue, 0)))) 
    plt.xlabel('Threshold (ε)')
    plt.ylabel('Euler entropy Sχ = ln |χ(ε)|')
    locs, labels = plt.xticks()
    plt.xticks(locs, list(locs/100))
    plt.xlim(0, maxvalue)
    plt.show()
    
    
def plotEuler_den(i, maxvalue):
    """Plotting the Euler entropy, i.e. the logarithm of the Euler 
    characteristics for a given density interval.
    
    Parameters
    ---------    
    i: numpy matrix
        connectivity matrix
    
    maxvalue: int
        from 0 to 100
    
    Returns
    -------
    A plot of the Euler entropy based on density
        
    """
    
    # Change to eulerchoice_dens if intended
    plt.plot(np.log(np.abs(Eulerchoice_dens(i, maxvalue, 0)))) 
    plt.xlabel('Density (d)')
    plt.ylabel('Euler entropy Sχ = ln |χ(d)|')
    locs, labels = plt.xticks()
    plt.xticks(locs, list(locs/100))
    plt.xlim(0, maxvalue)
    plt.show()

    

def Curv_density(d, i, verbose=False):
    """Compute nodal curvature (Knill's curvature) based on density
    
    Parameters
    ---------   
    d: float
        density value
        
    i: numpy matrix
        connectivity matrix
        
    Returns
    -------
    fden: float
        final density value for the graph
    
    curv: numpy array
        array with curvature values
        
    """
    
    def DIAGNOSTIC(*params) :
        if verbose : print(*params)
    DIAGNOSTIC('This function run over all nodes and computes the curvature of the nodes in the graph' )
    
    # This is the initial Graph
    #fden, 
    G_den = densthr(d,i) 
    
    # Enumerating all cliques of G up to a certain size
    temp = Kmax_all_cliques(G_den)
    
    # This lista is a vector V where each v_i is the number of cliques of size i
    lista = []
    
    # We suppose that the size of the cliques are smaller than 50, so we create an 
    # empty list of size 50 for the lista
    for i in G_den.nodes():
        # We start with empty scores for the curvature
        # creating a list of lists for each node - all empty for the 
        # scores for each size for each node
        lista.append([0] * 50) 
    
    DIAGNOSTIC('These are all cliques of the Network:')
    # THIS WILL PRINT ALL THE CLIQUES
    DIAGNOSTIC(temp)
    
    DIAGNOSTIC('We now print the curvature/clique score of each node in the network')
    
    # Now we run over all nodes checking if the node belongs to one clique or another
    Sc = []
    for node in G_den.nodes(): # now we run the script for each clique
        score = 0 # This is the initial score of the node in the participation rank
        for clique in temp:
            # Checking the size of the clique
            k = len(clique)
            # If this node is in the clique, we update the curvature
            if node in clique:
                score += 1 # If the node is in the clique raises the score
                # Increases the curvature score for a size k with a 
                # different weight due to Gauss-Bonnet theorem - is k-1 since 
                # len>0 and python starts from zero.
                lista[node][k-1] += (-1)**(k+1)*1/k 
        Sc.append(score)
        
        DIAGNOSTIC('The node '+str(node)+' has score ='+str(score))
    
    total=[]
    for elements in lista:
        # Summing the participation in all sizes, so that we can compute the 
        # curvature (TOTAL IS ACTUALLY THE CURVATURE - WITHOUT NORMALIZATION)
        total.append(sum(elements)) # This is good if one wants to normalize by the maximum
    DIAGNOSTIC(total)
    
    nor = sum(total) ####!!! not being used //REMOVE ?
    nor2 = max(total) ####!!! not being used //REMOVE ?
    # nt is normalized by the sum
    #nt2 is normalized by the max"
    nt = []
    nt2 = []
    
    # I just removed where one could find division by zero
    #for i in range(0,len(total)):
    #    nt.append(total[i]/nor)
    #    nt2.append(total[i]/nor2)
    most = np.argsort(-np.array(total))#
    
    #def showrank():
    for i in most:
            DIAGNOSTIC('the node ' + str(i)+ ' is in '+ str(total[i])+ ' cliques')
    #    return 
    #DIAGNOSTIC(showrank())
    
    DIAGNOSTIC('These are the most important nodes ranked according to the total clique score')
    DIAGNOSTIC(most)
    DIAGNOSTIC('These is the array nt')

    DIAGNOSTIC(nt)
    DIAGNOSTIC('These is the array nt2')

    DIAGNOSTIC(nt2)
    DIAGNOSTIC('These is the array lista')

    DIAGNOSTIC(lista)
    DIAGNOSTIC('The output is one vector normalizing the value from the maximum')
    #vector=10000*np.array(nt)
    # nor2 is the maximum- The output nt2 is in percentage - 
    # That means the max get 100 and the rest bet 0-100
    
    # curv gives the curvature  - put Sc instead of curv to get that 
    # the particiaption rank - notice that you can normalize in many ways"
    curv = []
    for i in range(0, len(lista)):
        # Summing up for a fixed node all the curvature scores gives the 
        # curvature of the nodes
        curv.append(sum(lista[i]))
        
    curv = np.array(curv)
    # Now, the curvature is not normalized!!!
    return curv#fden, curv


def Curv_thr(e, i, verbose=False):
    """Compute nodal curvature based on threshold
    
    Parameters
    ---------   
    e: float
        threshold value
        
    i: numpy matrix
        connectivity matrix
        
    Returns
    -------
    curv: numpy array
        array with curvature values
        
    """
    
    def DIAGNOSTIC(*params):
        if verbose : print(*params)
    DIAGNOSTIC('This function run over all nodes and computes the curvature of the nodes in the graph' )
    
    # This is the initial Graph
    G_thr = Graph_thresh(e,i)  
    temp = Kmax_all_cliques(G_thr)
    
    # This lista is a vector V where each v_i is the number of cliques of size i
    lista = []
    
    # We suppose that the size of the cliques are smaller than 20, 
    # so we create an empty list of size 20 for the lista
    for i in G_thr.nodes():
        lista.append([0] * 50) # creating a list of lists for each node - 
        # all empty for the scores for each size for each node
    
    DIAGNOSTIC('These are all cliques of the Network:')
    DIAGNOSTIC(temp)
    DIAGNOSTIC('We now print the curvature/clique score of each node in the network')
    
    # Now we run over all nodes checking if the node belongs to one clique or another
    Sc = []
    for node in G_thr.nodes(): # now we process for each clique
        score = 0 # This is the initial score of the node in the participation rank
        for clique in temp:
            k = len(clique)
            if node in clique:
                score += 1 # If the node is in the clique raises the score
                lista[node][k-1] += (-1)**(k+1)*1/k # Increases the curvature score 
                                                    #  for a size k with a different 
                                                    # weight due to Gauss-Bonnet theorem
        Sc.append(score)
        
        DIAGNOSTIC('The node '+str(node)+' has score ='+str(score))
    
    total=[]
    for elements in lista:
        total.append(sum(elements)) # This is good if one wants to normalize by the maximum
    DIAGNOSTIC(total)
    
    nor = sum(total) #### NOTICE THAT nor AND nor2 IS NOT BEING USED 
    #                   - HOWEVER, IF ONES DECIDES FOR NORMALIZATION- 
    #                   SHOULD TAKE THIS INTO ACCOUNT
    nor2 = max(total) 
    # nt is normalized by the sum
    #nt2 is normalized by the max"
    nt = []
    nt2 = []
    
    # I just removed where one could find division by zero
    #for i in range(0,len(total)):
    #    nt.append(total[i]/nor)
    #    nt2.append(total[i]/nor2)
    most=np.argsort(-np.array(total))#
    
    #def showrank():
    for i in most:
            DIAGNOSTIC('the node ' + str(i)+ ' is in '+ str(total[i])+ ' cliques')
    #    return 
    #DIAGNOSTIC(showrank())
    
    DIAGNOSTIC('These are the most important nodes ranked according to the total clique score')
    DIAGNOSTIC(most)
    DIAGNOSTIC('These is the array nt')

    DIAGNOSTIC(nt)
    DIAGNOSTIC('These is the array nt2')

    DIAGNOSTIC(nt2)
    DIAGNOSTIC('These is the array lista')

    DIAGNOSTIC(lista)
    DIAGNOSTIC('The output is one vector normalizing the value from the maximum')
    #vector=10000*np.array(nt)
    # nor2 is the maximum- The output nt2 is in percentage - 
    # That means the max get 100 and the rest bet 0-100
    
    # curv gives the curvature  - put Sc instead of curv to get that the 
    # particiaption rank - notice that you can normalize in many ways
    curv = []
    for i in range(0, len(lista)):
        curv.append(sum(lista[i]))# Summing up for a fixed node all the 
                                  #  curvature scores gives the curvature of the nodes
    curv = np.array(curv)
    
    return curv


def SaveEuler(individual, name, tresh):
    """Save Euler results
    
    Parameters
    ---------           
    individual: numpy matrix
        connectivity matrix
        
    name: str
        file name
        
    tresh: float
        threshold value
        
    Returns
    -------
    Files with results
        
    """
    
    values = (Eulerchoice_thr(individual,tresh,0)) # change to eulerchoice_dens if intended
    with open(name, 'w') as output:
        output.write(str(values))
        
        
        
def Participation_in_cliques(d, i, cl, verbose=False):
    """
    Returns a list with the participation rank in cliques of a fixed size
    inputs:
    d: density
    i: matrix
    cl: clique size
    """
    def DIAGNOSTIC(*params) :
        if verbose : print(*params)
        return

    # I want that the output is a vector analogous with the previous one, 
    # but for a fixed k, not for all k
    # COMPUTING THE CLIQUES
    G = densthr(d,i)  
    temp = Kmax_all_cliques(G) 
    DIAGNOSTIC('These are all cliques')
    DIAGNOSTIC(temp)
    "This lista is a vector V where each v_i is the number of cliques of size i"
    lista = []
    # We suppose that the size of the cliques are smaller than 50, 
    # so we create an empty list of size 50 for the lista
    for i in G.nodes():
        lista.append([0] * 50) 
    #Creating a list of lists - for the particiaption scores of the nodes
    "Now we run over all nodes checking if the node is in one clique or another"
    for node in G.nodes(): # now I have to do the process for is in clique
        score = 0 # This is the score of the node
        # running over all nodes in G
        for clique in temp:
            #running over all cliques enumerated before
            k = len(clique)
            #checking whether there is a node in the clique of size k 
            if node in clique:
                #INCLUDING THE SCORE FOR THE CLIQUE
                score += 1
                lista[node][k-1] += 1
       # print('the node '+str(node)+' has score ='+str(score))
    total = []
    for elements in lista:
        total.append(sum(elements))
    DIAGNOSTIC('This is the number of cliques each node is participating in')
    DIAGNOSTIC(total)
    DIAGNOSTIC(np.sum(total))
    nor = sum(total)
    nt = []
    for i in range(0,len(total)):
        nt.append(total[i]/nor)
    
    klist = []
    DIAGNOSTIC('Now lets plot the number of k-cliques each number is participating in')
    for i in G.nodes():
        klist.append(lista[i][cl-1])
        DIAGNOSTIC('the node '+str(i)+ ' has '+ str(cl) + ' - score =' + str(lista[i][cl-1]))
    
    mostk = np.argsort(-np.array(klist))
    
    
    DIAGNOSTIC('These are the most important nodes ranked according to the k-clique score')

    DIAGNOSTIC(mostk)
    for i in mostk:
            DIAGNOSTIC('the node ' +str(i)+ ' is in '+ str(klist[i])+ ' ' +str(cl)+ '-cliques')
    DIAGNOSTIC(klist)
    
    
    maxk = max(klist)
    totk = sum(klist)
    # We can do several choices on the particiaption rank: here we 
    # choose the percentage of all k-cliques
    return np.nan_to_num(100*np.array(klist)/totk)

def Betti_k(G, K_input, verbose=False):
    """#Function to compute the desired Betti number of a network 

    Parameters
    ----------
    Argument1: G is a network. It is in the shape of a network graph of 
               the library networkx G.graph()
    Argument2: K_input--> insert 1 if you want to compute Betti-1, 2 for
               Betti-2 etc.
    
    Returns
    -------
    out: value of the desired Betti number.

    """
    def DIAGNOSTIC(*params): # If verbose is True it will print all the DIAGNOSTIC
        if verbose: 
            print(*params)

    DIAGNOSTIC("Nodes in G: ", G.nodes())    
    DIAGNOSTIC("Edges in G: ", G.edges())
    DIAGNOSTIC("Number of nodes: {}, edges: {}".format(G.number_of_nodes(), G.number_of_edges()))
    
    # 1. Prepare maximal cliques
    
    # compute maximal cliques
     
    C = nx.find_cliques(G) # C now is the operator "find clique" 
                           # (to do the list I should do list(nx.find_cliques(G)))
    
    #Create list C with all the cliques
    #Sort each clique, convert it from list to tuple
    C = [tuple(sorted(c)) for c in C]
    DIAGNOSTIC("List with all maximal simplex C:",C)    
    DIAGNOSTIC("Number of maximal cliques: %i"%(len(C)))    
   
    # 2. Enumerate all simplices
    
    S = [] #List of dictionaries
    # S[k] is the dictionary which contain all k-simplices
    #S[k].keys() are simplex s (s is one of the k-simplex of the dictionary S[k])
    # S[k].values() are the ID of simplex s
    DIAGNOSTIC("I start the loop where I create the required Sk to then compute betti. Sk is a list with the k-simplex")
    #I set the range for the following loop
    if K_input==0:
        ini = 0
        fin = 2
    else:
        ini = K_input-1
        fin = K_input+2
        
    for k in range(ini,fin) : # k has 2 values for betti_0 and 3 values for betti1_2_3 
        
        Sk = sorted(set(c for mc in C for c in itertools.combinations(mc, k+1)))
       
        DIAGNOSTIC("list of %i-simplex S%i:"%(k,k), Sk)        
        # Check that each simplex is in increasing order
        assert(all((list(s) == sorted(s)) for s in Sk))        
        # Assign an ID to each simplex, in order
        S.append(dict(zip(Sk, range(0, len(Sk))))) # zip(Sk,range()) is 
                                                   # an object (composed by tuples) 
                                                   # where each element of Sk is 
                                                   # associated to a number.
                                                   # I then from the zip object 
                                                   # create the dictionary 
                                                   # where the key is the Sk element 
                                                   # and the value the number
                                                   # I put this dictionary in 
                                                   # the S list (list of dictionary)   
        DIAGNOSTIC("Number of %i-simplices: "%(k),len(Sk))    
    DIAGNOSTIC("S dictionary",S)      
    # The cliques are redundant now
    del C
    
    # 3. Construct the boundary operator
    
    # Boundary Matrix 
    D = [None, None] # List with the two different k-boundary operators
    
    if K_input==0:
        # D[0] is the zero matrix
        D[0]=(np.zeros((1, G.number_of_nodes())))# I create a matrix of size (1,#nodes)
        
    for k in range(1, len(S)):
        
        # I set the index of D[] and the number of nodes in each group for the combinatory part
        if K_input==0:
            index = k
            b = k
        else:
            index = k-1
            b = k+(K_input-1)
            
        D[index] = np.zeros( (len(S[k-1]), len(S[k])) ) # I create a matrix of size 
                                                        # (len(S[k-1]), len(S[k])

        for (ks, j) in S[k].items() :
            
            a = sorted(itertools.combinations(ks, b))
            #print("a",a)
            # Indices of all (k-1)-subsimplices s of the k-simplex ks
            # S is a list of dictionary with k different size 
            I = [S[k-1][s] for s in sorted(itertools.combinations(ks, b))] 
            #print("I",I)

            for i in range(0,len(I)):
                D[index][I[i]][j] = (-1)**(i)

        if D[index].shape[1]==0:
            DIAGNOSTIC("I can't create matrix D because I don't have the needed k-simplex")
        
        DIAGNOSTIC("Boundary matrix:")
        DIAGNOSTIC("D",D[index])
            
    DIAGNOSTIC("D_{} has shape {}".format(K_input, D[0].shape))
    DIAGNOSTIC("D_{} has shape {}".format(K_input+1, D[1].shape))
    # The simplices are redundant now
    del S
    
    # 4. Compute rank and dimker of the boundary operators
    
    # dim(Im)=Rank and dim(ker)=V-rank
    rank = [0 if d.shape[1]==0 else np.linalg.matrix_rank(d) for d in D] #dim(Im)
    ker = [(d.shape[1] - rank[n]) for (n, d) in enumerate(D)] #V - rank = dim(ker) ,rank=dim(Im)
    
    #The boundary operators are redundant now
    del D
    DIAGNOSTIC("ker:", ker)
    DIAGNOSTIC("rank:", rank)
    
    # 5. Compute the Betti number   
    
    # Betti number
    DIAGNOSTIC("Betti= ker[0]-rank[1]")
    B=ker[0]-rank[1]
    DIAGNOSTIC("End of computation\nBetti %i is:"%K_input,B)
    return B