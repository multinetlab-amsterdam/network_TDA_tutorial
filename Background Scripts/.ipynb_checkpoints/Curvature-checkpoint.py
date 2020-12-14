import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import preprocessing
import networkx as nx
import itertools

"This is the normalizaton for by taking the MaxAbs Scaler"

def normalize(matrix):# This is now the one we are going to use
    
    X_scale=np.copy(matrix) # This is the data I want to scale
    max_abs_scaler = preprocessing.MaxAbsScaler()# This is the one I can use for the HCP
    X_scale_maxabs = max_abs_scaler.fit_transform(X_scale)
   
    return X_scale_maxabs #X_train_minmax



"This is the main function"
def Euler_charac(G):    #
    # 1. Prepare maximal cliques
    
    #start_time = time.time()
    #main()
    
    C = nx.find_cliques(G)
    # If not provided, compute maximal cliques
    #if (C is None) : C = nx.find_cliques(G)
    
    # Sort each clique, make sure it's a tuple
    C = [tuple(sorted(c)) for c in C]
    
    
    S2= []
    for k in range(0, max(len(s) for s in C)) :
        # Get all (k+1)-cliques, i.e. k-simplices, from max cliques mc
        Sk = set(c for mc in C for c in itertools.combinations(mc, k+1))
        # Check that each simplex is in increasing order
        #assert(all((list(s) == sorted(s)) for s in Sk))
        # Assign an ID to each simplex, in lexicographic order
        #S.append(dict(zip(sorted(Sk), range(0, len(Sk)))))
        S2.append(len(Sk))
    tau=sum(S2) # Tau gives the total number of cliques
    kmax=len(S2) # Kmax is the maximum clique size one can find
    ec=0 # ec is the Euler characteristics
    for i in range(0,len(S2)):
        if i%2==0:
                ec+=S2[i]
        if i%2==1:
                ec+=-S2[i]
        #ec+=(-1)**(k % 2)*k
    #print((k))
    S2.insert(0,kmax)
    S2.insert(0,tau)
    S2.insert(0,ec)
    for i in range(kmax,30): # I want to include new elements after kmax with zero, to say that there are no simplicies with this size
        S2.insert(kmax+3,0) # The first guy is chi, the second is tau, the third is kmax
    # The output will be EC, tau, knax, clique_0,Clique_1,Clique_2, Clique_3, and so on...
    return S2




"Here we create the filtration - e is the threshold and i is the number of the individual"
def Graph_thresh(e,i): 
    
#this makes a graph from one slice only!!! Later, you have to ask the computer to make this many times
    #"Returns a graph that map all elements greater than e to the zero element"
    
    #"If you want to normalize just uncomment here"
    #ScdatanGA=np.array(normalize(Aut[i]))
    
    ScdatanGA=i
    
    
    teste1= (np.copy(np.abs(ScdatanGA))) # be careful to always assing  copy of data, othewise will change the data as well
    teste1[(np.copy(np.abs(ScdatanGA))) <= (1-e)] = 0.0
   
    temp=nx.from_numpy_matrix(teste1[:,:])
    
    return temp


"I'm doing here 100 slices in the interval 0-1, you can also go a bit more if necessary"
"i is the patient"
def Eulerange(i,max):
    Ec=[]
    for j in range(0,max):
        Ec.append(Euler_charac(Graph_thresh(j/100,i)))
    return(Ec)


"To see anotehr one parameter, just print [thresh][param], euler=0,total=1,max=2,3=vertices,4=edges,5=triangles, etc. etc. "

def Eulerchoice(i,max,k):
    "k=0 - means the euler characteristics- i is the individual"
    temp=Eulerange(i,max)
    output=[temp[i][k] for i in range(0,max)]
    return(output)

    
"if you want to plot directly copy and past the command bellow"

def plotEuler(i,max,k): #aquelas linhas, max tem que ser de 0 100
    plt.plot(np.log(np.abs(Eulerchoice(i,max,0))))
    plt.xlabel('Threshold (ε)')
    plt.ylabel('Euler entropy Sχ = ln |χ(ε)|')
    locs, labels = plt.xticks()
    plt.xticks(locs, list(locs/100))
    plt.xlim(0,max)
    plt.show()

    return



def Curv(e,i,verbose=False): # histogramas  
    def DIAGNOSTIC(*params) :
        if verbose : print(*params)
    DIAGNOSTIC('This function run over all nodes and computes the curvature of the nodes in the graph' )
    G = Graph_thresh(e,i)  # This is the initial Graph
    temp=list(nx.enumerate_all_cliques(G))
    "This lista is a vector V where each v_i is the number of cliques of size i"
    lista=[]
    "We suppose that the size of the cliques are smaller than 20, so we create an empty list of size 20 for the lista"
    for i in G.nodes():
        lista.append([0] * 20) # creating a list of lists for each node - all empty for the scores for each size for each node
    DIAGNOSTIC('These are all cliques of the Network:')
    DIAGNOSTIC(temp)
    "Now we run over all nodes checking if the node belongs to one clique or another"
    DIAGNOSTIC('We now print the curvature/clique score of each node in the network')
    Sc=[]
    for node in G.nodes(): # now I have to do the process for each clique
        score=0 # This is the initial score of the node in the participation rank
        for clique in temp:
            k=len(clique)
            if node in clique:
                score+=1 # If the node is in the clique raises the score
                lista[node][k-1]+=(-1)**(k+1)*1/k # Increases the curvature score for a size k with a different weight due to Gauss-Bonnet theorem
        DIAGNOSTIC('The node '+str(node)+' has score ='+str(score))
        Sc.append(score)
    total=[]
    for elements in lista:
        total.append(sum(elements)) # This is good if one wants to normalize by the maximum
    DIAGNOSTIC(total)
    nor=sum(total)
    nor2=max(total)
    "nt is normalized by the sum"
    "nt2 is normalized by the max"
    nt=[]
    nt2=[]
    # I just removed where one could find division by zero
    #for i in range(0,len(total)):
    #    nt.append(total[i]/nor)
    #    nt2.append(total[i]/nor2)
    most=np.argsort(-np.array(total))#
    
    #def showrank():
    for i in most:
            DIAGNOSTIC('the node ' +str(i)+ ' is in '+ str(total[i])+ ' cliques')
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
    " nor2 is the maximum- The output nt2 is in percentage - That means the max get 100 and the rest bet 0-100"
    curv=[]
    for i in range(0,len(lista)):
        curv.append(sum(lista[i]))# Summing up for a fixed node all the curvature scores gives the curvature of the nodes
    "curv gives the curvature  - put Sc instead of curv to get that the particiaption rank - notice that you can normalize in many ways"
    return np.array(curv)


"You can save one plot of the Euler characteristic here - choose the name of the file, the treshold and which individual you want"
"notice that for higher threshold this might not converge"

"Consider including a timer in the Euler algorithm"
def SaveEuler(individual,name,tresh):
    values=(Eulerchoice(individual,tresh,0))
    with open(name, "w") as output:
        output.write(str(values))