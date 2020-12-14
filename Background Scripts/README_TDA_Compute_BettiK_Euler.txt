The file "TDA_Compute_BettiK_Euler" contains Topological data analysis methods to analyse connectivity matrix.

The code is divided in two different part:
1) Computation of Betti-0 with an extreme precise algorithm.
   Visualisation of data.

2) Computation of Betti-k (also higher order Betti) with an algorithm conceptually different from the previous one.
   Computation of Euler characteristic using two different methods.
   Examples of application of the methods 
   Visualisation of the data
   
      
1) 
GOAL OF THE CODE:
Starting from connectivity matrix CX using an algorith we create the Dx matrices. 
In the Dx matrices all the information of the different clusters are stored.
The Betti-0 curves (# connected components vs filtration value) created with the following code are very precise, 
all the filtration value in the matrix are taken into account
We also illustrate how to visualize the Betti-0 curves to see how the connencted components are merging with the increasment
of the filtration value and we plot the results.

HOW TO USE THE CODE: 

You just need the matrix Cx you want to convert. It should be an Array type.
Then you just call the function "from_cx_to_dx":

Dx,a,conn= from_cx_to_dx(Your_Matrix)

The output are:
- the converted matrix Dx
- a= list containing the filtration values where changes in the network occur. 
- conn= list containing the number of connected components for each filtration value. 

This is all you need to create the new matrix and run the code!! 
You also now have all the information to plot the Betti-0 curves! 
(fig=plt.figure(figsize=(15, 10))
plt.plot(a,conn, marker='.'))

MORE DETAILED USE: 
If you want to run the code on several Cx matrices you can follow the example contained in the code (see code "#Code to create a big object with all the Cx matices"): 
txt_files_cx = glob.glob('Cx_matrix_YOUR_FOLDER/cx*.txt') HERE you put the folder containing your Cx matrices .txt files. If the name of your matrices files are cx_1,cx_2,cx_3... you just put cx* to take all of them.

Then you just need to run the part with the for loop (see code).
JUST REPLACE in these lines the desired folders:
g=g.replace('Cx_matrix_YOUR_FOLDER\\cx','Dx_matrix_YOUR_NEW_FOLDER\\dx')#You have to create an empty folder Dx_matrix_YOUR_NEW_FOLDER to run this part of the code.

l=l.replace('Cx_matrix_YOUR_FOLDER\\','Beta_YOUR_NEW_FOLDER\\') #You have to create an empty folder Beta_YOUR_NEW_FOLDER to run this part of the code. Here in each file (one for each initial Cx matrix) you store the filtration values list and correspondant connected components list)

In the code you have a visualisation of this part as well!

N.B.To run my examples and run the code as it is without changing anything you need the folder "Cx_matrix_var01" containg  my 10 Cx matrix. And you need to create two empty folders named "Dx_matrix_var01" and "Beta_var01".

2) 
GOAL OF THE CODE:
"""Computation of Betti numbers and Euler characteristic in a network"""

"""Also with this code we can compute Betti-0 but in a less precise way compared to the previous code, here we are not looking at ALL filtration values but slices.
Here we compute not just Betti-0, with the following code we can compute Betti-1,Betti-2....Betti-n and the Euler characteristic.
B_0 is the number of connected components
B_1 is the number of 2D holes
B_2 is the number of 3D holes 
B_3 is the number of 4D holes 
Euler characterisitc : alternative sum of the number of k-simplex
We compute the Euler characteristic in two different ways: 
-with an exact computation (taking into account all k-simplex)
-with an approximation (fixing the maximum dimension of the k-simplex we want to consider).
"""


HOW TO USE THE CODE:

EULER:
To compute the Euler characteristic of a network you just need G (network graph of the library networkx G.graph() ).
You call the Euler function "euler":

E=euler(G,verbose=True)
The output E is the euler characteristic.
 
If you want to use an approximation fixing the maximum number of k-simplex you want to consider you use the function "euler_k":

The imputs are G and the maximum number (e.g. 4)
S=euler_k(G,4,verbose=True)
The output S is a list S= EC, tau, kmax, clique_0,Clique_1,Clique_2, Clique_3...
the Euler characteristic is the first value of the list S[0]

BETTI NUMBERS:
To compute the Betti number use the function "Betti_k":
The imput of the function is the network G and the number of Betti you want to compute (e.g. 0 or 1 or 2...)
Betti_number= Betti_k(G,1,verbose=False)
The otuput is the number of kD-holes (Betti number) for that given network.

In the code you can find several examples of these computations on simple networks.

If you want to generate networks starting from your brain connectivity matrix:
THERE ARE TWO METHODS TO THRESHOLD A NETWORK when you have the connectivity matrix Cx: 
#METHOD 1). USING DIRECTY A VALUE OF a desired FILTRATION OR  
#METHOD 2). CHOOSE THE DENSITY OF CONNECTION YOU WANT TO CONSIDER. (e.g. 10% of strongest connections)
In the code both methods are illustrated. We also illustrated 1 application of method 1) and 2 application of method 2)

METHOD 1):
Use the function "filtr_thresh":

The imputs are--> a list Cx containing all the Cx matrices of the patients, the value of filtration you want to use, the index of patient you want to analys (his connectivity matrix)
N.B. In the code I illustrated how to create the Cx list containing all the Cx matrices of the different patients. (See code)
txt_files_cx = glob.glob('Cx_matrix_YOUR FOLDER/cx*.txt') # In Cx_matrix_YOUR FOLDER there are your .txt files, one for each matrix

G1,filtr= dens_thresh(Cx,0.1,3,verbose=False) # Call the function
The output are the new thresholded network G1 and the filtration value you used. 

Compute the Betti-1 number for that network:
B=Betti_k(G1,1, verbose=False)

DETAILED ANALYSIS of application of method 1:     
You can call the function "filtr_thresh" several times using many different values of filtrations so you obtain many networks. You can do it for each patient. You can store in one file (one for each patient) the number of Betti-k/ Euler you have for each created network.
In the code there is an illustration of this procedure. You need to create an empty folder "Betti_greater" to run that part of code.

N.B. To run my examples and run the code as it is without changing anything you need the folder "Cx_matrix_var01" containg  my 10 Cx matrix. And you need to create an empty folder named "Betti_greater".

METHOD 2):
Use the function "dens_thresh":

The imputs are--> a list Cx containing all the Cx matrices of the patients, the value of density you want to use (e.g 10% stronger connections), the index of patient you want to analys (his connectivity matrix)
N.B. In the code I illustrated how to create the Cx list containing all the Cx matrices of the different patients. (See code)
txt_files_cx = glob.glob('Cx_matrix_YOUR FOLDER/cx*.txt') # In Cx_matrix_YOUR FOLDER there are your .txt files, one for each matrix

G_fi,filtr_value= dens_thresh(Cx,0.1,3,verbose=False) # Call the function
The output are the new thresholded network G_fi and the filtration value correspondant to the desired density. 

Compute the Betti-1 number for that network:
B=Betti_k(G_fi,1, verbose=False)
     
DETAILED ANALYSIS of application of method 2:
You can call the function "dens_thresh" several times using many different values of density so you obtain many different networks. You can do it for each patient. You can store in one file (one for each patient) the number of Betti-k/ Euler you have for each created network.
In the code there is an illustration of this procedure. You need to create an empty folder "Betti_greater" to run that part of code.
In the code there is also the DATA VISUALISATION
N.B. To run my examples and run the code as it is without changing anything you need the folder "Cx_matrix_var01" containg  my 10 Cx matrix. And you need to create an empty folder named "Betti_greater" .

Also another example of application of method 2 is then illustrated:
-Creating a network with 25 nodes and increasind the density of connections and computin each time the desired quantities (Betti-1,Betti_2...)

Lastly there is an example of creation of random networks and comptation of the desired quantities.
To run it you need to create a empty folder "random".
