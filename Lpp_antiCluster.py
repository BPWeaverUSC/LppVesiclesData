import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rand
import itertools as it

from scipy.spatial import Voronoi, voronoi_plot_2d

############Begin Functions#################

def randpairs_cluster(N,R,cutoff):
    pairs = np.zeros((N,2))
    dif_x = np.zeros((N,1))
    dif_y = np.zeros((N,1))
    d_list = np.zeros((N,1))
    
    #First point and counter
    pairs[0][0] = R*rand.random()
    pairs[0][1] = R*rand.random()
    points = 1
    
    #BEGIN LOOP HERE
    while (points < N):
        #Candidate location
        px = R*rand.random()
        py = R*rand.random()
        
        #calculate distance to all points assigned so far
        for x in range(points):
            dif_x[x] = pairs[x][0] - px
            dif_y[x] = pairs[x][1] - py
            d_list[x] = dist(dif_x[x],dif_y[x])
        
        #Find the closest point already placed
        d_toSearch = d_list[:points]
        min_d_indx = np.argmin(d_toSearch)
        d = d_list[min_d_indx]
        
        #If the candidate point is within the cutoff range, change the location based on the generated random number.
        if (d < cutoff):
            alpha = cutoff/(d)
            
            X = (alpha-1.)*dif_x[min_d_indx]
            Y = (alpha-1.)*dif_y[min_d_indx]
            
            # In this version, points are moved out to the cutoff distance, instead of closer.
            px_new = -X + px #pairs[min_d_indx][0]
            py_new = -Y + py #pairs[min_d_indx][1]
        else:
            px_new = px
            py_new = py
        
        pairs[points][0] = px_new
        pairs[points][1] = py_new
        
        points = points + 1
        #END LOOP
        
    return pairs

def dist(modx,mody):
    d = ((modx)*(modx)+(mody)*(mody))**(0.5)

    return d

def pbc(mod,L):
    dif = mod

    if (dif >= L/2.0):
        dif -= L
    elif(dif <= -L/2.0):
        dif += L
    else:
        dif = dif

    return dif

#############End Functions####################

################################Main Body##########################################

size_r = 100
lpp = 140
pal = 30
ompA = 14

#if (lpp > 2*pal):
#    size = lpp - pal
#elif (lpp <= 2*pal):
#    size = pal

size = lpp + pal + ompA

###### In here, instead of assigning all of the locations randomly, points are forced to a minimum distance
cutoff = 10
p = randpairs_cluster(size,size_r,cutoff)

######
Pnew = np.zeros((size,2))

X = np.zeros(size)
Y = np.zeros(size)


for x in range(size):
    Pnew[x][0] = p[x][0]
    Pnew[x][1] = p[x][1]



for x in range(size):
    X[x] = Pnew[x][0]
    Y[x] = Pnew[x][1]

vor = Voronoi(Pnew)

vert = vor.vertices

Del = []

for x in range(len(vert)):
    if (vert[x][0] > size_r):
        Del.append(x)
    elif (vert[x][0] < 0):
        Del.append(x)
    elif (vert[x][1] > size_r):
        Del.append(x)
    elif (vert[x][1] < 0):
        Del.append(x)

vert = np.delete(vert, Del, axis=0)


circ_cen = np.zeros((len(vert),2))


dlist = np.zeros(( len(vert),(len(Pnew)+2) ))
distances = np.zeros(( len(vert),len(Pnew) ))
radii = np.zeros(len(vert))
ave_dist = np.zeros(len(vert))
find_rad = np.zeros(len(Pnew))

for x in range(len(vert)):
    dlist[x][0] = vert[x][0]
    dlist[x][1] = vert[x][1]
    
    for y in range(len(Pnew)):
        dif_x = vert[x][0] - X[y]
        dif_y = vert[x][1] - Y[y]
        
        #dif_x = pbc(dif_x,(size_r))
        #dif_y = pbc(dif_y,(size_r))
        
        #dlist[x][y+2] = dist(dif_x,dif_y)
        distances[x][y] = dist(dif_x,dif_y)
        #sum = dlist[x][y+2]
        #dlist[x][len(Pnew)+2] += (sum)/(len(Pnew))


distances.sort()


for x in range(len(distances)):
    radii[x] = distances[x][0]

Del = []

for x in range(len(vert)):
    if (vert[x][0] + radii[x] > size_r):
        Del.append(x)
    elif (vert[x][0] - radii[x] < 0):
        Del.append(x)
    elif (vert[x][1] + radii[x] > size_r):
        Del.append(x)
    elif (vert[x][1] - radii[x] < 0):
        Del.append(x)

vert = np.delete(vert, Del, axis=0)
radii = np.delete(radii,Del, axis=0)



sradii = np.sort(radii)


#################################################################################

#######################################################
circle = [None]*len(vert)
#################Plotting##############################
##

#mpl.rc('text', usetex = True)
plt.ylim(-0.5,size_r+0.5)
plt.xlim(-0.5,size_r+0.5)
##
plt.xlabel(r'$\mathrm{Position (nm)}$')
plt.ylabel(r'$\mathrm{Position (nm)}$')
##
##plt.title("All circular regions")
##plt.plot(vert,'x',)
##
for x in range(len(vert)):
    circle[x] = plt.Circle( (vert[x][0],vert[x][1]), radii[x], color='r', alpha=0.1)
##
fig, ax = plt.gcf(),plt.gca()
##
for x in range(len(vert)):
    ax.add_patch(circle[x])
##
ax.set_aspect('equal', adjustable='datalim')
##
####Use adjustable='box-forced' to make the plot area square-shaped as well.
ax.plot(X,Y,'.') #Uncomment for all circle plot
##ax.plot([0,0],[0,size_r],'g-')
##ax.plot([0,size_r],[0,0],'g-')
##ax.plot([0,size_r],[size_r,size_r],'g-')
##ax.plot([size_r,size_r],[size_r,0],'g-')
##
###voronoi_plot_2d(vor)
##
plt.rcParams.update({'font.size': 18})
#plt.show()
plt.savefig("Clustered.pdf", format="pdf")
#####################################################
