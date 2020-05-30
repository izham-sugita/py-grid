import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#read NACA0012 data file.
#The data file is generated from naca-C-grid-init.py
#The file naming format is:
# nacaXXXX-XXX-XXX-XXX-XXX-input.csv
# 1. The first XXXX : NACA 4-series airfoil number
# 2. The second XXX: The total nodes for xi-axis
# 3. The third XXX: The total nodes for eta-axis
# 4. The fourth XXX: The index number of the node at tail edge 1 : under
# 5. The fifth XXX: The index number of the node at tail edge 2 : over

filename = "naca0012-101-101-030-070-input.csv"
dat = pd.read_csv(filename)
naca = pd.DataFrame(dat)
airfoil = naca.to_numpy()

imax = filename[9:12]
jmax = filename[13:16]
itail_1 = filename[17:20]
itail_2 = filename[21:24]

imax = imax.lstrip("0")
jmax = jmax.lstrip("0")
itail_1 = itail_1.lstrip("0")
itail_2 = itail_2.lstrip("0")

imax = int(imax)
jmax = int(jmax)
itail_1 = int(itail_1)
itail_2 = int(itail_2)

print(imax, jmax, itail_1, itail_2)

#for drawing reference
all_nodes = len(airfoil)
print("Total nodes for the C-grid:")
print(all_nodes,"\n")


xg = np.ndarray((imax*jmax))
yg = np.ndarray((imax*jmax))

for i in range(imax*jmax):
    xg[i] = airfoil[i][0]
    yg[i] = airfoil[i][1]


xg = xg.reshape((imax, jmax))
yg = yg.reshape((imax, jmax))

#inner boundary plot    
inner_x = np.ndarray((imax))
inner_y = np.ndarray((imax))
for i in range(imax):
    j=0
    inner_x[i] = xg[i][j]
    inner_y[i] = yg[i][j]

#outer boundary plot
outer_x = np.ndarray((imax))
outer_y = np.ndarray((imax))

for i in range(imax):
    j = jmax-1
    outer_x[i] = xg[i][j]
    outer_y[i] = yg[i][j]

#C-type boundary input is the simplest
#not a very good initial distribution
xgn = np.ndarray((imax, jmax))
ygn = np.ndarray((imax, jmax))

xgn[:][:] = xg[:][:]
ygn[:][:] = yg[:][:]

zg = np.zeros_like(xg)

#check sanity
tail_edge = []
for i in range(imax):
    j = 0
    if xgn[i][j] == 1.0 and yg[i][j] == 0.0:
        tail_edge.append(i)

#Solving parabolic type equation 
eps = 1.0e-6
tol = 1.0e-6
itermax = 5000
iter = 0
steps = 100

#Heuristic; but >0.1 is the common practice for numerical stability.
A = 0.001
B = 0.001

tail1 = tail_edge[0]
tail2 = tail_edge[1]
print(tail1, tail2)

# xcenter moved to 0.5 because airfoil is 0.0<= c <= 1.0
xcenter = 0.5
ycenter = 0.0


while iter<itermax:

    # all section combine
    for i in range(0, imax-1):
        for j in range(0, jmax-2):
            if i == 0:
                x_eta2 = (xg[i+1][j]  - 2.0*xg[i][j]  + xg[0][j])
                y_eta2 = ( yg[i+1][j] - 2.0*yg[i][j]  + yg[0][j] )
            else:
                x_eta2 = (xg[i+1][j] -2.0*xg[i][j]  + xg[i-1][j]) 
                y_eta2 = ( yg[i+1][j] -2.0*yg[i][j]  + yg[i-1][j] )

            if j == 0:
                x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][0] - xg[i-1][j+1] + xg[i-1][0] )                
                y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][0] - yg[i-1][j+1] + yg[i-1][0] )
            else:
                x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][j-1] - xg[i-1][j+1] + xg[i-1][j-1] )
                y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][j-1] - yg[i-1][j+1] + yg[i-1][j-1] )

            #Calculate ratio
            # rad: based curvature length from surface to outer boundary
            r1 = np.sqrt( (xg[i][jmax-1] - xg[i][0] )**2 + (yg[i][jmax-1] - yg[i][0])**2  )
            r0 = np.sqrt( (xg[i][0] - xcenter )**2 + (yg[i][0] - ycenter )**2  )
            rvar = (j+1)*(r1 - r0)/(jmax-1) + r0
            ratio = np.log( (rvar+eps) / r0 ) / np.log( (r1+eps)/r0 )

            ratio = 1.0
            #Parabolic solver
            Sx = ratio*(xg[i][jmax-1] - xg[i][j] ) / float( (jmax-1) - j )
            xgn[i][j+1] = xg[i][j] + A*(x_eta2) - 2.0*B*x_reta + Sx

            Sy = ratio*(yg[i][jmax-1] - yg[i][j] ) / float( (jmax-1) - j )
            ygn[i][j+1] = yg[i][j] + A*(y_eta2) - 2.0*B*y_reta + Sy
            

    #Update periodic boundary condition
    for j in range(jmax):
        i = imax-1
        xgn[i][j] = xgn[0][j]
        ygn[i][j] = -ygn[0][j]
        
            
    if iter> 10:
        errx = 0.0
        erry = 0.0
        nodes = 0
        for i in range(1, imax-1):
            for j in range(1, jmax-1):
                errx += (xg[i][j] - xgn[i][j])**2
                erry += (yg[i][j] - ygn[i][j])**2
                nodes +=1

        errx = np.sqrt(errx/nodes)
        erry = np.sqrt(erry/nodes)

        if errx < tol and erry < tol:
            print("Solution converged at %d iteration"%(iter) )
            break

        #print out file, change to True if you want to keep frame data
        if False:
            if iter%steps == 0:
                print(iter, errx, erry)
                num = str(iter)
                filename = "./data/naca-C-grid-pb"+num.zfill(5)+".csv"
                df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
                df.to_csv(filename, index=False)

    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]
            
    iter +=1


print("Output file")
fileout = filename[0:24]+"-para-C-grid.csv" 
df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)

