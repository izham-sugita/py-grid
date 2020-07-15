import numpy as np
import pandas as pd

#importing data file
#dat = pd.read_csv("naca0012.txt", delimiter="\t")
#naca = pd.DataFrame(dat)

#filename = "naca6409-041.csv"
filename = "naca0012-041.csv"
dat = pd.read_csv(filename)
naca = pd.DataFrame(dat)

airfoil = naca.to_numpy()
etamax = int(len(airfoil))

rmax = etamax
print(rmax, etamax)

#etamax = 51

xg = np.ndarray((rmax,etamax), dtype = np.float32)
xgn = np.ndarray((rmax,etamax), dtype = np.float32)

yg = np.ndarray((rmax,etamax), dtype = np.float32)
ygn = np.ndarray((rmax,etamax), dtype = np.float32)

outer = 5.0
inner = 1.0

dr = (outer - inner)/(rmax-1)
dseta = 2.0*(np.pi)/(etamax-1)

xg[:][:] = 0.0
yg[:][:] = 0.0

#O-type boundary input is the simplest
#This is the part to change for geometry
#Boundary condition before eliptic equation solution
for i in range(rmax):
    j = 0
    xg[i][j] =  (inner + i*dr)*np.cos(j*dseta) + 0.5
    yg[i][j] =  (inner + i*dr)*np.sin(j*dseta)

    j = etamax-1
    xg[i][j] = xg[i][0]  #(inner + i*dr)*np.cos(j*dseta)
    yg[i][j] = yg[i][0]  #(inner + i*dr)*np.sin(j*dseta)


surf_y = np.ndarray((etamax))
temp = np.ndarray((etamax))

for j in range(etamax):
    temp[j] = airfoil[j][1]

surf_y = np.flip(temp)


for j in range(etamax):
    i = 0
    #grid data on the airfoil surface
    xg[i][j] = airfoil[j][0]
    yg[i][j] = surf_y[j]  #airfoil[j][1]

    i=1
    xg[i][j] = (0.5 + 0.3)*np.cos(j*dseta) + 0.5
    yg[i][j] = (0.5 + 0.3)*np.sin(j*dseta)
    
    i = rmax-1
    xg[i][j] = (inner + i*dr)*np.cos(j*dseta) + 0.5
    yg[i][j] = (inner + i*dr)*np.sin(j*dseta)

for i in range(2, rmax-1):
    for j in range(etamax):
        xg[i][j] = (0.8 + i*dr)*np.cos(j*dseta) + 0.5
        yg[i][j] = (0.8 + i*dr)*np.sin(j*dseta)
        
xgn[:][:] = xg[:][:]
ygn[:][:] = yg[:][:]

zg = np.zeros_like(xg)

#init data
filename ="naca-init.csv"
df = pd.DataFrame( {"x":xg.flatten(), "y":yg.flatten(), "z":zg.flatten()} )
df.to_csv(filename, index=False)

#Solving eliptic equation by iteration method
# remove dseta, dr in the elliptic equation solver
dseta = 1.0
dr = 1.0
eps = 1.0e-6
tol = 1.0e-6
itermax = 20000
iter = 0
steps = 100

while iter<itermax:

    for i in range(1, rmax-1): #start from 2 to avoid overlap
        for j in range(1, etamax-1):
            x_eta = 0.5*(xg[i][j+1]-xg[i][j-1])
            x_r = 0.5*(xg[i+1][j]-xg[i-1][j])
            y_eta = 0.5*(yg[i][j+1] - yg[i][j-1])
            y_r = 0.5*(yg[i+1][j] - yg[i-1][j] )

            a = x_eta**2 + y_eta**2
            b = x_eta*x_r + y_eta*y_r
            c = x_r**2 + y_r**2

            x_r2 = (xg[i+1][j]  + xg[i-1][j]) 
            x_eta2 = (xg[i][j+1]  + xg[i][j-1])

            y_r2 = ( yg[i+1][j]  + yg[i-1][j] ) 
            y_eta2 = ( yg[i][j+1]  + yg[i][j-1] )

            x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][j-1] - xg[i-1][j+1] + xg[i-1][j-1] )
            y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][j-1] - yg[i-1][j+1] + yg[i-1][j-1] )
                        
            xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta  )
            ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta  )            


    #reference residuals
    if iter > 10:
        errx = 0.0
        erry = 0.0
        nodes = 0
        for i in range(1, rmax-1):
            for j in range(1, etamax-1):
                errx += (xg[i][j] - xgn[i][j])**2
                erry += (yg[i][j] - ygn[i][j])**2
                nodes +=1

        errx = np.sqrt(errx/nodes)
        erry = np.sqrt(erry/nodes)

        if errx < tol and erry < tol:
            print("Solution converged at %d" %(iter) )
            break

        if iter%steps == 0:
            #print out file
            num = str(iter)
            filename = "./data/naca"+num.zfill(5)+".csv"
            df = pd.DataFrame( {"x":xgn.flatten(), "y": ygn.flatten(), "z":zg.flatten()} )
            df.to_csv(filename, index=False)
            print(iter, errx, erry)        

    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]
    #for i in range(2, rmax-1):
    #    for j in range(1, etamax-1):
    #        xg[i][j] = xgn[i][j]
    #        yg[i][j] = ygn[i][j]
            
    iter +=1

    
#output file
#for paraview x = rmax, y = etamax-1
filename = "naca-elip-final.csv"
df = pd.DataFrame( {"x": xg.flatten(), "y":yg.flatten(), "z": zg.flatten()} )
df.to_csv(filename, index=False)
