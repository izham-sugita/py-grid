import numpy as np
import pandas as pd

etamax = 51
rmax = 128

#read NACA0012 data file
#dat = pd.read_csv("naca0012.txt", delimiter="\t")
#dat = pd.read_csv("naca0012-201.csv")
filename = "naca6409-101.csv"
dat = pd.read_csv(filename)
naca = pd.DataFrame(dat)
airfoil = naca.to_numpy()
rmax = int(len(airfoil))

print( rmax, etamax )

xg = np.ndarray((rmax, etamax), dtype = np.float32)
xgn = np.ndarray((rmax, etamax), dtype = np.float32)
yg = np.ndarray((rmax, etamax), dtype = np.float32)
ygn = np.ndarray((rmax, etamax), dtype = np.float32)

#for paraview purpose
zg = np.zeros_like(yg)

outer = 5.0
inner = 1.0


dr = (outer - inner)/(etamax-1)
dseta = 2.0*(np.pi)/(rmax-1)

xg[:][:] = 0.0
yg[:][:] = 0.0

# xcenter moved to 0.5 because airfoil is 0.0<= c <= 1.0
xcenter = 0.5
ycenter = 0.0

#O-type boundary input is the simplest
#This is the part to change for geometry
#Boundary condition before eliptic equation solution
#Adjusted position
angle = 5.0 #degree
angle = (angle*np.pi)/180.0 #convert to rad
x_r = 0.5
y_r = 0.0
for i in range(rmax):
    j = 0
    # data on airfoil surface
    #xg[i][j] =  (inner + j*dr)*np.cos( -i*dseta) + xcenter
    #yg[i][j] =  (inner + j*dr)*np.sin( -i*dseta)

    #Original position
    xg[i][j] =  airfoil[i][0]
    yg[i][j] =  airfoil[i][1]

    #Adjusted position
    xg[i][j] = (xg[i][j] - x_r)*np.cos(angle) - (yg[i][j] - y_r)*np.sin(angle) + x_r  
    yg[i][j] = (yg[i][j] - y_r)*np.cos(angle) -(xg[i][j] - x_r)*np.sin(angle) + y_r
    
    j = 1
    xg[i][j] =  (inner + j*dr)*np.cos( -i*dseta) + xcenter
    yg[i][j] =  (inner + j*dr)*np.sin( -i*dseta) + ycenter
    
    j = etamax-1
    xg[i][j] = (inner + j*dr)*np.cos( -i*dseta) + xcenter
    yg[i][j] = (inner + j*dr)*np.sin( -i*dseta) + ycenter

print(xg[0][0], xg[rmax-1][0])
    
for j in range(1, etamax):

    i = 0
    xg[i][j] = (inner + j*dr)*np.cos( -i*dseta) + xcenter
    yg[i][j] = (inner + j*dr)*np.sin( -i*dseta) + ycenter

    i = rmax-1
    xg[i][j] = (inner + j*dr)*np.cos( -i*dseta) + xcenter
    yg[i][j] = (inner + j*dr)*np.sin( -i*dseta) + ycenter
    

#for i in range(1, rmax-1):
#    for j in range(2, etamax-1):
#        xg[i][j] = xg[i][0]
#        yg[i][j] = yg[i][0]

print(xg[0][0], xg[rmax-1][0])

#not a very good initial distribution
xgn[:][:] = xg[:][:]
ygn[:][:] = yg[:][:]

print(len(xg))
print(len(yg))

#pandas!!!
df = pd.DataFrame({'x': xg.flatten(), 'y': yg.flatten(), 'z': zg.flatten() } )
df.to_csv("naca-pb-init.csv", index=False)

#Solving eliptic equation by iteration method
# remove dseta, dr in the elliptic equation solver
eps = 1.0e-6
tol = 1.0e-6
itermax = 5000
iter = 0
steps = 100


A = 0.1*dr
B = 0.1*dr

while iter<itermax:

    for i in range(0, rmax-1):
        for j in range(0, etamax-2):

            if i == 0:
                x_eta2 = (xg[i+1][j]  - 2.0*xg[i][j]  + xg[rmax -1][j])
                y_eta2 = ( yg[i+1][j] - 2.0*yg[i][j]  + yg[rmax -1][j] )
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
            # rad: based curvature lenght from surface to outer boundary
            r1 = np.sqrt( (xg[i][etamax-1] - xg[i][0] )**2 + (yg[i][etamax-1] - yg[i][0])**2  )
            r0 = np.sqrt( (xg[i][0] - xcenter )**2 + (yg[i][0] - ycenter )**2  )
            rvar = (j+1)*dr + r0
            ratio = np.log( (rvar+eps) / r0 ) / np.log( (r1+eps)/r0 )
            
            #Parabolic solver
            #Sx = ratio*(xg[i][etamax-1] - xg[i][j] ) / float( (etamax-1) - j )
            Sx = (xg[i][etamax-1] - xg[i][j] ) / float( (etamax-1) - j )
            xgn[i][j+1] = xg[i][j] + A*(x_eta2) - 2.0*B*x_reta + Sx
            
            #Sy = ratio*(yg[i][etamax-1] - yg[i][j] ) / float( (etamax-1) - j )
            Sy = (yg[i][etamax-1] - yg[i][j] ) / float( (etamax-1) - j )
            ygn[i][j+1] = yg[i][j] + A*(y_eta2) - 2.0*B*y_reta + Sy
            

    #upgrade BC, periodic
    for j in range(etamax):
        xgn[rmax-1][j] = xgn[0][j]
        ygn[rmax-1][j] = ygn[0][j]


    if iter> 10:
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
            print("Solution converged at %d iteration"%(iter) )
            break

        #print out file
        #pandas!!!
        if iter%steps == 0:
            print(iter, errx, erry)
            num = str(iter)
            filename = "./data/naca-pb"+num.zfill(5)+".csv"
            df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
            df.to_csv(filename, index=False)

    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]
            
    iter +=1

    

#pandas!!!
fileout = filename[0:12]+"-grid.csv" 
df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)
