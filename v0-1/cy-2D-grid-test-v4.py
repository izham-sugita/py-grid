import numpy as np
import pandas as pd

etamax = 51
rmax = 51

#read NACA0012 data file
#dat = pd.read_csv("naca0012.txt", delimiter="\t")
#dat = pd.read_csv("naca0012-201.csv")
#filename = "naca6409-101.csv"
#dat = pd.read_csv(filename)
#naca = pd.DataFrame(dat)
#airfoil = naca.to_numpy()
#rmax = int(len(airfoil))

filename = "cy-2D-test-d-"
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
xcenter = 0.0
ycenter = 0.0
ans = 'n'
ans = input("Change cylinder center? y/n ")
if ans == 'y':
    xcenter = float(input("x-center? "))
    ycenter = float(input("y-center? "))
else:
    print("Continue default center, 0.0, 0.0 ")

#O-type boundary input is the simplest
#This is the part to change for geometry
#Boundary condition before eliptic equation solution
for i in range(rmax):
    j = 0
    # data on airfoil surface
    # xg[i][j] =  airfoil[i][0]
    # yg[i][j] =  airfoil[i][1]
    xg[i][j] =  (inner + j*dr)*np.cos( -i*dseta) + xcenter
    yg[i][j] =  (inner + j*dr)*np.sin( -i*dseta) + ycenter
    
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
    

#Initiate for elliptic solver only
for i in range(1, rmax-1):
    for j in range(2, etamax-1):
        xg[i][j] = xg[i][0]
        yg[i][j] = yg[i][0]

print(xg[0][0], xg[rmax-1][0])

#not a very good initial distribution
xgn[:][:] = xg[:][:]
ygn[:][:] = yg[:][:]

print(len(xg))
print(len(yg))

#pandas!!!
df = pd.DataFrame({'x': xg.flatten(), 'y': yg.flatten(), 'z': zg.flatten() } )
df.to_csv("cy-2D-ellip-init.csv", index=False)

#Solving eliptic equation by iteration method
# remove dseta, dr in the elliptic equation solver
eps = 1.0e-7 #get it as small as possible; just to avoid zero division
tol = 1.0e-4
itermax = 5000
iter = 0
steps = 100
omega = 2.0/3.0

#control term P, Q
#P = np.ndarray((rmax,etamax))
#Q = np.ndarray((rmax,etamax))

#Initial value
#P[:][:] = 0.0
#Q[:][:] = 0.0

while iter<itermax:

    for i in range(1, rmax-1):
        for j in range(1, etamax-1):

            x_eta = 0.5*(xg[i][j+1]-xg[i][j-1])
            x_r = 0.5*(xg[i+1][j]-xg[i-1][j])
            y_eta = 0.5*(yg[i][j+1] - yg[i][j-1])
            y_r = 0.5*(yg[i+1][j] - yg[i-1][j] )
            x_r2 = (xg[i+1][j]  + xg[i-1][j]) 
            x_eta2 = (xg[i][j+1]  + xg[i][j-1])
            y_r2 = ( yg[i+1][j]  + yg[i-1][j] )
            y_eta2 = ( yg[i][j+1]  + yg[i][j-1] )
            x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][j-1] - xg[i-1][j+1] + xg[i-1][j-1] )
            y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][j-1] - yg[i-1][j+1] + yg[i-1][j-1] )
            
            a = x_eta**2 + y_eta**2
            b = x_eta*x_r + y_eta*y_r
            c = x_r**2 + y_r**2
                        
            xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta  )
            ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta  )

            #add relaxation step
            xgn[i][j] = omega*xgn[i][j] + (1.0 - omega)*xg[i][j]
            ygn[i][j] = omega*ygn[i][j] + (1.0 - omega)*yg[i][j]

            
    for j in range(1, etamax-1):
        i = 0 #imax-1 is equal to i=0 index
        x_eta = 0.5*(xg[i][j+1]-xg[i][j-1])
        x_r = 0.5*(xg[1][j]-xg[rmax-2][j])
        y_eta = 0.5*(yg[i][j+1] - yg[i][j-1])
        y_r = 0.5*(yg[1][j] - yg[rmax-2][j] )
        
        x_r2 = (xg[1][j]  + xg[rmax-2][j]) 
        x_eta2 = (xg[i][j+1]  + xg[i][j-1])
        y_r2 = ( yg[1][j]  + yg[rmax-2][j] )
        y_eta2 = ( yg[i][j+1]  + yg[i][j-1] )

        x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][j-1] - xg[rmax-2][j+1] + xg[rmax-2][j-1] )
        y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][j-1] - yg[rmax-2][j+1] + yg[rmax-2][j-1] )

        a = x_eta**2 + y_eta**2
        b = x_eta*x_r + y_eta*y_r
        c = x_r**2 + y_r**2
         
        xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta  )
        ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta  )

        xgn[i][j] = omega*xgn[i][j] + (1.0 - omega)*xg[i][j]
        ygn[i][j] = omega*ygn[i][j] + (1.0 - omega)*yg[i][j]

    #Periodic boundary
    for j in range(1, etamax-1):
        i = rmax-1
        xgn[i][j] = xgn[0][j]
        ygn[i][j] = ygn[0][j]


    relax = 0.1
    for i in range(1, rmax-1):
        for j in range(1, etamax-3):
            x_r = 0.5*(xgn[i+1][j] - xgn[i-1][j])
            y_r = 0.5*(ygn[i+1][j] - ygn[i-1][j])
            x_r2 = xgn[i+1][j] - 2.0*xgn[i][j] + xgn[i-1][j]
            y_r2 = ygn[i+1][j] - 2.0*ygn[i][j] + ygn[i-1][j]

            g11 = x_r**2 + y_r**2
            sn = 0.01*np.exp(-(j-1)/(etamax-1))

            x_eta = sn * ( -y_r ) / np.sqrt( g11 )
            y_eta = sn * ( x_r ) / np.sqrt( g11 )
        
            x_eta2 = 0.5*( -7.0*xgn[i][j] + 8.0*xgn[i][j+1] - xgn[i][j+2]  ) - 3.0*x_eta
            y_eta2 = 0.5*( -7.0*ygn[i][j] + 8.0*ygn[i][j+1] - ygn[i][j+2]  ) - 3.0*y_eta

            g22 = x_eta**2 + y_eta**2
        

            P = -( x_r*x_r2 + y_r*y_r2 ) / (g11 +eps) \
            - (x_r*x_eta2 + y_r*y_eta2)/(g22+eps)                

            Q = -( x_eta*x_eta2 + y_eta*y_eta2 ) / (g22 +eps) \
            - (x_eta*x_r2 + y_eta*y_r2)/(g11+eps)
        
            xgn[i][j] =xgn[i][j] + relax*( -g22*P*x_r - g11*Q*x_eta )
            ygn[i][j] =ygn[i][j] + relax*( -g22*P*y_r - g11*Q*y_eta )
        
    
    #for control function P,Q
    #P, Q = 0; no control
    #P,Q: function from 6.10, computation blow-up.
    #relax = 0.0001
    #for i in range(1, rmax-1):
    #    for j in range(1, etamax-1):
    #        x_r = 0.5*(xgn[i+1][j] - xgn[i-1][j])
    #        x_eta = 0.5*(xgn[i][j+1] - xgn[i][j-1])
    #        y_r = 0.5*(ygn[i+1][j] - ygn[i-1][j])
    #        y_eta = 0.5*(ygn[i][j+1] - ygn[i][j-1])
    #        x_r2 = xgn[i+1][j] - 2.0*xgn[i][j] + xgn[i-1][j]
    #        y_r2 = ygn[i+1][j] - 2.0*ygn[i][j] + ygn[i-1][j]
    #        x_eta2 = xgn[i][j+1] - 2.0*xgn[i][j] + xgn[i][j-1]
    #        y_eta2 = ygn[i][j+1] - 2.0*ygn[i][j] + ygn[i][j-1]

     #       g11 = x_r**2 + y_r**2
     #       g22 = x_eta**2 + y_eta**2
            
     #       P = -( x_r*x_r2 + y_r*y_r2 ) / (g11 +eps) \
     #       - (x_r*x_eta2 + y_r*y_eta2)/(g22+eps)                
     #       Q = -( x_eta*x_eta2 + y_eta*y_eta2 ) / (g22 +eps) \
     #       - (x_eta*x_r2 + y_eta*y_r2)/(g11+eps)

     #       xgn[i][j] =xgn[i][j] + relax*( -g22*P*x_r - g11*Q*x_eta )
     #       ygn[i][j] =ygn[i][j] + relax*( -g22*P*y_r - g11*Q*y_eta )
        

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
        if iter%steps == 0:
            print(iter, errx, erry)
            num = str(iter)
            filestep = "./data/cy-2D-ellip-"+num.zfill(5)+".csv"
            df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
            df.to_csv(filestep, index=False)

    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]
            
    iter +=1


#pandas!!!
fileout = filename[0:12]+"-ellip-grid.csv" 
df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)
