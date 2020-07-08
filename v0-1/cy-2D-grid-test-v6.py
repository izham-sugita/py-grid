import numpy as np
import pandas as pd

etamax = 31
rmax = 31
ans ="n"
ans = input("Change the value of rmax, etamax? ")
if ans == "y":
    rmax = int(input("Enter rmax value: ") )
    etamax = int(input("Enter etamax value: ") )
else:
    print("rmax = ", rmax)
    print("etamax = ", etamax)

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

outer = 8.0
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
tol = 1.0e-6
itermax = 8000
iter = 0
steps = 100
omega = 2.0/3.0

#source term
p0 = np.ndarray((rmax))
q0 = np.ndarray((rmax))

p1 = np.ndarray((rmax))
q1 = np.ndarray((rmax))

p0[:]  = 0.0
q0[:]  = 0.0

p1[:]  = 0.0
q1[:]  = 0.0

wp = 0.0001
wq = 0.0001

p0_n = np.zeros_like(p0)
q0_n = np.zeros_like(q0)

p1_n = np.zeros_like(p0)
q1_n = np.zeros_like(q0)


while iter<itermax:

    for i in range(1, rmax-1):
        #calculating boundary area source term
        j=0
        x_xi = 0.5*( xg[i+1][j] - xg[i-1][j] )
        y_xi = 0.5*( yg[i+1][j] - yg[i-1][j] )
        c = x_xi**2 + y_xi**2 

        sn = 0.01 #set small mesh size at solid boundary
        x_eta = sn * ( - y_xi ) / np.sqrt( c )
        y_eta = sn * ( x_xi ) / np.sqrt( c )

        J = 1.0/( x_xi*y_eta - x_eta*y_xi )
        a = x_eta**2 + y_eta**2
        b = x_xi*x_eta + y_xi*y_eta
        c = x_xi**2 + y_xi**2

        x_xi2 = (xg[i+1][j] - 2.0*xg[i][j] + xg[i-1][j] )
        y_xi2 = (yg[i+1][j] - 2.0*yg[i][j] + yg[i-1][j] )
        x_eta2 = 0.5*( -7.0*xg[i][j] + 8.0*xg[i][j+1] - xg[i][j+2]  ) - 3.0*x_eta
        y_eta2 = 0.5*( -7.0*yg[i][j] + 8.0*yg[i][j+1] - yg[i][j+2]  ) - 3.0*y_eta

        x_xieta = 0.25*(xg[i+1][j+1] - xg[i+1][0] - xg[i-1][j+1] + xg[i-1][0] )
        y_xieta = 0.25*(yg[i+1][j+1] - yg[i+1][0] - yg[i-1][j+1] + yg[i-1][0] )

        R1 = -J*J*( a*x_xi2 -2.0*b*x_xieta + c*x_eta2 )
        R2 = -J*J*( a*y_xi2 -2.0*b*y_xieta + c*y_eta2 )

        p0_n[i] = p0[i] + wp* ( J*( y_eta*R1 - x_eta*R2 ) - p0[i] )
        q0_n[i] = q0[i] + wq* ( J*( -y_xi*R1 + x_xi*R2 ) - q0[i] )

        #outer boundary control function
        #j=jmax-1 #real outer boundary
        j=etamax-1
        x_xi = 0.5*( xg[i+1][j] - xg[i-1][j] )
        y_xi = 0.5*( yg[i+1][j] - yg[i-1][j] )
        c = x_xi**2 + y_xi**2 

        sn = 1.0 #set bigger mesh size at outer boundary
        x_eta = sn * (  -y_xi ) / np.sqrt( c ) 
        y_eta = sn * ( x_xi ) / np.sqrt( c )

        J = 1.0/( x_xi*y_eta - x_eta*y_xi )
        a = x_eta**2 + y_eta**2
        b = x_xi*x_eta + y_xi*y_eta
        c = x_xi**2 + y_xi**2

        x_xi2 = (xg[i+1][j] - 2.0*xg[i][j] + xg[i-1][j] )
        y_xi2 = (yg[i+1][j] - 2.0*yg[i][j] + yg[i-1][j] )
        x_eta2 = 0.5*( -7.0*xg[i][j] + 8.0*xg[i][j-1] - xg[i][j-2]  ) + 3.0*x_eta
        y_eta2 = 0.5*( -7.0*yg[i][j] + 8.0*yg[i][j-1] - yg[i][j-2]  ) + 3.0*y_eta

        # j+1 = jmax; out of bound, j+1 == jmax-1 for simplicity
        x_xieta = 0.25*(xg[i+1][etamax-1] - xg[i+1][j-1] - xg[i-1][etamax-1] + xg[i-1][j-1] )
        y_xieta = 0.25*(yg[i+1][etamax-1] - yg[i+1][j-1] - yg[i-1][etamax-1] + yg[i-1][j-1] )

        R1 = -J*J*( a*x_xi2 -2.0*b*x_xieta + c*x_eta2 )
        R2 = -J*J*( a*y_xi2 -2.0*b*y_xieta + c*y_eta2 )

        p1_n[i] = p1[i] + wp* ( J*( y_eta*R1 - x_eta*R2 ) - p1[i] )
        q1_n[i] = q1[i] + wq* ( J*( -y_xi*R1 + x_xi*R2 ) - q1[i] )
        

    p0_n[0] = p0_n[1]
    p0_n[rmax-1] = p0_n[rmax-2]
    q0_n[0] = q0_n[1]
    q0_n[rmax-1] = q0_n[rmax-2]

    p1_n[0] = p1_n[1]
    p1_n[rmax-1] = p1_n[rmax-2]
    q1_n[0] = q1_n[1]
    q1_n[rmax-1] = q1_n[rmax-2]

    for i in range(1, rmax-1):
        #for j in range(1, etamax-1, 1):
        for j in range(etamax-2, 0, -1):

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

            #Sorenson control function
            J = (x_r*y_eta - x_eta*y_r)
            J2 = J*J
            a1 = 0.7
            a2 = 0.7
            pp = p0_n[i]*np.exp( -a1*(j-0) ) + p1_n[i]*np.exp( -a2*(etamax-1 - j) )
            qq = q0_n[i]*np.exp( -a1*(j-0) ) + q1_n[i]*np.exp( -a2*(etamax-1 - j) )
            Sx = -J2*( pp*x_r + qq*x_eta )
            Sy = -J2*( pp*y_r + qq*y_eta )
                        
            xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta -Sx  )
            ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta -Sy )

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

         #Sorenson control function
        J = (x_r*y_eta - x_eta*y_r)
        J2 = J*J
        a1 = 0.7
        a2 = 0.7
        pp = p0_n[i]*np.exp( -a1*(j-0) ) + p1_n[i]*np.exp( -a2*(etamax-1 - j) )
        qq = q0_n[i]*np.exp( -a1*(j-0) ) + q1_n[i]*np.exp( -a2*(etamax-1 - j) )
        Sx = -J2*( pp*x_r + qq*x_eta )
        Sy = -J2*( pp*y_r + qq*y_eta )
         
        xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta -Sx  )
        ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta -Sy  )

        xgn[i][j] = omega*xgn[i][j] + (1.0 - omega)*xg[i][j]
        ygn[i][j] = omega*ygn[i][j] + (1.0 - omega)*yg[i][j]

    #Periodic boundary
    for j in range(1, etamax-1):
        i = rmax-1
        xgn[i][j] = xgn[0][j]
        ygn[i][j] = ygn[0][j]

    #enforced Dirichlet orthogonality
    #relax = 0.01
    #for i in range(1, rmax-1):
    #    for j in range(1, etamax-1):
    #        x_r = 0.5*(xgn[i+1][j] - xgn[i-1][j])
    #        y_r = 0.5*(ygn[i+1][j] - ygn[i-1][j])
    #        x_r2 = xgn[i+1][j] - 2.0*xgn[i][j] + xgn[i-1][j]
    #        y_r2 = ygn[i+1][j] - 2.0*ygn[i][j] + ygn[i-1][j]

    #        g11 = x_r**2 + y_r**2

    #        x_eta_0 = xg[i][j] - xg[i][j-1]
    #        y_eta_0 = yg[i][j] - yg[i][j-1]
            
    #        sn = (-y_r*x_eta_0 + x_r*y_eta_0)

     #       x_eta = sn * ( -y_r ) / np.sqrt( g11 ) # -y_r is original
     #       y_eta = sn * ( x_r ) / np.sqrt( g11 )  # +x_r is original
        
            #x_eta2 = 0.5*( -7.0*xgn[i][j] + 8.0*xgn[i][j+1] - xgn[i][j+2]  ) - 3.0*x_eta
            #y_eta2 = 0.5*( -7.0*ygn[i][j] + 8.0*ygn[i][j+1] - ygn[i][j+2]  ) - 3.0*y_eta

     #       x_eta2 =  xgn[i][j+1] -2.0*xgn[i][j] - xgn[i][j-1]
     #       y_eta2 =  ygn[i][j+1] - 2.0*ygn[i][j] - ygn[i][j-1]
            
     #       g22 = x_eta**2 + y_eta**2

     #       P = -( x_r*x_r2 + y_r*y_r2 ) / (g11 +eps) \
     #       - (x_r*x_eta2 + y_r*y_eta2)/(g22+eps)                

     #       Q = -( x_eta*x_eta2 + y_eta*y_eta2 ) / (g22 +eps) \
     #       - (x_eta*x_r2 + y_eta*y_r2)/(g11+eps)
        
     #       xgn[i][j] =xgn[i][j] + relax*( -g22*P*x_r - g11*Q*x_eta )
     #       ygn[i][j] =ygn[i][j] + relax*( -g22*P*y_r - g11*Q*y_eta )

    #relax =0.01
    #for j in range(1, etamax-1):
    #    i=0
    #    x_r = 0.5*(xgn[1][j] - xgn[rmax-2][j])
    #    y_r = 0.5*(ygn[1][j] - ygn[rmax-2][j])
    #    x_r2 = xgn[1][j] - 2.0*xgn[i][j] + xgn[rmax-2][j]
    #    y_r2 = ygn[1][j] - 2.0*ygn[i][j] + ygn[rmax-2][j]
    #    g11 = x_r**2 + y_r**2
    #    x_eta_0 = xg[i][j] - xg[i][j-1]
    #    y_eta_0 = yg[i][j] - yg[i][j-1]
    #    sn = (-y_r*x_eta_0 + x_r*y_eta_0)
    #    x_eta = sn * ( -y_r ) / np.sqrt( g11 ) # -y_r is original
    #    y_eta = sn * ( x_r ) / np.sqrt( g11 )  # +x_r is original

        #x_eta2 = 0.5*( -7.0*xgn[i][j] + 8.0*xgn[i][j+1] - xgn[i][j+2]  ) - 3.0*x_eta
        #y_eta2 = 0.5*( -7.0*ygn[i][j] + 8.0*ygn[i][j+1] - ygn[i][j+2]  ) - 3.0*y_eta

    #    x_eta2 =  xgn[i][j+1] -2.0*xgn[i][j] - xgn[i][j-1]
    #    y_eta2 =  ygn[i][j+1] -2.0*ygn[i][j] - ygn[i][j-1]

    #    g22 = x_eta**2 + y_eta**2
    #    P = -( x_r*x_r2 + y_r*y_r2 ) / (g11 +eps) \
    #    - (x_r*x_eta2 + y_r*y_eta2)/(g22+eps)                
    #    Q = -( x_eta*x_eta2 + y_eta*y_eta2 ) / (g22 +eps) \
    #    - (x_eta*x_r2 + y_eta*y_r2)/(g11+eps)
    #    xgn[i][j] =xgn[i][j] + relax*( -g22*P*x_r - g11*Q*x_eta )
    #    ygn[i][j] =ygn[i][j] + relax*( -g22*P*y_r - g11*Q*y_eta )

    #Periodic boundary
    #for j in range(1, etamax-1):
    #    i = rmax-1
    #    xgn[i][j] = xgn[0][j]
    #    ygn[i][j] = ygn[0][j]


        

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
        #    num = str(iter)
        #    filestep = "./data/cy-2D-ellip-"+num.zfill(5)+".csv"
        #    df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
        #    df.to_csv(filestep, index=False)
            
            
    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]

    #update
    p0[:] = p0_n[:]
    q0[:] = q0_n[:]

    p1[:] = p1_n[:]
    q1[:] = q1_n[:]
            
    iter +=1


#pandas!!!
fileout = filename[0:12]+"-ellip-grid.csv" 
df = pd.DataFrame({'x': xg.flatten(), 'y': yg.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)

#coordinates for cell center
#using list data structure for ease of accounting
cell_x = []
cell_y = []
cell_z = []

for i in range(0, rmax-1):
    for j in range(0, etamax-1):
        cell_center_x = 0.25*(xg[i][j] + xg[i+1][j] + xg[i+1][j+1] + xg[i][j+1] )
        cell_center_y = 0.25*(yg[i][j] + yg[i+1][j] + yg[i+1][j+1] + yg[i][j+1] )
        cell_x.append(cell_center_x)
        cell_y.append(cell_center_y)
        cell_z.append(0.0)

fileout = "cell_center.csv"
dict = {'x':cell_x, 'y':cell_y, 'z':cell_z }
df = pd.DataFrame(dict)
df.to_csv(fileout, index=False)
