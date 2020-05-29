import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#read NACA0012 data file
#filename = "naca0012-101-101-030-070-input.csv"
#filename = "naca0012-501-125-150-350-input.csv"
#filename="naca0012-207-051-083-123-input.csv"
filename = "naca0012-101-025-030-070-input.csv"
#filename = "naca0012-101-301-030-070-input.csv"
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
print("Total nodes in xi-axis")
print(all_nodes,"\n")

xg = np.ndarray((imax*jmax))
yg = np.ndarray((imax*jmax))

for i in range(imax*jmax):
    xg[i] = airfoil[i][0]
    yg[i] = airfoil[i][1]

#print(xg, yg)

xg = xg.reshape((imax, jmax))
yg = yg.reshape((imax, jmax))

#for i in range(imax):
#    j = 0
#    print(xg[i][j], yg[i][j])

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

#parabolic iteration has a minimum: the number of mesh in
#the propagating direction.
itermax = 2*jmax + 1

iter = 0
steps = 100

A = 0.005
B = 0.005

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
            #r1 = np.sqrt( (xg[i][jmax-1] - xg[i][0] )**2 + (yg[i][jmax-1] - yg[i][0])**2  )
            #r0 = np.sqrt( (xg[i][0] - xcenter )**2 + (yg[i][0] - ycenter )**2  )
            #rvar = (j+1)*(r1 - r0)/(jmax-1) + r0
            #ratio = np.log( (rvar+eps) / r0 ) / np.log( (r1+eps)/r0 )
            #Sx = ratio*(xg[i][jmax-1] - xg[i][j] ) / float( (jmax-1) - j )
            #Sy = ratio*(yg[i][jmax-1] - yg[i][j] ) / float( (jmax-1) - j )
            
            #Parabolic solver
            Sx = (xg[i][jmax-1] - xg[i][j] ) / float( (jmax-1) - j )
            xgn[i][j+1] = xg[i][j] + A*(x_eta2) - 2.0*B*x_reta + Sx

            Sy = (yg[i][jmax-1] - yg[i][j] ) / float( (jmax-1) - j )
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

        #print out file
        #pandas!!!
        if iter%steps == 0:
            print(iter, errx, erry)
            #num = str(iter)
            #filename = "./data/naca-C-grid-pb"+num.zfill(5)+".csv"
            #df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
            #df.to_csv(filename, index=False)

    #update
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]
            
    iter +=1

print("Output data from parabolic equation solver")
fileout = filename[0:24]+"-C-grid-para.csv"
df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)

print("Initiate grid smoothing")
print("Solving elliptic equation with source term P0, Q0")
#using elliptic equation as smoothing
iter = 0
steps = 100
itermax = 3001
omega = 0.5
tol = 1.0e-4

#source term
p0 = np.ndarray((imax))
q0 = np.ndarray((imax))

p1 = np.ndarray((imax))
q1 = np.ndarray((imax))

p0[:]  = 0.0
q0[:]  = 0.0

p1[:]  = 0.0
q1[:]  = 0.0

wp = 0.00005
wq = 0.00005

p0_n = np.zeros_like(p0)
q0_n = np.zeros_like(q0)

p1_n = np.zeros_like(p0)
q1_n = np.zeros_like(q0)

while iter < itermax:

    for i in range(1, imax-1):
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
        j=jmax-1
        x_xi = 0.5*( xg[i+1][j] - xg[i-1][j] )
        y_xi = 0.5*( yg[i+1][j] - yg[i-1][j] )
        c = x_xi**2 + y_xi**2 

        sn = 5.0 #set small mesh size at solid boundary
        x_eta = sn * ( - y_xi ) / np.sqrt( c )
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
        x_xieta = 0.25*(xg[i+1][jmax-1] - xg[i+1][j-1] - xg[i-1][jmax-1] + xg[i-1][j-1] )
        y_xieta = 0.25*(yg[i+1][jmax-1] - yg[i+1][j-1] - yg[i-1][jmax-1] + yg[i-1][j-1] )

        R1 = -J*J*( a*x_xi2 -2.0*b*x_xieta + c*x_eta2 )
        R2 = -J*J*( a*y_xi2 -2.0*b*y_xieta + c*y_eta2 )

        p1_n[i] = p1[i] + wp* ( J*( y_eta*R1 - x_eta*R2 ) - p1[i] )
        q1_n[i] = q1[i] + wq* ( J*( -y_xi*R1 + x_xi*R2 ) - q1[i] )
        

    p0_n[0] = p0_n[1]
    p0_n[imax-1] = p0_n[imax-2]
    q0_n[0] = q0_n[1]
    q0_n[imax-1] = q0_n[imax-2]

    p1_n[0] = p1_n[1]
    p1_n[imax-1] = p1_n[imax-2]
    q1_n[0] = q1_n[1]
    q1_n[imax-1] = q1_n[imax-2]
    
    
    #elliptic main loop
    for i in range(1, imax-1):
        for j in range(1, jmax-1):
            
            x_eta = 0.5*(xg[i][j+1]-xg[i][j-1])
            x_r = 0.5*(xg[i+1][j]-xg[i-1][j])
            y_eta = 0.5*(yg[i][j+1] - yg[i][j-1])
            y_r = 0.5*(yg[i+1][j] - yg[i-1][j] )
            
            x_r2 = (xg[i+1][j]  + xg[i-1][j]) 
            x_eta2 = (xg[i][j+1]  + xg[i][j-1])
            y_r2 = ( yg[i+1][j]  + yg[i-1][j] ) 
            y_eta2 = ( yg[i][j+1]  + yg[i][j-1] )
            
            #d2x_dxi2 = (xg[i+1][j] -2.0*xg[i][j] + xg[i-1][j])
            #d2x_deta2 = (xg[i][j+1] -2.0*xg[i][j]  + xg[i][j-1])
            #d2y_dxi2 = ( yg[i+1][j] -2.0*yg[i][j]  + yg[i-1][j] )
            #d2y_deta2 = ( yg[i][j+1] -2.0*yg[i][j]  + yg[i][j-1] )

            x_reta = 0.25*(xg[i+1][j+1] - xg[i+1][j-1] - xg[i-1][j+1] + xg[i-1][j-1] )
            y_reta = 0.25*(yg[i+1][j+1] - yg[i+1][j-1] - yg[i-1][j+1] + yg[i-1][j-1] )

            J = ( x_r*y_eta - x_eta*y_r ) #Jacobian
            J2 = J*J
            
            a = x_eta**2 + y_eta**2
            b = x_eta*x_r + y_eta*y_r
            c = x_r**2 + y_r**2

            #calculate source term
            a1 = 0.7
            a2 = 0.7
            
            #Sx = -J2*( p0_n[i]*np.exp(-a1*(j-0) )*x_r + q0_n[i]*np.exp(-a2*(j-0) )*x_eta )
            #Sy = -J2*( p0_n[i]*np.exp(-a1*(j-0) )*y_r + q0_n[i]*np.exp(-a2*(j-0) )*y_eta )

            pp = p0_n[i]*np.exp( -a1*(j-0) ) + p1_n[i]*np.exp( -a1*(jmax-2 - j) )
            qq = q0_n[i]*np.exp( -a2*(j-0) ) + q1_n[i]*np.exp( -a2*(jmax-2 - j) )
            
            Sx = -J2*( pp*x_r + qq*x_eta )
            Sy = -J2*( pp*y_r + qq*y_eta )
            
            xgn[i][j] = ( 0.5/(a + c + eps ) )*( a*x_r2 + c*x_eta2 - 2.0*b*x_reta - Sx )
            ygn[i][j] = ( 0.5/(a + c + eps ) )*( a*y_r2 + c*y_eta2 - 2.0*b*y_reta - Sy )

            #add relaxation step
            xgn[i][j] = omega*xgn[i][j] + (1.0 - omega)*xg[i][j]
            ygn[i][j] = omega*ygn[i][j] + (1.0 - omega)*yg[i][j]
            

    #upgrade BC, periodic
    #for j in range(etamax):
    #    xgn[rmax-1][j] = xgn[0][j]
    #    ygn[rmax-1][j] = ygn[0][j]
        

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
        #print(iter)

        if errx < tol and erry < tol:
            print("Solution converged at %d iteration"%(iter) )
            break

    if iter%steps == 0:
        print(iter, errx, erry)
        

    #update all
    xg[:][:] = xgn[:][:]
    yg[:][:] = ygn[:][:]

    #update outer boundary so that its moves with orthogonality
    #for i in range(imax):
    #    xg[i][jmax-1] = xgn[i][jmax-2]
    #    yg[i][jmax-1] = ygn[i][jmax-2]
    
    #update p0, q0
    p0[:] = p0_n[:]
    q0[:] = q0_n[:]

    #update p1, q1
    p1[:] = p1_n[:]
    q1[:] = q1_n[:]
    
            
    iter +=1


print("Output file")
#pandas!!!
fileout = filename[0:24]+"-C-grid-para-ellip.csv" 
df = pd.DataFrame({'x': xgn.flatten(), 'y': ygn.flatten(), 'z': zg.flatten() } )
df.to_csv(fileout, index=False)

df = pd.DataFrame( {"p0_n": p0_n.flatten(), "q0_n": q0_n.flatten() } )
df.to_csv("source-debug.csv", index=False)
