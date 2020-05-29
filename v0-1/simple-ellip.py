import numpy as np
import pandas as pd

imax = 41
jmax = 41

L = 1.0

#for initial grid
dx = L / (imax-1)
dy = L / (jmax-1)

x = np.ndarray((imax,jmax))
y = np.ndarray((imax, jmax))

z = np.zeros_like(x)
xn = np.zeros_like(x)
yn = np.zeros_like(y)

#Initial grid
for i in range(imax):
    for j in range(jmax):
        x[i][j] = i*dx
        y[i][j] = j*dy


xn[:][:] = x[:][:]
yn[:][:] = y[:][:]
        
eps = 1.0e-5
itermax = 200
steps = 100
iter = 0

P = np.zeros_like(x)
Q = np.zeros_like(x)

p1 = np.ndarray((imax))
p2 = np.ndarray((imax))
q1 = np.ndarray((imax))
q2 = np.ndarray((imax))

while iter < itermax+1:

    for i in range(1, imax-1):
        j = 1 #inner boundary
        #calculating metrics
        x_xi = 0.5*( x[i+1][j] - x[i-1][j] )
        y_xi = 0.5*( y[i+1][j] - y[i-1][j] )

        sn = 0.01 # = np.sqrt( x_eta**2 + y_eta**2)
        x_eta = sn * ( -y_xi ) / np.sqrt( x_xi**2 + y_xi**2 +eps )
        y_eta = sn * ( x_xi ) / np.sqrt( x_xi**2 + y_xi**2 + eps )
        
        x_xieta = 0.25*(x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
        y_xieta = 0.25*(y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )

        x_xi2 = (x[i+1][j] -2.0*x[i][j] + x[i-1][j])
        y_xi2 = ( y[i+1][j] -2.0*y[i][j]  + y[i-1][j] )

        x_eta2 = 0.5*( -7.0*x[i][j-1] +8.0*x[i][j]  - x[i][j+1] ) - 3.0*x_eta
        y_eta2 = 0.5*( -7.0*y[i][j-1] +8.0*y[i][j]  - y[i][j+1] ) - 3.0*y_eta

        a = x_eta**2 + y_eta**2
        b = x_xi*x_eta + y_xi*y_eta
        c = x_xi**2 + y_xi**2
        
        #Jacobian
        J = x_xi*y_eta - x_eta*y_xi
        J2 = J**2
        R11 = -J2*( a*x_xi2 - 2.0*b*x_xieta + c*x_eta2  )
        R12 = -J2*( a*y_xi2 - 2.0*b*y_xieta + c*y_eta2  )
        p1[i] =   J*( y_eta*R11 - x_eta*R12 )
        q1[i] =  J*( -y_xi*R11 + x_xi*R12 )

        #outer boundary
        j = jmax-2
        #calculating metrics
        x_xi = 0.5*( x[i+1][j] - x[i-1][j] )
        y_xi = 0.5*( y[i+1][j] - y[i-1][j] )
        sn = 2.0 # = np.sqrt( x_eta**2 + y_eta**2), must be >= 1.0
        x_eta = sn * ( -y_xi ) / np.sqrt( x_xi**2 + y_xi**2 +eps )
        y_eta = sn * ( x_xi ) / np.sqrt( x_xi**2 + y_xi**2 + eps )
        
        x_xieta = 0.25*(x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
        y_xieta = 0.25*(y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )

        x_xi2 = (x[i+1][j] -2.0*x[i][j] + x[i-1][j])
        y_xi2 = ( y[i+1][j] -2.0*y[i][j]  + y[i-1][j] )

        x_eta2 = 0.5*( -7.0*x[i][jmax-1] +8.0*x[i][j]  - x[i][j-1] ) + 3.0*x_eta
        y_eta2 = 0.5*( -7.0*y[i][jmax-1] +8.0*y[i][j]  - y[i][j-1] ) + 3.0*y_eta

        a = x_eta**2 + y_eta**2
        b = x_xi*x_eta + y_xi*y_eta
        c = x_xi**2 + y_xi**2
        
        #Jacobian
        J = x_xi*y_eta - x_eta*y_xi
        J2 = J**2
        R11 = -J2*( a*x_xi2 - 2.0*b*x_xieta + c*x_eta2  )
        R12 = -J2*( a*y_xi2 - 2.0*b*y_xieta + c*y_eta2  )
        p2[i] =  J*( y_eta*R11 - x_eta*R12 )
        q2[i] =  J*( -y_xi*R11 + x_xi*R12 )
        

    #filling in
    p1[0] = p1[1]
    q1[0] = q1[1]
    p1[imax-1] = p1[imax-2]
    q1[imax-1] = q1[imax-2]

    #filling in
    p2[0] = p2[1]
    q2[0] = q2[1]
    p2[imax-1] = p2[imax-2]
    q2[imax-1] = q2[imax-2]

    #filing in P[i][j], Q[i][j]
    for i in range(imax):
        for j in range(jmax): # j=0 up to jmax-1
            a1 = 0.7
            a2 = 0.7
            P[i][j] = p1[i]*np.exp(-a1*j) + p2[i]*np.exp(-a2*(jmax-1 - j) )
            Q[i][j] = q1[i]*np.exp(-a1*j) + q2[i]*np.exp(-a2*(jmax-1 - j) )
    
    
    # grid
    for i in range(1, imax-1):
        for j in range(1, jmax-1):
            x_xi = 0.5*( x[i+1][j] - x[i-1][j] )
            y_xi = 0.5*( y[i+1][j] - y[i-1][j] )

            x_eta =0.5*( x[i][j+1] - x[i][j-1] )
            y_eta =0.5*( y[i][j+1] - y[i][j-1] )

            x_a = (x[i+1][j]  + x[i-1][j]) 
            x_b = (x[i][j+1]  + x[i][j-1])
            y_a = ( y[i+1][j]  + y[i-1][j] ) 
            y_b = ( y[i][j+1]  + y[i][j-1] )
            
            x_xi2   =  ( x[i+1][j] -2.0*x[i][j] + x[i-1][j])
            x_eta2 =  ( x[i][j+1] -2.0*x[i][j]  + x[i][j-1])
            y_xi2    = ( y[i+1][j] -2.0*y[i][j]  + y[i-1][j] )
            y_eta2  = ( y[i][j+1] -2.0*y[i][j]  + y[i][j-1] )

            x_xieta = 0.25*(x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1] )
            y_xieta = 0.25*(y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1] )

            J = x_xi*y_eta - x_eta*y_xi #Jacobian

            Sx = - (1.0/(J*J +eps ) )*( P[i][j]*x_xi + Q[i][j]*x_eta  )
            Sy = - (1.0/(J*J +eps ) )*( P[i][j]*y_xi + Q[i][j]*y_eta  )
         
            a = x_eta**2 + y_eta**2
            b = x_xi*x_eta + y_xi*y_eta
            c = x_xi**2 + y_xi**2

            xn[i][j] =0.5/( a + c + eps )*( a*x_a + c*x_b - 2.0*b*x_xieta - Sx )
            yn[i][j] =0.5/( a + c + eps )*( a*y_a + c*y_b - 2.0*b*y_xieta - Sy )

    if iter > 0:
        #error calculation
        errx = 0.0
        erry = 0.0
        nodes = 0
        for i in range(1, imax-1):
            for j in range(1, jmax-1):
                errx +=  (xn[i][j] - x[i][j])**2
                erry +=  (yn[i][j] - y[i][j])**2 
                nodes +=1

        errx = np.sqrt( errx / nodes )
        erry = np.sqrt( erry / nodes )

    if iter>0 and iter%steps == 0:
        print(iter, errx, erry)

    x[:][:] = xn[:][:]
    y[:][:] = yn[:][:]
    iter +=1

    
            
df = pd.DataFrame({"x": x.flatten(), "y": y.flatten(), "z": z.flatten()})
filename = "laplace-simple.csv"
df.to_csv(filename, index=False)
