import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd #to write out .csv file; I am getting lazy

x = np.arange(0.0, 1.0, 0.1)
y = np.arange(1.0, 2.0, 0.1)

df = pd.DataFrame( {'x': x, 'y': y} )
#print(df)

#life made easier by pandas.
df.to_csv('simple.csv', index=False)

#Now to the real business
imax = 21
jmax = 21

imax = int(input("Enter imax ") )
jmax = imax

#This is just us.. human
xg = np.ndarray((imax,jmax), dtype=float)
yg = np.ndarray((imax,jmax), dtype=float)
zg = np.ndarray((imax,jmax), dtype=float)

xc = 0.0
yc = 0.0
A = 3.0
B = 3.0
L = 1.0

#a1 = beta + (2.0*alpha+1.0)*(i*dx) - 2.0*alpha
#a2 = beta  - (2.0*alpha+1.0)*(i*dx) + 2.0*alpha
#b1 = beta + 1.0
#b2 = beta - 1.0
#x[i] = alpha + (1.0 - alpha)*( np.log(a1/a2) / np.log(b1/b2 )  ) 


dxi = 1.0 / (imax-1)
deta = 1.0 /(jmax-1)
alpha = 0.5
beta = 2.0
c = 1.0
for i in range(imax):
    for j in range(jmax):
        #xg[i][j] = L*(i*dxi) + A*(xc - L*(i*dxi) )*(1.0 - (i*dxi) )*(i*dxi)
        #yg[i][j] = L*(j*deta) + B*(xc - L*(j*deta) )*(1.0 - (j*deta) ) * (j*deta)

        #a1 = beta + (2.0*alpha+1.0)*(i*dxi) - 2.0*alpha
        #a2 = beta  - (2.0*alpha+1.0)*(i*dxi) + 2.0*alpha
        #b1 = beta + 1.0
        #b2 = beta - 1.0
        #xg[i][j] = alpha + (1.0 - alpha)*( np.log(a1/a2) / np.log(b1/b2 )  )

        #a1 = beta + (2.0*alpha+1.0)*(j*deta) - 2.0*alpha
        #a2 = beta  - (2.0*alpha+1.0)*(j*deta) + 2.0*alpha
        #b1 = beta + 1.0
        #b2 = beta - 1.0
        #yg[i][j] = alpha + (1.0 - alpha)*( np.log(a1/a2) / np.log(b1/b2 )  )

        a1 = ( np.exp(beta*i*dxi) - 1.0 ) /abs( np.exp(beta*c) - 1.0 )
        a2 = -( np.exp(-beta*i*dxi) - 1.0 ) /abs( np.exp(-beta*c) - 1.0 )
        alpha = (i*dxi)/c
        xg[i][j] = (1.0 - alpha)*a1 + alpha*a2

        a1 = ( np.exp(beta*j*deta) - 1.0 ) /abs( np.exp(beta*c) - 1.0 )
        a2 = -( np.exp(-beta*j*deta) - 1.0 ) /abs( np.exp(-beta*c) - 1.0 )
        alpha = (j*deta)/c
        yg[i][j] = (1.0 - alpha)*a1 + alpha*a2
        
        zg[i][j] = 0.0

df = pd.DataFrame( {'x': xg.flatten(), 'y': yg.flatten(), 'z': zg.flatten()} )
df.to_csv('test-grid.csv', index=False)

#plot immediately
#poor quality
X, Y = np.meshgrid( xg.flatten(), yg.flatten() )
plt.plot(X, Y, 'k-', linewidth=0.5)
plt.plot( np.transpose(X), np.transpose(Y) , 'k-', linewidth=0.5)
plt.show()

