import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

imax = 51
imax = int( input("Enter total nodes on the top of airfoil ") )

#Input airfoil series XXX
series_number = input("Enter NACA XXXX 4-digit numbers ")

m1 = float(series_number[0])
p2 = float(series_number[1])

t34 = float(series_number[2:4])

yt = np.ndarray((imax))
x = np.ndarray((imax))

#maximum chamber
# 100*m = X
# NACA X*** ( 4-digits)
m = 0.00
#update serial number
m1 = 0.01*m1
m = m + m1

#location of maximum chamber
# 10*p = X
# NACA *X*** (4-digits)
p = 0.0
p2 = 0.1*p2
p = p + p2

#chord length
c = 1.0

#maximum thickness as a fraction of chord i.e. 12 is 12/100
# 100*t = XX
# NACA **XX (4-digits)
t = 0.00
t34 = 0.01*t34
t = t + t34

#x-coordinates
#dx = c / (imax)
#x = np.arange( 0.0, c+dx, dx)
#print(len(x))
#print(x)

dx = c / (imax -1)
#x-values can be stretched; focus distribution on either end
A = -1.0 # A >= 0.0 will get the uniform distribution
xc = 1.0 # xc: the concentration point

for i in range(imax):
    x[i] = c*(i*dx) + A*(xc - c*(i*dx) )*(1.0 - (i*dx) )*(i*dx)


#surface curve
a0 = 0.2969
a1 = -0.1260
a2 = -0.3516
a3 = 0.2843
a4 = -0.1015

yt[:] = 5.0*t*c*( a0*np.sqrt(x[:]/c) + a1*(x[:]/c) + a2*(x[:]/c)**2 + a3*(x[:]/c)**3 \
                  + a4*(x[:]/c)**4 )

#for centerline
eps = 1.0e-6
yc = np.ndarray((imax))
ycdash = np.ndarray((imax))
for i in range(imax):
    if x[i] < (p*c):
        pe = p + eps
        yc[i] = m*( x[i]/(pe*pe) )*( 2.0*p - (x[i]/c) )
        ycdash[i] = (2.0*m)/(pe*pe) *(p - (x[i]/c) )
    elif x[i] >= (p*c):
        yc[i] = m*( ( c - x[i] ) / (1.0 - p)**2 )*( 1.0 +  (x[i]/c) - 2.0*p )
        ycdash[i] = (2.0*m)/( (1.0 - p)**2 ) * ( p - (x[i]/c) )

#defining upper and lower coordinate
xup = np.ndarray((imax))
xlow = np.ndarray((imax))
yup = np.ndarray((imax))
ylow = np.ndarray((imax))

for i in range(imax):
    seta = np.arctan(ycdash[i])
    xup[i] = x[i] - yt[i]*np.sin(seta)
    xlow[i] = x[i] + yt[i]*np.sin(seta)

    yup[i] = yc[i] + yt[i]*np.cos(seta)
    ylow[i] = yc[i] - yt[i]*np.cos(seta)

print(xup)
print()
print(np.flip(xlow) )

tot_nodes = 2*(imax-1) + 1

xnaca = np.ndarray( (tot_nodes) ) #twice the length
ynaca = np.ndarray( (tot_nodes) ) #twice the length

print("Total nodes on the surface: %d"%(tot_nodes) )

#for dummy check
xnaca[:] = 5.0 #out of bound initial value
ynaca[:] = 5.0 #out of bound initial value

tempx = np.ndarray((imax))
tempy = np.ndarray((imax))
tempx = np.flip(xlow)
tempy = np.flip(ylow)

#lower coordinates
for i in range(imax-1):
    xnaca[i] = tempx[i]
    ynaca[i] = tempy[i]

#upper coordinates
for i in range(imax):
    xnaca[imax-1 + i] = xup[i]
    ynaca[imax-1 + i] = yup[i]

#calibrate to zero
ynaca[0] = 0.0
ynaca[tot_nodes-1] = 0.0
    
z = np.zeros_like(xnaca)

tk = str(t)
series = str( int(100*m) )+str( int(10.0*p) )+tk[2]+tk[3]
nodes = str(tot_nodes)
filename ="naca"+series+"-"+nodes.zfill(3)+".csv"
df = pd.DataFrame({'x':xnaca.flatten(), 'y':ynaca.flatten(), 'z':z.flatten()} )
df.to_csv(filename, index=False)

xmin = -0.1
xmax = 1.1
ymin = -0.3
ymax = 0.3
plt.axis( (xmin, xmax, ymin, ymax) )
plt.plot( xnaca.flatten(), ynaca.flatten(), 'bo-')
plt.show()
