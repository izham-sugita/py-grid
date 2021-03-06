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

#imax is the NACA airfoil total surface coordinates
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

#calibrate to tail edge reference point of (1.0, 0.0)
xnaca[0] = 1.0
ynaca[0] = 0.0

xnaca[tot_nodes-1] = 1.0
ynaca[tot_nodes-1] = 0.0
    
z = np.zeros_like(xnaca)

#Tail grid BC
#The tail edge zone total nodes can be increased
#in this setting tailmax = 2*imax
coeff = 1.5
tailmax = int(coeff*imax)

#xt = np.ndarray( (imax) )
#yt = np.ndarray( (imax) )

xt = np.ndarray( (tailmax) )
yt = np.ndarray( (tailmax) )

xu = np.zeros_like(xt)
yu = np.zeros_like(yt)

#temp = np.ndarray((imax))
#temp = np.ndarray( (tailmax) )
temp = np.ndarray( (len(xt)) )

tail_pos_x = xnaca[0]
tail_pos_y = ynaca[0]
dom_end = 20.0*c
length = dom_end - tail_pos_x
dl = length / (tailmax-1) #don't change this position    


#test sanity
#stretching parameter range 0<=temp<=1.0
# xb and xe were selected just to avoid division by zero
xb = 1.0
xe = 100.0
dls = (xe-xb) / (tailmax-1)
for i in range(len(temp)):
    ii = len(temp)-1 - i #reverse direction
    temp[i] =  np.log( (xe/xb)/( (ii*dls+xb)/ xb )  ) 
    yt[i] = 0.0

#stdize
std_one = max(temp)
for i in range(len(temp)):
    temp[i] = (temp[i]/std_one)*length + tail_pos_x 

#df = pd.DataFrame( {"x":temp, "y":yt} )
#filedebug="stretched-debug.csv"
#df.to_csv(filedebug, index=False)

#original tail loop
for i in range( len(temp) ):
    #temp[i] = i*dl + tail_pos_x #uniform distribution
    yt[i] = 0.0

xt = np.flip(temp)
xu[:] = temp[:]
yu[:] = yt[:]

for i in range( len(xt) ):
    print(xt[i], xu[i])

#joining the array
# xt + xnaca + xu
# yt + ynaca + yu

testline_x = np.concatenate( (xt[0:tailmax-1], \
                              xnaca, xu[1:tailmax]) )

testline_y = np.concatenate( (yt[0:tailmax-1],  \
                              ynaca, yu[1:tailmax]) )

print(len(testline_x))
bound_nodes = len(testline_x)

#Debug data
#df = pd.DataFrame( {'x': testline_x.flatten(), 'y':testline_y.flatten()} )
#filedebug = "line-debug.csv"
#df.to_csv(filedebug, index=False)

#tail-edge BC nodes index
#very important for boundary condition setting
#use list for the index
tail_bc_index = []
i = 0
while i < bound_nodes:
    if testline_x[i] == 1.0 and testline_y[i] == 0.0:
        tail_bc_index.append(i)
    i +=1

print(tail_bc_index)


#outer boundary
print(len(xnaca), len(ynaca))
tail_pos_x = xnaca[0]
tail_pos_y = ynaca[0]
dseta = np.pi /( len(xnaca)-1 )
semi_circle_x = np.ndarray( (len(xnaca)) )
semi_circle_y = np.ndarray( (len(xnaca)) )
rad = 10.0*c
for i in range(len(xnaca)):
    semi_circle_x[i] = rad*np.cos( 1.5*np.pi - i*dseta ) + tail_pos_x
    semi_circle_y[i] = rad*np.sin( 1.5*np.pi - i*dseta ) + tail_pos_y 

#df = pd.DataFrame({"x": semi_circle_x, "y": semi_circle_y})
#filedebug = "semi_circle-debug.csv"
#df.to_csv(filedebug, index=False)

#xt, xu, yt, yu
#the outer boundary line x-coordinate = inner x-coordinate
#only the value of yt and yu
#with yt being the lower section, yt = const. = - rad (= -10.0c )
#and yu = const. = + rad (= 10.0c)

outer_low_x = np.ndarray( (len(xt)) )
outer_up_x = np.ndarray( (len(xu)) )
outer_low_y = np.ndarray( (len(yt) ) )
outer_up_y = np.ndarray( (len(yu)) )

#y-coordinate
outer_low_y[:] = -rad
outer_up_y[:] = rad

#x-coordinate
outer_low_x[:] = xt[:]
outer_up_x[:] = xu[:]

#concatenate lines
outer_x_line = np.concatenate( (outer_low_x[0:tailmax-1],  \
                                semi_circle_x, outer_up_x[1:tailmax]) )

outer_y_line = np.concatenate( (outer_low_y[0:tailmax-1],  \
                                semi_circle_y, outer_up_y[1:tailmax]) )

print("Outer boundary total nodes:")
print( len(outer_x_line) )

#sanity check
sanity = len(outer_x_line) - len(testline_x)
if sanity == 0:
    print("Boundary line is consistent!")

#closure section xt = xu = 20.0c
#or at testline_x[0] = testline_x[bound_nodes-1]
#need to set the j-index
# i-index range: 0<= i <= bound_nodes-1
#temporary jmax = bound_nodes
jmax = int( 0.25*bound_nodes )
print("Default jmax: ", jmax)
ans = input("Change default value? y/n ")
if ans == 'y':
    jmax = int( input("Enter jmax ") )
else:
    print("using default jmax: ", jmax)

print("xi-axis total nodes: ", bound_nodes)
print("eta-axis total nodes: ", jmax)

#setting up the coordinates array
cl_low_x = np.ndarray((jmax))
cl_low_y = np.ndarray((jmax))

cl_up_x = np.ndarray((jmax))
cl_up_y = np.ndarray((jmax))

#x-coordinate is constant, testline_x[0] = testline_x[bound_nodes-1]
cl_low_x[:] = testline_x[0]
cl_up_x[:] = testline_x[0]

#y-spacing for cl_low_y and cl_up_y
dy = (10.0*c) / (jmax-1)
for j in range(jmax):
    cl_low_y[j] = -j*dy + testline_y[0]
    cl_up_y[j] = j*dy  + testline_y[0]

#df = pd.DataFrame( {"x": cl_low_x, "y": cl_low_y} )
#filedebug = "cl_low_line_debug.csv"
#df.to_csv(filedebug, index=False)
    
tk = str(t)
series = str( int(100*m) )+str( int(10.0*p) )+tk[2]+tk[3]
nodes = str(tot_nodes)
filename ="naca"+series+"-"+nodes.zfill(3)+"-C-grid-init.csv"
#df = pd.DataFrame({'x':xnaca.flatten(), 'y':ynaca.flatten(), 'z':z.flatten()} )
#df.to_csv(filename, index=False)

xmin = -1.0
xmax = 20.5
ymin = -0.3
ymax = 0.3

#plt.axis( (xmin, xmax, ymin, ymax) )
#plt.plot( xnaca.flatten(), ynaca.flatten(), 'bo-')

plt.plot( testline_x.flatten(), testline_y.flatten(), 'bo-')

#plt.plot( semi_circle_x.flatten(), semi_circle_y.flatten(), 'k-')

plt.plot( outer_x_line.flatten(), outer_y_line.flatten(), 'ro-')

plt.plot( cl_low_x.flatten(), cl_low_y.flatten(),  'ko-')
plt.plot( cl_up_x.flatten(), cl_up_y.flatten(),  'ko-')

plt.show()

#final output requires (testline_x, testline_y) -> inner boundary
#and (outer_x_line, outer_y_line) ->outer boundary.
#the in-between output can be use just as a reference or for elliptic PDE
#initial condition.
# closure line (cl_low_x, cl_low_y)
# closure line (cl_up_x, cl_up_y)

#xi =bound_nodes
#eta =jmax

x2d = np.ndarray( (bound_nodes, jmax) )
y2d = np.ndarray( (bound_nodes, jmax) )

#data sanity check
x2d[:][:] = 0.5
y2d[:][:] = 0.0

for i in range(bound_nodes):
    j = 0
    x2d[i][j] = testline_x[i]
    y2d[i][j] = testline_y[i]

    j = jmax-1
    x2d[i][j] = outer_x_line[i]
    y2d[i][j] = outer_y_line[i]

for j in range(jmax):
    i = 0
    x2d[i][j] = cl_low_x[j]
    y2d[i][j] = cl_low_y[j]

    i = bound_nodes-1
    x2d[i][j] = cl_up_x[j]
    y2d[i][j] = cl_up_y[j]

df = pd.DataFrame( {"x": x2d.flatten(), "y": y2d.flatten() } )
str_imax = str(bound_nodes)
str_jmax = str(jmax)
str_tail_1 = str(tail_bc_index[0])
str_tail_2 = str(tail_bc_index[1])

#stdizing
str_imax = str_imax.zfill(3)
str_jmax = str_jmax.zfill(3)
str_tail_1 = str_tail_1.zfill(3)
str_tail_2 = str_tail_2.zfill(3)

fileout = "naca"+series_number+"-"+str_imax+"-"+str_jmax+"-" \
          +str_tail_1+"-"+str_tail_2+ "-input.csv"

df.to_csv(fileout, index=False)
