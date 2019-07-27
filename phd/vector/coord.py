import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def coord_real(i,j,xmin, ymin, deltax, deltay):
    return xmin + j*deltax, ymin + i*deltay

def coord(i,j,xmin, ymin, deltax, deltay):
    return complex(xmin + j*deltax, ymin + i*deltay)

def mult(cpoint):
    x = cpoint.real
    y = cpoint.imag
    return x*y

def vel(cpoint):
    x = cpoint.real
    y = cpoint.imag
    return complex(x, y)


def plotvelocity(ar,br,vx,vy):
    X, Y = np.meshgrid(ar,br)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.quiver(X[::3,::3],Y[::3,::3], vx[::3,::3],vy[::3,::3])
    plt.show()
    plt.close()

def plotscalar(ar, br, a):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cs = ax.contourf(ar,br,a, num_lines)
    cs2 = ax.contour(ar,br,a, num_lines, linewidths=0.6, colors='black')
    ax.clabel(cs2, inline=1,fontsize=8)
    fig.colorbar(cs)
    plt.show()
    plt.close()


num_lines = 10
nx = 40
ny = 40

xmin = 10.
ymin = 10.
xmax = 20.
ymax = 25.
#include end points
deltax = (xmax - xmin)/(nx-1)
deltay = (ymax - ymin)/(ny-1)
a  = np.zeros([ny,nx])
vx = np.zeros([ny,nx])
vy = np.zeros([ny,nx])

for i in xrange(ny):
    for j in xrange(nx):
        print coord(i,j,xmin,ymin,deltax, deltay), i, j
        a[i,j] = mult(coord(i,j,xmin,ymin,deltax, deltay))
        vx[i,j] = vel(coord(i,j,xmin,ymin,deltax, deltay)).real
        vy[i,j] = vel(coord(i,j,xmin,ymin,deltax, deltay)).imag

ar = np.linspace(xmin,xmax,nx)
br = np.linspace(ymin,ymax,ny)

print nx, len(ar), ny, len(br)
    
print a

plotvelocity(ar,br,vx,vy)

plotscalar(ar,br,a)


#dont include endpoints
#deltax = (xmax - xmin)/(nx)
#deltay = (ymax - ymin)/(ny)
