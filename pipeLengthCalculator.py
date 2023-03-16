import math
import matplotlib.pyplot as plt
import pandas as pd

def layLoops(vz,wa,la,l_edge):
    d=2*vz
    wlmin = min(wa, la)
    nLoops2 = math.ceil(wlmin/d)
    nLoops = math.ceil(wlmin/(2*d))
    lastLoop = nLoops2 >= 2*nLoops
    rest=nLoops*2*d - wlmin
    if rest > 3*vz:
        lastLoop = True
        nLoops -= 1
    elif rest > vz:
        lastLoop = False
    nPoints = 4*nLoops
    if lastLoop:
        nPoints += 4*nLoops
    else:
        nPoints += 4*(nLoops-1)
    x=[0]*nPoints
    y=[0]*nPoints
    for i in range(0,nLoops):
        x[4*i] = i*d
        y[4*i] = (i-1)*d
        x[4*i+1] = i*d
        y[4*i+1] = la-i*d
        x[4*i+2] = wa-i*d
        y[4*i+2] = la-i*d
        x[4*i+3] = wa-i*d
        y[4*i+3] = i*d
        # back
        if lastLoop or i < nLoops-1:
            x[nPoints-1-(4*i)] = i*d+vz
            y[nPoints-1-(4*i)] = (i-1)*d+vz
            x[nPoints-1-(4*i+1)] = i*d+vz
            y[nPoints-1-(4*i+1)] = la-i*d-vz
            x[nPoints-1-(4*i+2)] = wa-i*d-vz
            y[nPoints-1-(4*i+2)] = la-i*d-vz
            x[nPoints-1-(4*i+3)] = wa-i*d-vz
            y[nPoints-1-(4*i+3)] = i*d+vz
            
    y[0] += d-l_edge
    if lastLoop:
        x[4*nLoops] = x[4*nLoops+3]+vz
        y[4*nLoops] -= vz
    else:
        x[4*nLoops-1] = x[4*nLoops-4]+vz
        y[4*nLoops-1] -= vz
    y[nPoints-1] += vz-l_edge
    return x,y,4*nLoops

def calcPipeLength(x,y,r):
    length = 0
    for i in range(len(x)-1):
        dx = x[i+1]-x[i]
        dy = y[i+1]-y[i]
        partLength = math.sqrt(pow(dx,2)+pow(dy,2))
        length = length + partLength
    length = length+(len(x)-2)*(-2*r+math.pi*r/2) # quarter perimeter per point
    return length


# Parameters
w_edge =  0.2 # space to wall in x-direction
l_edge =  0.2 # space to wall in y-direction
vz = 0.1 # m pipe distance
r = 0.3 if(vz>=0.2) else 0.1 # m radius of the pipe edges
wa = 3.9-2*w_edge # space in x-direction w/o edges (total width of the heat exchanger)
la = 2.3-2*l_edge # space in y-direction w/o edges (total length of the heat exchanger)

# Calculation
[x,y,rlStart] = layLoops(vz,wa,la,l_edge)
print("Length of pipe = " + str(round(calcPipeLength(x,y,r),2)) + " m")

"""
# Plot
fig, ax = plt.subplots()
ax.plot(x[0:rlStart], y[0:rlStart], color="red")
ax.plot(x[rlStart-1:rlStart+1], y[rlStart-1:rlStart+1], color="black")
ax.plot(x[rlStart:], y[rlStart:], color="blue")
ax.set(aspect=1)
plt.show()

# Save points in a table
table = pd.DataFrame({"x":x,"y":y})
table.to_csv("layLoops_xy_d"+str(vz)+".csv",index=False)
"""
