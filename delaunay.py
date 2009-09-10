from numpy import *
from structures import *

def voronoi(x,y):
	n, = x.shape
	p0 = array([x[0], y[0]])
	p1 = array([x[1], y[1]])
	
	L = Line(p0,p1).bisector()
	
	T0 = Tile(array([L.p0x, L.p1x]), array([L.p0y, L.p1y]), center = array([x[0], y[0]]), infinite = True)
	T1 = Tile(array([L.p1x, L.p0x]), array([L.p1y, L.p0y]), center = array([x[1], y[1]]), infinite = True, neighbours = array([T0]))
	T0.neighbours = array([T1])
	
	T0.plot()
	T1.plot()
	
	#for i in range(0,n):
