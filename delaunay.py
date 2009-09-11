from numpy import *
from structures import *

#x = array([0,0,1,2])
#y = array([0,1,1,0])
def voronoi(x,y):
	n, = x.shape
	P = c_[x,y].T

	plot(x,y,'go')

	p0 = P[:,0]
	p1 = P[:,1]

	L = Line(p0,p1).bisector()

	T0 = Tile(c_[L.p0, L.p1], center = p0, infinite = (True,True))
	T1 = Tile(c_[L.p1, L.p0], center = p1, infinite = (True,True),neighbours = [T0])

	T0.neighbours.append(T1)

	T0.plot()
	T1.plot()

	T = [T0,T1]

	p = P[:,2]

	# Find which tile p is in by checking which tile's center it is closest to
	r = Inf
	Ti = None
	for Tk in T:
		rk = norm(Tk.center-p)
		if rk < r:
			r = rk
			Ti = Tk
	# Find the perpendicular bisector to the line between p and Ti's center
	L = Line(p, Ti.center).bisector()
	# Then find the two intersections between the boundary of Ti and this line (or one intersection if Ti
	# is infinite)
	(s1, s2, i1, i2) = Ti.intersections(L)
	if s2 == None:
		q1 = L.interpolate(s1)
		e1 = copy(q1)
		e2 = L.p0
		if allclose(e2, q1):
			e2 = L.p1

		T2 = Tile(c_[e1, e2], p, infinite = (False,True), neighbours = copy(Ti))
		
		# Adjust Ti edges
		Ti.p[:,i1] = q1
		Ti.add_point(e2,0)
		Ti.plot()
		
		for Tk in Ti.neighbours:
			# Find the perpendicular bisector to the line between p and Ti's center
			L2 = Line(p, Tk.center).bisector()
			# Check if intersect...
			(s1, s2, i1, i2) = Tk.intersections(L2)
			if s1 != None:
				q1 = L2.interpolate(s1)
				if allclose(q1,Tk.p[:,0]):
					T2.add_point(L2.p0,0)
				else:
					e1 = copy(q1)
					if s1 < 0: # L2.p0 is on the wrong side of q1
						e2 = L2.interpolate(2*s1)
					else:
						e2 = L2.p0
					#e2 = Line(q1,L2.p0).interpolate(-1)
					T2.add_point(e1, 0)
					T2.add_point(e2, 0)
					print Tk.p
					Tk.p[:,i1+1] = q1
					Tk.add_point(e2,Tk.n)
			

	T2.plot()

	T.append(T2)

	#figure()
	for Tk in T:
		figure()
		Tk.plot()
	#	
#
##	for Tk in Ti.neighbours:
##		
###	Line(T0.center, p).plot()
##		
##	for Tk in T:
##		Line(Tk.center, p).bisector().plot()
##	
##	#for i in range(0,n):
