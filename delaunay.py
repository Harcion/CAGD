from numpy import *
from structures import *

#x = array([0,0,1,2])
#y = array([0,1,1,0])
def voronoi(x,y):
	N, = x.shape
	P = c_[x,y].T

	plot(x,y,'go')

	p0 = P[:,0]
	p1 = P[:,1]

	L = Line(p0,p1).bisector()

	T0 = Tile(c_[L.p0, L.p1], center = p0, infinite = True)
	T1 = Tile(c_[L.p1, L.p0], center = p1, infinite = True,neighbours = [T0])

	T0.neighbours.append(T1)

	T0.plot()
	T1.plot()

	T = [T0,T1]
	
	if N < 3:
		return
	
	p = P[:,2]
	done = False
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
	L.plot()
	# Then find the intersection between the boundary of Ti and this line
	(p1, p2, i1, i2) = Ti.intersections(L)
	(s1,t1) = p1
	if p2 == None:
		infinite = True
		edge = Ti.edge(i1)
		s = s1
		t = t1
		# Check if p is to the left or right of Ti.center
		Pc = edge.projection(Ti.center)
		Pp = edge.projection(p)
		tc = (Pc[0]-edge.p0[0])/(edge.p1[0]-edge.p0[0])
		tp = (Pp[0]-edge.p0[0])/(edge.p1[0]-edge.p0[0])
		#  If tp > tc then p is to the right of Ti.center. We cannot have equality if p != Ti.center.
		# Left is the bad case
		start_index = 0
		if tp > tc:
			orientation = 'right'
			e1 = L.interpolate(s)
			Tnew = Tile(c_[L.p0, e1], center = p, infinite = True, neighbours = [Ti])
		else:
			orientation = 'left'
			start_index = -1
			for k in range(len(Ti.neighbours)-1, -1,-1 ):  # Goes from last Ti.neighbour to 0
				Tk = Ti.neighbours[k]
				Lk = Line(p, Tk.center).bisector()
				(p1, p2, i1, i2) = Tk.intersections(Lk)
				if p1 == None:
					# Go to the next neighbour and check for intersection there instead
					k += 1
				else:
					start_index = k  # We will start the walk-around in this tile
					(s,t) = p1
					s0 = min(-1, -2*abs(s))
					e0 = Lk.interpolate(s0)
					e1 = Lk.interpolate(s)
					Tnew = Tile(c_[e0,e1], center = p, infinite = True, neighbours = [Tk])
					Lk.plot()
					print Tk is T0
					break
	else:
		# Pick the "right" edge if there are two intersections (i.e. not left)
		# This way the orientation will be 'right' automatically
		orientation = 'right'
		(s2,t2) = p2
		if s2 > s1:
			edge = Ti.edge(i2)
			s = s2
		else:
			edge = Ti.edge(i1)
			s = s1
		
	
	e1 = L.interpolate(s)  # Intersection
#	m = (Ti.center + p)/2.  # Other point on the new edge
	
#	
#	if tp > tc:
#		e2 = edge.interpolate(max(2*abs(t),1))
#		Tnew = Tile(c_[e1,m,e2], center = p, infinite = True, neighbours = [Ti])
#		if t < 0:
#			Ti.p[:,1] = e1
#			Ti.p[:,0] = edge.interpolate(2*t)
#		elif t < 1:
#			Ti.p[:,1] = e1
#		else:
#			Ti.p[:,0] = e1
#		Ti.add_point(m,Ti.n)
#	else:
#		e2 = edge.interpolate(min(-2*abs(t),-1))
#		Tnew = Tile(c_[e2,e1,m], center = p, infinite = True, neighbours = [Ti])
#		if t < 1:
#			Ti.p[:,0] = e1
#		else:
#			Ti.p[:,0] = e1
#			Ti.p[:,1] = edge.interpolate(2*t)
#		Ti.add_point(m,0)
	
	ek = copy(e1)
	
#	if done:
#		start_index = len(Ti.neighbours)+1 # Skip the for-loop

	for k in range(start_index,len(Ti.neighbours)):
		Tk = Ti.neighbours[k]

		Lk = Line(p, Tk.center).bisector()
	#	Lk.plot()
		
		(p1, p2, i1, i2) = Tk.intersections(Lk)
		if p2 == None:
			# Go to the next neighbour and check for intersection there instead
			k += 1
		else:
			(s1,t1) = p1
			(s2,t2) = p2
			p1 = Lk.interpolate(s1)
			p2 = Lk.interpolate(s2)
			if allclose(ek, p1):
				ek = p2
			else:
				ek = p1
			Tnew.add_point(ek, Tnew.n)
			Tnew.neighbours.append(Tk)

	if not infinite: # Then we should reach Ti again
		Lk = Line(p, Tk.center).bisector()
		(p1, p2, i1, i2) = Ti.intersections(Lk)
		(s,t) = p2
		p2 = Lk.interpolate(s)
		Tnew.add_point(p2, Tnew.n)
	else:
		# Add a point along the last line which is infinite
	#	Lk.plot()		
		if orientation == 'right':
			Tnew.add_point(Lk.p1, Tnew.n)
			# Also add a point along the first line to make it a line
			Tnew.add_point(L.p0, 0)
		else:
			Tnew.add_point(L.p1, Tnew.n)
	
	figure()
	plot(x,y,'go')
	Tnew.plot()
	print Tnew.p
#		
#	for Tk in Ti.neighbours:
#		L2 = Line(p, Tk.center).bisector()
#		L2.plot()
#		m = (Tk.center + p)/2.  # Other point on the new edge
#		print m
#
#		if tp > tc:
#			if allclose(m, e1):
#				m = L2.p0
#				print m
#			T2.p[:,0] = m
#			Tk.p[:,-1] = e1
#			Tk.add_point(m, Tk.n)
#			
#		else:
#			if allclose(m, e1):
#				m = L2.p1
#			T2.p[:,-1] = m
#			Tk.p[:,0] = e1
#			Tk.add_point(m,0)
#	
#	T.append(T2)
#	figure()
#	plot(x,y,'go')
#	for Tk in T:
#		#figure()
#		Tk.plot()
#	i = 0
#	for Tk in T:
#		figure()
#		plot(x,y,'go')
#		Tk.plot()
#		title(i)
#		i += 1
#	e2 = L.p0				# One more point on the line
#
#	if allclose(e2, e1):	# If coincide, pick the other one which is farther away
#		e2 = L.p1
#	plot(e2[0],e2[1],'bx')		
#	print Ti.p
#	# Check where the intersection is in relation to edge.p0, edge.p1 and move the appropriate point
#	t = p1[1]
#	print t
#	if t < 0:
#		# Adjust T2 edge
#		T2.p[:,1] = e1
#		T2.p[:,0] = T2.edge(0).interpolate(2*t)
#		# Add new T2 edge
#		T2.add_point(e2, T2.n)
#		# Adjust Ti edges
#		Ti.p[:,0] = e1
#		Ti.add_point(e2,0)
#	elif t < 1:
#		# Adjust T2 edge
#		T2.p[:,1] = e1
#		# Add new T2 edge
#		T2.add_point(e2, T2.n)
#		print Ti.p
#		# Adjust Ti edges
#		Ti.p[:,0] = e1
#		Ti.add_point(e2,0)
#		print "here"
#		print Ti.p
#	else:
#		# Adjust T2 edge
#		T2.p[:,1] = e1
#		# Add new T2 edge
#		T2.add_point(e2, T2.n)
#		# Adjust Ti edges
#		Ti.p[:,0] = e1
#		Ti.p[:,1] = Ti.edge(0).interpolate(2*t)
#		Ti.add_point(e2,0)
#	
#	
#	Tk = Ti.neighbours[0]
#	# Find the perpendicular bisector to the line between p and Ti's center
#	L2 = Line(p, Tk.center).bisector()
#	# Will have same intersection e1 so just pick the midpoint to get another point on the bisector
#	e2 = (Tk.center+p)/2.
#	if allclose(e1,e2):
#		e2 = Line(p,Tk.center).bisector().p0
#	# Adjust T2 edge
#	T2.p[:,0] = e2
#	if t < 0:
#		# Adjust Tk edge
#		Tk.p[:,1] = e1
#		Tk.add_point(e2,Tk.n)
#	elif t < 1:
#		# Adjust Tk edge
#		Tk.p[:,1] = e1
#		Tk.add_point(e2,Tk.n)
#	else:
#		# Adjust Tk edge
#		Tk.p[:,1] = e1
#		Tk.p[:,0] = Ti.edge(0).interpolate(2*t)
#		Tk.add_point(e2,Tk.n)
#
#	T2.plot()
#
#	T.append(T2)
#
#	#figure()
#	for Tk in T:
#		figure()
#		Tk.plot()
#	#	
#


