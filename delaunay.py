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
		return T
	
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
	# Then find the intersection between the boundary of Ti and this line
	(p1, p2, i1, i2) = Ti.intersections(L)
	if p1 == None: # No intersections - parallell!
		pass
	elif p2 == None: # One intersection
		infinite = True
		(s,t) = p1
		# Check which neighbour is on the other side of the edge
		r = Inf
		for k in range(0, len(Ti.neighbours)):
			rk = norm(Ti.neighbours[k].center - e1)
			if rk < r:
				rk = r
				start_index = k

		# Check if p is to the left or right of Ti.center
		edge = Ti.edge(i1)
		Pc = edge.projection(Ti.center)
		Pp = edge.projection(p)
		tc = (Pc[0]-edge.p0[0])/(edge.p1[0]-edge.p0[0])
		tp = (Pp[0]-edge.p0[0])/(edge.p1[0]-edge.p0[0])
		# If tp > tc then p is to the right of Ti.center. We cannot have equality if p != Ti.center.
		# Left is the "bad" case, we cannot move clockwise at once since that will take us to infinity
		# (and beyond), so we have to rewind and find the first tile having just one intersection with
		# the line as well, this will be the one where we are coming back from infinity and the one
		# we should start moving clockwise from.
		if tp > tc:
			orientation = 'right'
			e1 = L.interpolate(s)
			Tnew = Tile(c_[L.p0, e1], center = p, infinite = True, neighbours = [Ti])
			
			#TODO: Adjust Ti edges
		else:
			orientation = 'left'
			start_index = -1
			index = range(start_index,len(Ti.neighbours))
			index.extend(0, start_index)
			for k in index:  # This will go backwards from the first neighbour past the one we should end at
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
					break
			#TODO: Adjust Ti edges
	else:
		# Two intersections
		# Pick the "right" edge to start with (i.e. not left)
		# This way the orientation will be 'right' automatically
		orientation = 'right'
		infinite = False
		(s1, t1) = p1
		(s2, t2) = p2
		p1 = L.interpolate(s1)
		p2 = L.interpolate(s2)
		if s2 > s1:
			e1 = p1
			e2 = p2
			edge1 = i1
			edge2 = i2
		else:
			e1 = p2
			e2 = p1
			edge1 = i2
			edge2 = i1
		Tnew = Tile(c_[e1,e2], center = p, infinite = False, neighbours = [Ti])
		# Adjust edges of Ti:
		Ti.p = c_[Ti.p[:,0:edge1], e2, e1, Ti.p[:,edge2+1:]]
		Ti.n = Ti.p.shape[1]
	

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
	
	index = range(start_index,len(Ti.neighbours))
	index.extend(range(0,start_index))
	for k in index:
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
			if allclose(e1,p1) or allclose(e1,p2):
				# We are back where we started
				break			
			if allclose(ek, p1):
				ek2 = p2
				edge1 = i2	# Index of the edge which the line intersects at ek2
				edge2 = i1	# Index of the edge which the line intersects at ek1				
			else:
				ek2 = p1
				edge1 = i1
				edge2 = i2
			Tnew.add_point(ek2, Tnew.n)
			Tnew.neighbours.append(Tk)
			
			# Adjust edges of Tk:
			# Edges with index between edge1 and edge2 should be removed and exchanged with one edge
			# from ek2 to ek. (Since the intersecting line is moving clockwise.)
			Tk.p = c_[Tk.p[:,0:edge1], ek2, ek, Tk.p[:,edge2+1:]]
			Tk.n = Tk.p.shape[1]
			
			ek = copy(ek2)

	if infinite:
		# Add a point along the last line, which is infinite
		if orientation == 'right':
			Tnew.add_point(Lk.p1, Tnew.n)
		else:
			Tnew.add_point(L.p1, Tnew.n)
	
	figure()
	plot(x,y,'go')
	Tnew.plot()
	print Tnew.p