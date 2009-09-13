from numpy import *
from pylab import *

def point_sort(z, coord = 0):
	""" Sort the points in z according to one of their coordinates, the first one by default. """
	arg = argsort(z[coord])
	sorted = take(z, arg, 1)
	return sorted

def collinear(p0, p1, p2):
	if allclose(p0,p1) or allclose(p0,p2) or allclose(p1,p2):
		return True
	v1 = p1-p0
	v1 = v1/norm(v1)
	v2 = p2-p0
	v2 = v2/norm(v2)
	return allclose(v1, v2)

def intersect(l1, l2, s0 = 0, s1 = 1, t0 = 0, t1 = 1):
	""" Check whether two lines intersect between the points that define them. """
	A = c_[l1.p1-l1.p0, l2.p0-l2.p1]
	b = l2.p0-l1.p0
	try:
		(s,t) = solve(A,b)
	except LinAlgError:
		# Singular matrix => same line or no intersection
		# Check if they are the same:
		if collinear(l1.p0, l1.p1, l2.p0):
			return (0,0)
		else:
			return (None,None)
	if (s < s0) or (s > s1) or (t < t0) or (t > t1):
		return (None,None)
	return (s,t)


#
#class Point:
#	def __init__(self, x, y, z='2D'):
#		self.x = x
#		self.y = y
#		if z == '2D':
#			self.dim = 2
#			self.z = 0
#		else:
#			self.dim = 3
#			self.z = z
#	
#	def __getitem__(self, i):
#		if i == 0:
#			return self.x
#		elif i == 1:
#			return self.y
#		elif i == 2 and self.dim == 3:
#			return self.z
#		else:
#			raise IndexError('Invalid Index	')
#		
#	def a(self):
#		if self.dim == 2:
#			return array([self.x, self.y])
#		else:
#			return array([self.x, self.y, self.	z])

class Line:
	""" Defines a line through the points p0 and p1. """
	def __init__(self, p0, p1):
		self.p0 = p0
		self.p0x = p0[0]
		self.p0y = p0[1]
		self.p1 = p1
		self.p1x = p1[0]
		self.p1y = p1[1]

	def begin(self):
		return self.p0
	def end(self):
		return self.p1

	def length(self):
		return norm(self.p1-self.p0)

	def interpolate(self, t):
		return array([(1-t)*self.p0x + t*self.p1x,  (1-t)*self.p0y + t*self.p1y])

	def points(self,num_points = 50, t0 = 0, t1 = 1):
		t = linspace(t0,t1,num_points)
		p = zeros((num_points, self.p0.size))
		for i in range(0,num_points):
			p[i] = self.interpolate(t[i])
		return p

	def projection(self,p):
		""" Return the orthogonal projection of the point p onto this line. """
		v = self.p1-self.p0
		v = v/norm(v)
		p = p - self.p0
		return v*dot(v,p)+self.p0

	def bisector(self):
		""" Return the perpendicular bisector. (It will have the same length as this line.) """
		midpoint = 0.5*(self.p0 + self.p1)
		dir = self.p1-self.p0
		bisector_dir = array([dir[1], -dir[0]])
		p0new = midpoint - 0.5*bisector_dir
		p1new = midpoint + 0.5*bisector_dir
		return Line(p0new, p1new)

	def plot(self, t0 = 0, t1 = 1, color = 'k'):
		p0 = self.interpolate(t0)
		p1 = self.interpolate(t1)
		plot((p0[0],p1[0]), (p0[1],p1[1]), color)


class Triangle:
	""" Defines a triangle with vertices p0, p1 and p2. Neighbouring triangles are stored in the variable neighbours. """
	def __init__(self, p0, p1, p2, neighbours = None):
		self.p0 = p0
		self.p1 = p1
		self.p2 = p2
		self.neighbours = neighbours
	def edge0(self):
		return Line(self.p0, self.p1)
	def edge1(self):
		return Line(self.p1, self.p2)
	def edge2(self):
		return Line(self.p2, self.p0)
	def area(self):
		""" Calculate the area of the triangle. """
		v0 = self.p1-self.p0
		v1 = self.p2-self.p0
		return norm(cross(v0,v1))/2.
	
	def plot(self, color = 'k'):
		self.edge0().plot()
		self.edge1().plot()
		self.edge2().plot()
		
class Polygon:
	""" Defines a convex polygon with n corners. The x- and y-coordinates for each corner are stored
		in p[0,:] and p[1,:] respectively. """
	def __init__(self, p):
		self.p = p #c_[x,y].T
		self.n = p.shape[1] #x.size

	def x(self):
		return self.p[0]
	def y(self):
		return self.p[1]

	def edge(self,i):
		""" Return the i:th edge, i.e. the line between the i:th corner and the (i+1):th corner. """
		return Line(self.p[:,i], self.p[:,i+1])
	
	def add_point(self, pnew, i):
		""" Add the point pnew after the i:th corner and before the (i+1):th corner. """
		if self.n == 0:
			self.p = pnew
		else:		
			self.p = c_[self.p[:,:i], pnew, self.p[:,i:]]
		self.n += 1
		
	def remove_point(self, i):
		""" Remove the i:th corner. """
		self.p = delete(self.p, i, 1)
		self.n -= 1

	def xmin(self):
		return min(self.p[0])
	def xmax(self):
		return max(self.p[0])
	def ymin(self):
		return min(self.p[1])
	def ymax(self):
		return max(self.p[1])

	def in_bbox(self, p):
		if (p[0] < self.xmin()) or (p[0] > self.xmax()) or (p[1] < self.ymin()) or (p[1] > self.ymax()):
			return False
		else:
			return True

	def contains(self, p):
		""" Check whether the polygon contains the point p or not. """
		# We check for each edge in the polygon whether p is "to the right" of the edge
		# by checking if the crossproduct of the vector from one corner to p and the vector from
		# the same corner to the next corner is positive.
		for i in range(0, self.n-1):
			if cross(p-self.p[:,i], self.p[:,i+1]-self.p[:,i]) < 0:
				return False
		if cross(p-self.p[:,-1], self.p[:,0]-self.p[:,-1]) < 0:
			return False
		return True

	def area(self):
		""" Calculate the area of the polygon. """
		# We divide the polygon into n-2 triangles and sum their areas
		A = 0
		for i in range(2, self.n):
			t = Triangle(self.p[:,0], self.p[:,i-1], self.p[:,i])
			A += t.area()
		return A
		
	def plot(self, color = 'k'):
		for i in range(0,self.n-1):
			Line(self.p[:,i], self.p[:,i+1]).plot(color = color)
		Line(self.p[:,-1], self.p[:,0]).plot(color = color)

class Tile(Polygon):
	""" Defines a tile for Delaunay triangulations. Center is the defining point of the Tile.
		If the argument infinite is true then the first and last point are not connected and the
		relevant lines are assumed to continue indefinitely. """
	def __init__(self, p, center, infinite, neighbours = []):
		Polygon.__init__(self, p)
		self.center = center
		self.infinite = infinite
		self.neighbours = neighbours

	def intersections(self, L):
		""" Find the intersections p, q of the line L and this tile. Returns (None, None) if there 
			is no intersection, (p, None) if there is just one intersection and (p,q) otherwise. """
		(p,q) = (None, None)
		(ip, iq) = (None, None)
		for i in range(1, self.n-2):
			e = self.edge(i)
			(s,t) = intersect(L, e, s0 = -Inf, s1 = Inf)
			if s != None:
				if p == None:
					#p = L.interpolate(s)
					p = (s,t)
					ip = i
				else:
					#q = L.interpolate(s)
					q = (s,t)
					iq = i
					break

		if self.infinite:
			if self.n > 2:
				for i in [0, self.n-2]:
					e = self.edge(i)
					if i == 0:
						(s,t) = intersect(L, e, s0 = -Inf, s1 = Inf, t0 = -Inf, t1 = 1)
					else:
						(s,t) = intersect(L, e, s0 = -Inf, s1 = Inf, t0 = 0, t1 = Inf)
					if s != None:
						if p == None:
							p = (s,t)
							ip = i
						else:
							q = (s,t)
							iq = i
							break
			else:
				e = self.edge(0)
				(s,t) = intersect(L, e, s0 = -Inf, s1 = Inf, t0 = -Inf, t1 = Inf)
				if s != None:
					if p == None:
						p = (s,t)
						ip = 0
					else:
						q = (s,t)
						iq = 0
		return (p,q,ip,iq)
			

	def plot(self, color = 'k'):
	#	plot(self.center[0], self.center[1], 'ro')
		for i in range(1,self.p.shape[1]-1):
			Line(self.p[:,i], self.p[:,i+1]).plot(color = color)
	#	plot(self.p[0,:], self.p[1,:],'gx')
		if not self.infinite:
			Line(self.p[:,0], self.p[:,1]).plot(color = color)
			Line(self.p[:,-1], self.p[:,0]).plot(color = color)
		else:
			Line(self.p[:,0], self.p[:,1]).plot(t0 = -10, color = color)
			Line(self.p[:,-2], self.p[:,-1]).plot(t1 = 10, color = color)			