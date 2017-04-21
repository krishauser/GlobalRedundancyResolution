import pkg_resources
pkg_resources.require("klampt>=0.6.2")
if pkg_resources.get_distribution("klampt").version >= '0.7':
	#Klampt v0.7.x
	from klampt.math import vectorops,so3,se3
	KLAMPT_VERSION = 0.7
else:
	#Klampt v0.6.x
	from klampt import vectorops,so3,se3
	KLAMPT_VERSION = 0.6
import math
import numpy as np
import sys

def toUtf8(input):
	"""Given a possibly hierarchical structure with unicode values (such as an object returned from json.load())
	converts all unicode strings to UTF-8 strings."""
	if isinstance(input, dict):
		return {toUtf8(key): toUtf8(value)
				for key, value in input.iteritems()}
	elif isinstance(input, list):
		return [toUtf8(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input

def linf_distance(a,b,spin=False):
	"""Return the L-infinity distance between configurations a and b where
	each joint is assumed to be revolute.  If spin=True, it assumes joints
	can spin 360 degrees."""
	if spin:
		return max(abs(v-w)%math.pi for (v,w) in zip(a,b))
	else:
		return max(abs(v-w) for (v,w) in zip(a,b))

def workspace_distance(a,b):
	"""Returns a distance in R^n or SE(3) between points a and b.  If len(a)<=3, 
	this will do cartesian distance.  Otherwise, len(a)==6 is required and
	this will do SE(3) distance (with a max orientation distance of 1)."""
	if len(a) <= 3:
		return vectorops.distance(a,b)
	assert len(a)==6,"Can only interpolate in R2, R3, or SE(3)"
	dt = vectorops.distance(a[:3],b[:3])
	dR = so3.distance(so3.from_moment(a[3:]),so3.from_moment(b[3:]))
	return dt + dR/math.pi

def workspace_interpolate(a,b,u):
	"""Interpolates either in R^n or SE(3) between points a and b.  If len(a)<=3, 
	this will do cartesian interpolation.  Otherwise, len(a)==6 is required and
	this will do SE(3) interpolation."""
	if len(a) <= 3:
		return vectorops.interpolate(a,b,u)
	assert len(a)==6,"Can only interpolate in R2, R3, or SE(3)"
	Ra,Rb = so3.from_moment(a[3:]),so3.from_moment(b[3:])
	R = so3.interpolate(Ra,Rb,u)
	return vectorops.interpolate(a[:3],b[:3],u) + so3.moment(R)

def mmult(A,b):
	"""Matrix-vector multiplication with raw Python objects"""
	v = []
	for ai in A:
		assert (len(ai) == len(b))
		v.append(vectorops.dot(ai,b))
	return v

def segment_point_distance(a,b,x):
	"""Returns (distance, parameter) pair between segment (a,b) and point x"""
	d = vectorops.sub(b,a)
	u = vectorops.dot(vectorops.sub(x,a),d)/vectorops.dot(d,d)
	u = min(1.0,max(u,0.0))
	d = vectorops.distance(x,vectorops.interpolate(a,b,u))
	return (d,u)

def segment_segment_intersection(seg1,seg2,umin=0,umax=1,vmin=0,vmax=1):
	"""Given 2 2D segments (a1,b1) and (a2,b2) returns whether the portion of
	seg1 in the range [umin,umax] overlaps seg2 in the range [vmin,vmax].

	Returns the point of intersection or None if no intersection exists
	"""
	a1,b1 = seg1
	a2,b2 = seg2
	assert len(a1)==2
	assert len(b1)==2
	assert len(a2)==2
	assert len(b2)==2
	d = vectorops.sub(b1,a1)
	a,b = vectorops.sub(a2,a1),vectorops.sub(b2,a1)
	#now it's centered at a1
	#u*d = a + v*(b-a)
	e = vectorops.sub(b2,a2)
	A = np.array([d,vectorops.mul(e,-1)]).T
	try:
		uv = np.dot(np.linalg.inv(A),a)
	except np.linalg.linalg.LinAlgError:
		return None;
	u,v = uv
	if u < umin or u > umax:
		return None
	if v < umin or v > umax:
		return None
	return vectorops.interpolate(a1,b1,u)

def segment_triangle_intersection(seg,tri,umin=0,umax=1):
	"""Given a 3D segment (p,q) and a triangle (a,b,c), returns the point of
	intersection if the region of the segment in interval [umin,umax] intersects the
	triangle. Otherwise returns None"""
	p,q = seg
	a,b,c = tri
	d = vectorops.sub(b,a)
	e = vectorops.sub(c,a)
	n = vectorops.cross(d,e)
	if vectorops.norm(n) < 1e-7: #degenerate
		return None
	n = vectorops.unit(n)
	ofs = vectorops.dot(n,a)
	#Solve for n^T(p + u(q-p)) = ofs
	denom = vectorops.dot(n,vectorops.sub(q,p))
	numer = ofs - vectorops.dot(n,p)
	if abs(denom) < 1e-7: #near parallel
		return None
	u = numer / denom
	if u < umin or u > umax:
		return None
	#now check whether point of intersection lies in triangle
	x = vectorops.madd(p,vectorops.sub(q,p),u)
	xloc = vectorops.sub(x,a)

	#solve for coordinates [r,s,t] such that xloc = r*n + s*d + t*e
	try:
		rst = np.dot(np.linalg.inv(np.array([n,d,e]).T),xloc)
	except np.linalg.linalg.LinAlgError:
		return None;
	r,s,t = rst
	assert abs(r) < 1e-4
	if s >= 0 and t >= 0 and s + t <= 1:
		return x
	return None

def robot_average(robot,configs,weights=None):
	"""Returns an "average" configuration from a set of configurations,
	somewhat properly taking into account spin joints and floating joints.
	"""
	assert len(configs) > 0
	qres = configs[0]
	if weights == None:
		for (i,q) in enumerate(configs):
			if i == 0: continue
			qres = robot.interpolate(qres,q,1.0/(i+1))
	else:
		assert len(weights) == len(configs)
		sumw = weights[0]
		for (i,(q,w)) in enumerate(zip(configs,weights)):
			if i == 0: continue
			if w == 0: continue
			qres = robot.interpolate(qres,q,w/(sumw+w))
			sumw += w
	return qres

def set_cover(sets,universe=None):
	"""
	Performs greedy set-cover.

	In: 
	- sets: a dictionary mapping set names to lists of items in the universe
	- universe: a list or set of items, or None if this is to be determined automatically.

	Return value: a dictionary mapping set names to the items that are *first* covered
	using the greedy algorithm.
	If a set is not included in the set cover, then it is not included in the result.
	"""
	newsets = dict()
	for k,v in sets.iteritems():
		if isinstance(v,(list,tuple)):
			newsets[k] = set(v)
		else:
			newsets[k] = v
	sets = newsets
	if universe == None:
		#auto-determine universe
		universe = set()
		for k,v in sets.iteritems():
			universe |= v
	else:
		universe = set(universe)
	#precompute element coverse
	covers = {}
	for e in universe:
		covers[e] = []
	for k,v in sets.iteritems():
		for e in v:
			assert e in covers,"Set "+str(k)+" contains element not in universe"
			covers[e].append(k)
	#compute sizes
	sizes = {}
	for k,v in sets.iteritems():
		sizes[k] = len(v)
	result = {}
	while len(universe) > 0:
		best = max((v,k) for (k,v) in sizes.iteritems())[1]
		covered = universe & sets[best]
		result[best] = covered
		del sizes[best]
		for e in covered:
			for s in covers[e]:
				if s in sizes:
					sizes[s] -= 1
		universe = universe - sets[best]
	return result

def ring(G,S,k=1):
	"""Computes the k-ring of the set of node S in the graph G.  Given an undirected graph
	and set of nodes S, the 1-ring of S in G is the set of nodes neighboring one node in S,
	not including S.  The k-ring is defined recursively as the 1-ring of the (k-1)-ring.

	The return value is a set. It also runs fastest when S is a set.
	"""
	if k>1:
		return ring(G,ring(G,S,k-1))
	res = set()
	for v in S:
		res |= set(G.neighbors(v))
	return res - set(S)

def triangles(G):
	"""Returns the list of all triangles (triples of mutually connected nodes) in G.
	Takes time O(nd^2) where n is the number of nodes in the graph and d is its degree"""
	res = []
	for i in G.nodes():
		Ni = set([j for j in G.neighbors(i) if j > i])
		for j in Ni:
			for k in Ni & set(G.neighbors(j)):
				if k < j: continue
				res.append((i,j,k))
	return res

def tetrahedra(G):
	"""Returns the list of all tetrahedra (4-tuples of mutually connected nodes) in G.
	Takes time O(nd^3) where n is the number of nodes in the graph and d is its degree"""
	res = []
	for i in G.nodes():
		Ni = set([j for j in G.neighbors(i) if j > i])
		for j in Ni:
			Nij = Ni & set([k for k in G.neighbors(j) if k > j])
			for k in Nij:
				for l in set(G.neighbors(k)) & Nij:
					if l < k: continue
					res.append((i,j,k,l))
	return res

class ProgressUpdater:
	"""A class that progressively prints out a progress bar to console ranging from 0 - 100%."""
	def __init__(self,n=None,frequency=1):
		self.cnt = 0
		self.n = n
		self.frequency = frequency
	def update(self):
		self.cnt += 1
		if self.n == None:
			if self.cnt % self.frequency == 0:
				print '  ...'
		else:
			if (100*self.cnt) / (self.n*self.frequency) != (100*(self.cnt-1)) / (self.n*self.frequency):
				print "  %d%%.."%((100*self.cnt) / self.n,),
				if self.cnt == self.n: print
				else: sys.stdout.flush()
	def done(self):
		pass


def test_set_cover():
	universe = range(14)
	sets = {'A':range(7),
		'B':range(7,14),
		'C':[0,1,2,3,7,8,9,10],
		'D':[4,5,11,12],
		'E':[6,13],
		'F':[1,4,7,12],
		'G':[16,17]}
	#print set_cover(sets,universe)
	print set_cover(sets,None)

if __name__ == '__main__':
	print "Testing set cover"
	test_set_cover()