from klampt.math import vectorops,so3
import math
import networkx as nx

def so3_grid(N):
	"""Returns a nx.Graph describing a grid on SO3 with resolution N. 
	The resulting graph has ~4N^3 nodes

	The node attribute 'params' gives the so3 element.
	"""
	assert N >= 1
	G = nx.Graph()
	R0 = so3.from_quaternion(vectorops.unit((1,1,1,1)))
	namemap = dict()
	for ax in range(4):
		for i in xrange(N+1):
			u = 2*float(i)/N - 1.0
			u = math.tan(u*math.pi*0.25)
			for j in xrange(N+1):
				v = 2*float(j)/N - 1.0
				v = math.tan(v*math.pi*0.25)
				for k in xrange(N+1):
					w = 2*float(k)/N -1.0
					w = math.tan(w*math.pi*0.25)
					idx = [i,j,k]
					idx = idx[:ax] + [N] + idx[ax:]
					idx = tuple(idx)
					if idx in namemap:
						continue
					else:
						namemap[idx] = idx
						flipidx = tuple([N-x for x in idx])
						namemap[flipidx] = idx
						q = [u,v,w]
						q = q[:ax] + [1.0] + q[ax:]
						#print idx,"->",q
						q = vectorops.unit(q)
						R = so3.mul(so3.inv(R0),so3.from_quaternion(q))
						G.add_node(idx,params=R)
	for ax in range(4):
		for i in xrange(N):
			for j in xrange(N):
				for k in xrange(N):
					baseidx = [i,j,k]
					idx = baseidx[:ax] + [N] + baseidx[ax:]
					idx = tuple(idx)
					for n in range(3):
						nidx = baseidx[:]
						nidx[n] += 1
						nidx = nidx[:ax] + [N] + nidx[ax:]
						nidx = tuple(nidx)
						G.add_edge(namemap[idx],namemap[nidx])
	return G

def so3_staggered_grid(N):
	"""Returns a nx.Graph describing a staggered grid on SO3 with resolution N. 
	The resulting graph has ~8N^3 nodes

	The node attribute 'params' gives the so3 element.
	"""
	import itertools
	assert N >= 1
	G = nx.Graph()
	R0 = so3.from_quaternion(vectorops.unit((1,1,1,1)))
	#R0 = so3.identity()
	namemap = dict()
	for ax in range(4):
		for i in xrange(N+1):
			u = 2*float(i)/N - 1.0
			u = math.tan(u*math.pi*0.25)
			u2 = 2*float(i + 0.5)/N - 1.0
			u2 = math.tan(u2*math.pi*0.25)
			for j in xrange(N+1):
				v = 2*float(j)/N - 1.0
				v = math.tan(v*math.pi*0.25)
				v2 = 2*float(j+ 0.5)/N - 1.0
				v2 = math.tan(v2*math.pi*0.25)
				for k in xrange(N+1):
					w = 2*float(k)/N -1.0
					w = math.tan(w*math.pi*0.25)
					w2 = 2*float(k + 0.5)/N -1.0
					w2 = math.tan(w2*math.pi*0.25)
				
					idx = [i,j,k]
					idx = idx[:ax] + [N] + idx[ax:]
					idx = tuple(idx)
					if idx in namemap:
						pass
					else:
						namemap[idx] = idx
						flipidx = tuple([N-x for x in idx])
						namemap[flipidx] = idx
						q = [u,v,w]
						q = q[:ax] + [1.0] + q[ax:]
						q = vectorops.unit(q)
						R = so3.mul(so3.inv(R0),so3.from_quaternion(q))
						G.add_node(idx,params=R)

					if i < N and j < N and k < N:
						#add staggered point
						sidx = [i+0.5,j+0.5,k+0.5]
						sidx = sidx[:ax] + [N] + sidx[ax:]
						sidx = tuple(sidx)
						sfidx = tuple([N-x for x in sidx])
						namemap[sidx] = sidx
						namemap[sfidx] = sidx
						q = [u2,v2,w2]
						q = q[:ax] + [1.0] + q[ax:]
						q = vectorops.unit(q)
						R = so3.mul(so3.inv(R0),so3.from_quaternion(q))
						G.add_node(sidx,params=R)
	for ax in range(4):
		for i in xrange(N):
			for j in xrange(N):
				for k in xrange(N):
					baseidx = [i,j,k]
					idx = baseidx[:ax] + [N] + baseidx[ax:]
					idx = tuple(idx)
					for n in range(3):
						nidx = baseidx[:]
						nidx[n] += 1
						nidx = nidx[:ax] + [N] + nidx[ax:]
						nidx = tuple(nidx)
						G.add_edge(namemap[idx],namemap[nidx])

					#edges between staggered point and grid points
					baseidx = [i+0.5,j+0.5,k+0.5]
					idx = baseidx[:ax] + [N] + baseidx[ax:]
					idx = tuple(idx)
					for ofs in itertools.product(*[(-0.5,0.5)]*3):
						nidx = vectorops.add(baseidx,ofs)
						nidx = [int(x) for x in nidx]
						nidx = nidx[:ax] + [N] + nidx[ax:]
						nidx = tuple(nidx)
						#print "Stagger-grid edge",idx,nidx
						G.add_edge(namemap[idx],namemap[nidx])

					#edges between staggered points -- if it goes beyond face, 
					#go to the adjacent face
					for n in range(3):
						nidx = baseidx[:]
						nidx[n] += 1
						if nidx[n] > N:
							#swap face
							nidx[n] = N
							nidx = nidx[:ax] + [N-0.5] + nidx[ax:]
							nidx = tuple(nidx)
							#if not G.has_edge(namemap[idx],namemap[nidx]):
							#	print "Stagger-stagger edge",idx,nidx,"swap face"
						else:
							nidx = nidx[:ax] + [N] + nidx[ax:]
							nidx = tuple(nidx)
							#if not G.has_edge(namemap[idx],namemap[nidx]):
							#	print "Stagger-stagger edge",idx,nidx
						G.add_edge(namemap[idx],namemap[nidx])
	return G


def so3_grid_test(N=5,staggered=True):
	from klampt import vis
	from klampt.model import trajectory
	if staggered:
		G = so3_staggered_grid(N)
	else:
		G = so3_grid(N)
	vispt = [1,0,0]
	for n in G.nodes():
		R = G.node[n]['params']
		trans = so3.apply(R,vispt)
		#trans = so3.axis_angle(R)[0]
		#trans = vectorops.unit(so3.quaternion(R)[:3])
		vis.add(str(n),[R,trans])
		vis.hideLabel(str(n))
	#draw edges?
	minDist = float('inf')
	maxDist = 0.0
	for i,j in G.edges():
		Ri = G.node[i]['params']
		Rj = G.node[j]['params']
		tmax = 9
		times = range(tmax+1)
		milestones = []
		for t in times:
			u = float(t)/tmax
			trans = so3.apply(so3.interpolate(Ri,Rj,u),vispt)
			milestones.append(trans)
		vis.add(str(i)+'-'+str(j),trajectory.Trajectory(times,milestones))
		vis.hideLabel(str(i)+'-'+str(j))
		dist = so3.distance(Ri,Rj)
		if dist > maxDist:
			maxDist = dist
			print "Max dispersion at",i,j,":",dist
		if dist < minDist:
			minDist = dist
	print "Number of points:",G.number_of_nodes()
	print "Min/max dispersion:",minDist,maxDist
	vis.run()
