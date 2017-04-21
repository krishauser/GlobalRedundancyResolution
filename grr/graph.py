from klampt import *
from klampt.model import ik,trajectory,collide
from klampt.math import vectorops,so3,se3
from ikdb.ikproblem import IKProblem,IKSolverParams
from utils.so3_grid import *
import networkx as nx
import scipy.spatial
import itertools
import random
import math
from utils.utils import *
from ikdb.utils import mkdir_p
from utils.nearestneighbors import *
from utils.metric import *
from collections import deque
import os

class MeshTopology:
	def __init__(self,verts=None,faces=None):
		self.verts = verts
		self.faces = faces if faces is not None else []
		self.incidentVerts = dict()
		self.incidentFaces = dict()
		self.computeTopology()

	def computeTopology(self):
		self.incidentFaces = dict()
		self.incidentVerts = dict()
		for i,f in enumerate(self.faces):
			for v in f:
				if v not in self.incidentFaces:
					self.incidentFaces[v] = []
					self.incidentVerts[v] = []
				self.incidentFaces[v].append(i)
				for w in f:
					if w != v:
						self.incidentVerts[v].append(w)

	def addFace(self,face):
		index = len(self.faces)
		self.faces.append(face)
		for v in face:
			if v not in self.incidentFaces:
				self.incidentFaces[v] = []
				self.incidentVerts[v] = []
			self.incidentFaces[v].append(index)
			for w in face:
				if w != v:
					self.incidentVerts[v].append(w)

	def mergeVertex(self,v,w,fairing=False):
		"""Merges vertex v to w.  If fairing = True, tries to reshape
		the mesh so that the vertices on edges are closer (only works in triangle meshes)"""
		assert self.incidentFaces.get(v,None) is not None and len(self.incidentFaces[v]) >= 1
		if len(self.faces[self.incidentFaces[v][0]])==2 and False:
			if1 = self.incidentFaces[v][0]
			if2 = self.incidentFaces[v][1]
			f1 = self.faces[if1]
			f2 = self.faces[if2]
			v1 = (f1[0] if f1[1] == v else f1[1])
			v2 = (f2[0] if f2[1] == v else f2[1])
			self.faces[if2] = None
			self.faces[if1] = [v1,v2]
			self.incidentFaces[v2].remove(if2)
			self.incidentFaces[v2].append(if1)
			self.incidentVerts[v1].remove(v)
			self.incidentVerts[v2].remove(v)
			self.incidentVerts[v1].append(v2)
			self.incidentVerts[v2].append(v1)
			del self.incidentVerts[v]
			del self.incidentFaces[v]
		else:
			#all faces with v but not w get v replaced by w
			#all faces with both v and w get removed
			#print "...Modifying",len(incidentFaces[v]),"faces"
			modified = []
			for f in self.incidentFaces[v]:
				if w in self.faces[f]:
					#print "  Deleting face",f,"=",faces[f]
					for u in self.faces[f]:
						if u != v:
							self.incidentFaces[u].remove(f)
					self.faces[f] = None
				else:
					#replace v with w
					#print "  Face",f,"=",faces[f],
					self.faces[f] = tuple([(w if u==v else u) for u in self.faces[f]])
					modified.append(f)
					assert v not in self.faces[f]
					#print "becomes",faces[f]
					for u in self.faces[f]:
						if u != w:
							self.incidentVerts[u].append(w)
							self.incidentVerts[w].append(u)
					self.incidentFaces[w].append(f)
			for u in self.incidentVerts[v]:
				assert v in self.incidentVerts[u]
				self.incidentVerts[u].remove(v)
			#sanity check
			for f in modified:
				assert v not in self.faces[f]
				for u in self.faces[f]:
					assert f in self.incidentFaces[u]
					for w in self.faces[f]:
						if u != w:
							assert u in self.incidentVerts[w]
			#do edge fairing amongst incident faces.
			#1. detect any pairs of faces with shared edges
			while fairing:
				assert self.verts is not None
				change = False
				for i,f1 in enumerate(modified):
					for (j,f2) in enumerate(modified[:i]):
						#print f1,f2
						shared = set(self.faces[f1]) & set(self.faces[f2])
						if len(shared) == 2:
							#shared, determine whether edge flip fixes
							unshared1 = [u for u in self.faces[f1] if u not in shared]
							unshared2 = [u for u in self.faces[f2] if u not in shared]
							unshared = unshared1+unshared2
							assert len(unshared1)==1 and len(unshared)==2
							shared = list(shared)
							if vectorops.distance(self.verts[unshared[0]],self.verts[unshared[1]]) < vectorops.distance(self.verts[shared[0]],self.verts[shared[1]]):
								self.edgeFlip(f1,f2,shared[0],shared[1])
								change = True
								break
					if change:
						break
				if not change:
					break
				#sanity check
				for f in modified:
					assert v not in self.faces[f]
					for u in self.faces[f]:
						assert f in self.incidentFaces[u]
						for w in self.faces[f]:
							if u != w:
								assert u in self.incidentVerts[w]
			del self.incidentFaces[v]
			del self.incidentVerts[v]

		def edgeFlip(self,f1,f2,v1,v2):
			"""For a triangle mesh and an edge (v1,v2) common to faces f1 and f2, flips the edge."""
			assert len(self.faces[f1]) == 3
			assert len(self.faces[f2]) == 3
			assert v1 in self.faces[f1]
			assert v2 in self.faces[f1]
			assert v1 in self.faces[f2]
			assert v2 in self.faces[f2]
			#do a flip
			unshared1 = [u for u in self.faces[f1] if u not in shared]
			unshared2 = [u for u in self.faces[f2] if u not in shared]
			unshared = unshared1+unshared2
			assert len(unshared1)==1 and len(unshared)==2
			self.faces[f1] = tuple(sorted([v1,unshared[0],unshared[1]]))
			self.faces[f2] = tuple(sorted([v2,unshared[0],unshared[1]]))
			self.incidentFaces[v2].remove(f1)
			self.incidentFaces[v1].remove(f2)
			self.incidentFaces[unshared[1]].append(f1)
			self.incidentFaces[unshared[0]].append(f2)
			self.incidentVerts[v1].remove(v2)
			self.incidentVerts[v2].remove(v1)
			self.incidentVerts[unshared[0]].append(unshared[1])
			self.incidentVerts[unshared[1]].append(unshared[0])
			#sanity check
			assert unshared[0] in self.incidentVerts[unshared[1]]
			assert unshared[1] in self.incidentVerts[unshared[0]]
			assert v1 in self.incidentVerts[unshared[0]]
			assert unshared[0] in self.incidentVerts[v1]
			assert v1 in self.incidentVerts[unshared[1]]
			assert unshared[1] in self.incidentVerts[v1]
			assert v2 in self.incidentVerts[unshared[0]]
			assert unshared[0] in self.incidentVerts[v2]
			assert v2 in self.incidentVerts[unshared[1]]
			assert unshared[1] in self.incidentVerts[v2]
			for u in self.faces[f1]:
				assert f1 in self.incidentFaces[u]
				for w in self.faces[f1]:
					if u != w:
						assert u in self.incidentVerts[w]
			for u in self.faces[f2]:
				assert f2 in self.incidentFaces[u]
				for w in self.faces[f2]:
					if u != w:
						assert u in self.incidentVerts[w]

class RedundancyResolutionGraph:
	"""A graph storing a map from workspace points to configurations.  This also stores information
	about the redundancy resolution problem.

	The key member is Gw, which stores a roadmap in workspace along with its
	associated configurations.  Nodes in Gw store 'params' the workspace parameter and
	'config' the configuration (may be None or empty, if no node is assigned).  Edges in Gw
	may store a flag 'connected' which may be True if there is a valid interpolating path
	between configurations associated with the ends of the edge.

	Other members include
	- world: if not None, contains the robot and the obstacles used for collision avoidance
	- robot: the robot used for the redundancy resolution
	- domain: the workspace domain under consideration
	- ikTemplate: stores the IK problem under consideration.
	- ikSolverParams: defines parameters for the IK solver.
	- nnw: (internally used) nearest neighbor structure for seeking workspace graph nodes
	- wRadius: (internally used) max radius of the workspace graph
	"""
	def __init__(self,world_or_robot):
		if isinstance(world_or_robot,RobotModel):
			self.world = None
			self.robot = world_or_robot
		else:
			self.world = world_or_robot
			self.robot = world_or_robot.robot(0)
		self.spinJoints = False
		self.domain = None
		self.Gw = nx.Graph()
		self.ikTemplate = IKProblem()
		self.ikTemplate.setJointLimits()
		self.ikSolverParams = IKSolverParams()
		self.ikSolverParams.startRandom = True
		#self.ikSolverParams.numIters = 20
		self.nnw = NearestNeighbors(L2Metric,'kdtree')
		self.wRadius = None
		self.resolvedCount = 0

	def hasResolution(self):
		return self.resolvedCount != 0

	def readJsonSetup(self,jsonObj):
		"""Reads the ik problem, the domain, and other workspace graph options from the
		given json object (examples of the format is given in the problems/ directory)"""
		pdef = jsonObj
		orientation = pdef.get('orientation','free')
		domain = pdef['domain']
		eelocal = pdef.get('eelocal',(0,0,0))
		useCollision = pdef.get('useCollision',True)
		link = pdef.get('link',None)
		links = pdef.get('links',None)
		fixedOrientation = pdef.get('fixedOrientation',None)
		if isinstance(fixedOrientation,str):
			#assume some klamp't code
			fixedOrientation = eval(fixedOrientation)
		self.domain = domain

		#setup defaults
		robot = self.robot
		if link == None:
			link = robot.link(robot.numLinks()-1)
		else:
			link = robot.link(link)
		assert link.index >= 0,"Invalid link %s specified in JSON file"%(str(link),)
		if links:
			self.ikTemplate.activeDofs = [robot.link(l).index for l in links]
			assert all(v >= 0 for v in self.ikTemplate.activeDofs),"Invalid active DOF specified in JSON file"

		if orientation == 'fixed':
			if fixedOrientation == None:
				fixedOrientation = link.getTransform()[0]

		movinglinks = []
		stilllinks = []
		movingrobotgeoms = []
		stillrobotgeoms = []
		moving = [False]*robot.numLinks()
		if self.ikTemplate.activeDofs is None:
			moving = [True]*robot.numLinks()
		else:
			subset = set(self.ikTemplate.activeDofs)
			for i in range(robot.numLinks()):
				if i in subset: moving[i] = True
				else:
					p = robot.link(i).getParent()
					if p >= 0 and moving[p]: moving[i]=True
		for i in xrange(robot.numLinks()):
			g = robot.link(i).geometry()
			if not g.empty():
				if moving[i]:
					movingrobotgeoms.append(g)
					movinglinks.append(i)
				else:
					stillrobotgeoms.append(g)
					stilllinks.append(i)
		worldgeoms = []
		if self.world:
			#we'll do collision testing with the world. extract all geometries from the world
			for i in xrange(self.world.numRigidObjects()):
				g = self.world.rigidObject(i).geometry()
				if not g.empty():
					worldgeoms.append(g)
			
		def collisionFree(q):
			robot.setConfig(q)
			#if robot.selfCollides():
			#	return False
			for x in collide.self_collision_iter(movingrobotgeoms,lambda i,j:robot.selfCollisionEnabled(movinglinks[i],movinglinks[j])):
				return False
			for x in collide.group_collision_iter(movingrobotgeoms,stillrobotgeoms,lambda i,j:robot.selfCollisionEnabled(movinglinks[i],stilllinks[j])):
				return False
			for x in collide.group_collision_iter(movingrobotgeoms,worldgeoms):
				return False
			return True

		#now set up the problem
		if orientation == 'fixed':
			self.addIKConstraint(link,eelocal,orientation=fixedOrientation)
		else:
			self.addIKConstraint(link,eelocal,orientation=orientation)
			if orientation == 'variable':
				if len(domain[0])==3:
					active = [1 if (a!=b) else 0 for (a,b) in zip(*domain)]
					if sum(active)==2:
						#pick the one non-variable axis
						rmin = [0,0,0]
						rmax = [0,0,0]
						for i,a in enumerate(active):
							if a==0:
								rmin[i] = -math.pi
								rmax[i] = math.pi
						domain = (domain[0]+rmin,domain[1]+rmax)
					else:
						domain = (domain[0]+[-math.pi]*3,domain[1]+[-math.pi]*3)
				#need to change the domain
				self.domain = domain
		if useCollision:
			self.ikTemplate.feasibilityTest = collisionFree


	def save(self,folder):
		mkdir_p(folder)
		nx.write_gpickle(self.Gw,os.path.join(folder,"resolution_graph.pickle"))
		return

	def load(self,folder):
		self.Gw = nx.read_gpickle(os.path.join(folder,"resolution_graph.pickle"))
		#setup domain
		bmin = self.Gw.node[0]['params'][:]
		bmax = self.Gw.node[0]['params'][:]
		self.resolvedCount = 0
		for i,d in self.Gw.nodes_iter(data=True):
			if d.get('config',None) is not None:
				self.resolvedCount += 1
			x = d['params']
			for j,v in enumerate(x):
				bmin[j] = min(bmin[j],v)
				bmax[j] = max(bmax[j],v)
		self.domain = (bmin,bmax)
		self.buildNearestNeighbors()

	def buildNearestNeighbors(self):
		#setup nearest neighbors
		assert self.domain != None,"called buildNearestNeighbors without a domain"
		bmin,bmax = self.domain
		if len(bmax) == 6:
			self.nnw = NearestNeighbors(workspace_distance,'bruteforce')
		else:
			self.nnw = NearestNeighbors(L2Metric,'kdtree')
		for n,d in self.Gw.nodes_iter(data=True):
			self.nnw.add(d['params'],n)

	def printStats(self):
		print "Roadmap has",self.Gw.number_of_nodes(),"workspace nodes and",self.Gw.number_of_edges(),"edges"
		numResolved = 0
		for (i,d) in self.Gw.nodes_iter(data=True):
			if d.get('config',None) is not None:
				numResolved += 1
		numConnected = 0
		for (i,j,d) in self.Gw.edges_iter(data=True):
			if d.get('connected',False):
				numConnected += 1
		print numResolved,"nodes are resolved and",numConnected,"edge are resolved"
		sumdistances = 0
		sumwsdistances = 0
		for i,j,d in self.Gw.edges_iter(data=True):
			if d.get('connected',False):
				sumdistances += self.robot.distance(self.Gw.node[i]['config'],self.Gw.node[j]['config'])
				sumwsdistances += workspace_distance(self.Gw.node[i]['params'],self.Gw.node[j]['params'])
		if numConnected != 0:
			print "  Average configuration space distance:",sumdistances / numConnected
			print "  Average configuration space distance / workspace distance:",sumdistances / sumwsdistances

	def verify(self):
		numInvalid = 0
		for iw,jw in self.Gw.edges():
			if self.isResolvedEdge(iw,jw):
				if not self.validEdge(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
					print "Warning, edge",iw,jw,"is marked as connected, but it is actually not valid"
					numInvalid += 1
					del self.Gw.edge[iw][jw]['connected']
		print "Verification found",numInvalid,"/",self.Gw.number_of_edges()+numInvalid,"edges infeasible"

	def clear(self):
		"""Reinitializes workspace graph."""
		self.Gw = nx.Graph()
		self.nnw = NearestNeighbors(L2Metric,'kdtree')
		self.resolvedCount = 0

	def clearResolution(self):
		"""Reinitializes resolution, keeping workspace graph the same"""
		for i,d in self.Gw.nodes_iter(data=True):
			if 'config' in d:
				del d['config']
		for i,j,d in self.Gw.edges_iter(data=True):
			if 'connected' in d:
				del d['connected']
		self.resolvedCount = 0

	def addIKConstraint(self,link,eelocal,orientation='free'):
		"""Orientation can be free, variable, or a so3 element"""
		if orientation == 'free':
			self.ikTemplate.addObjective(ik.objective(link,local=eelocal,world=[0,0,0]))
		elif orientation == 'variable':
			T = se3.identity()
			obj = ik.objective(link,R=T[0],t=T[1])
			obj.setFixedPosConstraint(eelocal,[0,0,0])
			self.ikTemplate.addObjective(obj)
		else:
			assert isinstance(orientation,(list,tuple))
			#assume it's a fixed orientation matrix
			obj = ik.objective(link,R=orientation,t=[0,0,0])
			obj.setFixedPosConstraint(eelocal,[0,0,0])
			self.ikTemplate.addObjective(obj)

	def setIKProblem(self,x):
		"""Given a workspace setting x, sets the self.ikTemplate problem so that
		it corresponds with the workspace parameters x"""
		n = 0
		for obj in self.ikTemplate.objectives:
			np = obj.numPosDims()
			nr = obj.numRotDims()
			assert n + np <= len(x)
			assert np == 3 or np == 0,"Can only handle fixed or free position constraints"
			assert nr == 3 or nr == 0,"Can only handle fixed or free rotation constraints"
			if np == 3:
				local,world = obj.getPosition()
				obj.setFixedPosConstraint(local,x[n:n+3])
			n += np
			if nr == 3 and len(x) >= 3 + n:
				assert n + nr <= len(x)
				R = obj.getRotation()
				R = so3.from_moment(x[n:n+3])
				obj.setFixedRotConstraint(R)
				n += nr
		assert n == len(x),"Workspace parameters must be exactly the right number of variables in template IK constraint"

	def getIKParameters(self,q):
		"""Given a configuration q, gets the workspace parameters that matches q exactly."""
		x = []
		qOrig = self.robot.getConfig()
		self.robot.setConfig(q)
		for obj in self.ikTemplate.objectives:
			link = self.robot.link(obj.link())
			np = obj.numPosDims()
			nr = obj.numRotDims()
			assert np == 3 or np == 0,"Can only handle fixed or free position constraints"
			assert nr == 3 or nr == 0,"Can only handle fixed or free rotation constraints"
			if np == 3:
				local,world = obj.getPosition()
				x += link.getWorldPosition(local)
			#TODO: handle rotations
		self.robot.setConfig(qOrig)
		return x

	def solve(self,qinit,x):
		"""Solves the IK problem at x with the starting config qinit.  Returns None if
		failed, or the solution if successful."""
		self.setIKProblem(x)
		self.robot.setConfig(qinit)
		self.ikSolverParams.startRandom = False
		return self.ikTemplate.solve(self.robot,self.ikSolverParams)

	def addWorkspaceNode(self,x):
		"""Adds a node to the workspace graph with params x.  Returns the index"""
		iw = self.Gw.number_of_nodes()
		self.Gw.add_node(iw,params=x)
		self.nnw.add(x,iw)
		return iw

	def addWorkspaceEdge(self,i,j):
		"""Adds an edge to the workspace graph"""
		if i > j:
			self.addWorkspaceEdge(j,i)
			return
		assert i < self.Gw.number_of_nodes()
		assert j < self.Gw.number_of_nodes()
		self.Gw.add_edge(i,j)

	def connectWorkspaceNeighbors(self,iw,radius=None,k=None):
		"""Connects a workspace node with its neighbors, either in a radius or the k-nearest
		neighbors.
		"""
		x = self.Gw.node[iw]['params']
		if radius != None:
			knn = self.nnw.neighbors(x,radius)
			if len(x) > 3:
				#rotational, check connections on opposite side of unit-sphere
				rx = x[3:]
				assert len(rx) <= 3
				if vectorops.norm(rx) > 0.01:
					rxnew = vectorops.madd(rx,vectorops.unit(rx),-2*math.pi)
					xnew = x[:3]+rxnew
					knn += self.nnw.neighbors(xnew,radius)
		else:
			assert k != None,"Need either k or radius to be specified"
			knn = self.nnw.knearest(x,k)
			if len(x) > 3:
				#rotational, check connections on opposite side of unit-sphere
				rx = x[3:]
				assert len(rx) <= 3
				if vectorops.norm(rx) > 0.01:
					rxnew = vectorops.madd(rx,vectorops.unit(rx),-2*math.pi)
					xnew = x[:3]+rxnew
					knn += self.nnw.knearest(xnew,k)
				#now subselect only the k nearest in terms of workspace distance
				sortlist = [(workspace_distance(pt,x),pt,n) for pt,n in knn]
				knn = [(pt,n) for (d,pt,n) in sorted(sortlist)]
				knn = knn[:k]
		for pt,n in knn:
			if n == iw:
				continue
			self.addWorkspaceEdge(iw,n)
	
	def removeWorkspaceEdge(self,iw,jw):
		assert iw != jw
		if iw > jw:
			self.removeWorkspaceEdge(jw,iw)
			return
		assert iw < self.Gw.number_of_nodes()
		assert jw < self.Gw.number_of_nodes()
		self.Gw.remove_edge(iw,jw)

	def setConfig(self,i,q):
		assert i < self.Gw.number_of_nodes()
		if self.Gw.node[i].get('config',None) is not None:
			self.resolvedCount -= 1
		self.Gw.node[i]['config'] = q
		if q is not None:
			self.resolvedCount += 1
		else:
			for j in self.Gw.neighbors(i):
				if self.Gw.edge[i][j].get('connected',False):
					del self.Gw.edge[i][j]['connected']

	def markConnected(self,iw,jw,connected=True):
		"""Adds a configuration space edge corresponding to workspace/configuration space node index (iw,iq) to
		index (jw,jq).
		"""
		if iw > jw:
			self.markConnected(jw,iw,connected)
			return
		assert iw < self.Gw.number_of_nodes()
		assert jw < self.Gw.number_of_nodes()
		if connected:
			assert self.Gw.node[iw].get('config',None) is not None
			assert self.Gw.node[jw].get('config',None) is not None
			#TEMP: sanity check
			#if not self.validEdge(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
			#	print "Warning, edge",iw,jw,"is being marked connected, but the edge is not valid"
		self.Gw.edge[iw][jw]['connected'] = connected

	def isResolvedNode(self,iw):
		return self.Gw.node[iw].get('config',None) is not None

	def isResolvedEdge(self,iw,jw):
		return self.Gw.edge[iw][jw].get('connected',False) == True

	def ikError(self,q,x):
		"""Returns the ik error (norm of residual) at config q for workspace parameter x"""
		self.setIKProblem(x)
		s = ik.solver(self.ikTemplate.objectives)
		self.robot.setConfig(q)
		r = s.getResidual()
		return vectorops.norm(r)

	def validEdgeSimple(self,qa,qb,wa,wb,epsilon=5e-2):
		ea = self.ikError(qb,wa)
		eb = self.ikError(qa,wb)
		if self.wRadius == None:
			self.autoSetWorkspaceRadius()
		emax = self.wRadius*0.5
		#midpoint fast reject test
		x = self.robot.interpolate(qa,qb,0.5)
		if self.ikError(x,workspace_interpolate(wa,wb,0.5)) > emax:
			#print "Reject deviation from line too large at midpoint"
			return False
		#if self.ikError(x,wa) > ea or self.ikError(x,wb) > eb:
		#	return False
		eaold = 0
		ebold = eb
		Ndivs = int(math.ceil(self.robot.distance(qa,qb)/epsilon))
		path = []
		for i in range(Ndivs):
			u = float(i+1)/(Ndivs+1)
			x = self.robot.interpolate(qa,qb,u)
			eu = self.ikError(x,workspace_interpolate(wa,wb,u))
			eau = self.ikError(x,wa)
			ebu = self.ikError(x,wb)
			#if eau < eaold:
			#	#print "Reject monotonic increase from a at",u
			#	return False
			#if ebu > ebold:
			#	#print "Reject monotonic decrease from b at",u,",",ebold,"->",ebu
			#	return False
			if eu > emax:
				#print "Reject deviation from line too large at",u
				return False
			#if eau > ea or ebu > eb:
			#	return False
			eaold = eau
			ebold = ebu
			path.append(x)
		if self.ikTemplate.feasibilityTest != None:
			for q in path:
				if not self.ikTemplate.feasibilityTest(q):
					return False
		return True

	def validEdgeLinear(self,qa,qb,wa,wb,epsilon=1e-3):
		#ea = self.ikError(qb,wa)
		#eb = self.ikError(qa,wb)
		qprev = qa
		#TODO determine a suitable lipschitz constant for singularity avoidance
		#right now this assumes that each DOF contributes about 1 unit to the lipschitz constant
		nlinks = (self.robot.numLinks() if self.ikTemplate.activeDofs is None else len(self.ikTemplate.activeDofs))
		epsilon *= nlinks
		Ndivs = int(math.ceil(self.robot.distance(qa,qb)/epsilon))
		discontinuityThreshold = epsilon*10

		#this does a bisection technique and should be faster for infeasible edges
		queue = deque()
		queue.append((0,Ndivs+1,qa,qb))
		while len(queue)>0:
			ia,ib,q0,q1 = queue.popleft()
			if ib == ia+1:
				d = self.robot.distance(q0,q1)
				if d > discontinuityThreshold:
					#print "Discontinuity threshold exceeded",d,discontinuityThreshold
					return False
				continue
			im = (ia+ib)/2
			u = float(im)/(Ndivs+1)
			x = self.robot.interpolate(qa,qb,u)
			qm = self.solve(x,workspace_interpolate(wa,wb,u))
			if qm == None:
				#print "Unable to solve",x
				return False
			d = self.robot.distance(q0,qm)
			if d > discontinuityThreshold*(ib-ia):
				#print "Discontinuity threshold exceeded mid-segment",d,discontinuityThreshold*(ib-ia)
				return False
			d = self.robot.distance(qm,q0)
			if d > discontinuityThreshold*(ib-ia):
				#print "Discontinuity threshold exceeded mid-segment",d,discontinuityThreshold*(ib-ia)
				return False
			if self.ikTemplate.feasibilityTest is not None:
				if not self.ikTemplate.feasibilityTest(qm):
					#print "Midpoint is in collision"
					return False
			queue.append((ia,im,q0,qm))
			queue.append((im,ib,qm,q1))
		"""
		for i in range(Ndivs):
			u = float(i+1)/(Ndivs+1)
			x = self.robot.interpolate(qa,qb,u)
			q = self.solve(x,workspace_interpolate(wa,wb,u))
			if q == None:
				#print "Unable to solve",x
				return False
			path.append(q)
			d = self.robot.distance(q,qprev)
			if d > discontinuityThreshold:
				#print "Discontinuity threshold exceeded",d,discontinuityThreshold
				return False
			qprev = q
		#check continuity at final configuration
		d = self.robot.distance(qb,qprev)
		if d > discontinuityThreshold:
			#print "Discontinuity threshold exceeded",d,discontinuityThreshold
			return False
		#test feasibility of path
		if self.ikTemplate.feasibilityTest != None:
			for q in path:
				if not self.ikTemplate.feasibilityTest(q):
					return False
		"""
		return True


	def interpolateEdgeLinear(self,qa,qb,wa,wb,u):
		x = self.robot.interpolate(qa,qb,u)
		q = self.solve(x,workspace_interpolate(wa,wb,u))
		if q is None:
			return self.robot.interpolate(qa,qb,u)
		return q

	def validEdgeBisection(self,qa,qb,wa,wb,epsilon=5e-2,c=0.9):
		d0 = self.robot.distance(qa,qb)
		if d0 <= epsilon:
			return True
		wm = workspace_interpolate(wa,wb,0.5)
		qm = self.robot.interpolate(qa,qb,0.5)
		q = self.solve(qm,wm)
		if q is None:
			return False
		d1 = self.robot.distance(qa,q)
		d2 = self.robot.distance(q,qb)
		if max(d1,d2) > c*d0: return False
		if d1 > epsilon and not self.validEdgeBisection(qa,q,wa,wm,epsilon,c):
			return False
		if d2 > epsilon and not self.validEdgeBisection(q,qb,wm,wb,epsilon,c):
			return False
		return True

	def interpolateEdgeBisection(self,qa,qb,wa,wb,u,epsilon=5e-2,c=0.9):
		d0 = self.robot.distance(qa,qb)
		if d0 <= epsilon:
			return self.robot.interpolate(qa,qb,u)
		wm = workspace_interpolate(wa,wb,0.5)
		qm = self.robot.interpolate(qa,qb,0.5)
		q = self.solve(qm,wm)
		if q is None:
			return self.robot.interpolate(qa,qb,u)
		d1 = self.robot.distance(qa,q)
		d2 = self.robot.distance(q,qb)
		if max(d1,d2) > c*d0: return self.robot.interpolate(qa,qb,u)
		if u < 0.5:
			return self.interpolateEdgeBisection(qa,q,wa,wm,u*2,epsilon,c)
		else:
			return self.interpolateEdgeBisection(q,qb,wm,wb,(u-0.5)*2,epsilon,c)

	def validEdge(self,qa,qb,wa,wb):
		#return self.validEdgeBisection(qa,qb,wa,wb)
		return self.validEdgeLinear(qa,qb,wa,wb)
		#return self.validEdgeSimple(qa,qb,wa,wb)


	def sampleWorkspace(self,N,k='auto',method='random'):
		"""Samples N workspace points in the domain [bmin,bmax] using some method and connects them
		together.

		method can be one of
		- random: samples at random
		- grid: uniform grid
		- staggered_grid: uniform grid with staggering (which makes nicer shaped simplices)

		k = 'auto' automatically connects the appropriate neighbors in the workspace graph.
		Otherwise, k-nearest neighbors is used to connect the graph.

		Accepts len(bmin) = 2, 3, or 6. If len(bmin) is 2 or 3, the workspace is Cartesian space.
		If len(bmin) == 6, the workspace is SE3, and bmin[3:6] and bmax[3:6] are assumed to be [-pi,pi].
		All rotations in SO(3) are sampled.  (The distance metric used for k-nearest neighbors is...?)
		"""
		assert self.domain is not None,"Need to set domain before calling sampleWorkspace"
		(bmin,bmax) = self.domain 
		print "Sampling",N,"workspace points from",bmin,"->",bmax,"with method",method
		dims = len([1 for (a,b) in zip(bmin,bmax) if a!=b])
		if method == 'random':
			for i in xrange(N):
				x = [random.uniform(a,b) for (a,b) in zip(bmin,bmax)]
				if len(bmin) == 6:
					#sample rotation uniformly
					x[3:6] = so3.sample()
				self.addWorkspaceNode(x)
			if k == 'auto':
				k = int((1+1.0/dims)*math.e*math.log(N))
				#k = 4*dims
			#k-nearest neighbors
			c = ProgressUpdater(self.Gw.number_of_nodes(),(1 if len(bmin) > 3 else 10))
			for i,d in self.Gw.nodes_iter(data=True):
				c.update()
				self.connectWorkspaceNeighbors(i,k=k+1)
			c.done()
		elif method == "grid":
			add_so3 = False
			if len(bmin) == 6:
				add_so3 = True
				bmin = bmin[:3]
				bmax = bmax[:3]

			vol = (1 if not add_so3 else pow(math.pi,3.0))
			for (a,b) in zip(bmin,bmax):
				if b > a:
					vol *= (b-a)
			cellvol = float(vol)/N
			celldim = math.pow(cellvol,1.0/dims)
			assert celldim > 0
			divs = [max(int(math.ceil((b-a)/celldim)),1) for (a,b) in zip(bmin,bmax)]
			for cell in itertools.product(*[range(n) for n in divs]):
				params = [float(i)/float(div) for i,div in zip(cell,divs)]
				x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
				self.addWorkspaceNode(x)
			print "Grid discretization added",self.Gw.number_of_nodes(),"workspace nodes"
			if k == 'auto':
				celldim = max((b-a)/div for (a,b,div) in zip(bmin,bmax,divs))
				for i in self.Gw.nodes_iter():
					self.connectWorkspaceNeighbors(i,radius=celldim*1.05)
			else:
				for i,d in self.Gw.nodes_iter(data=True):
					self.connectWorkspaceNeighbors(i,k=k+1)
			if add_so3:
				print "Creating staggered rotation graph with N=",int(math.ceil(1.0/celldim))
				#Rgraph = so3_staggered_grid(int(math.ceil(1.0/celldim)))
				#create rotation graph
				Rgraph = so3_grid(int(math.ceil(1.0/celldim)))
				print "  ",Rgraph.number_of_nodes(),"nodes created"
				#form cross product of self.Gw and Rgraph
				CpGraph = nx.cartesian_product(self.Gw,Rgraph)
				print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
				#renumber, convert params from pair to 6D [translation,rotation moment] vectors
				self.Gw = nx.Graph()
				nodemap = dict()
				for n,d in CpGraph.nodes(data=True):
					t = d['params'][0]
					R = d['params'][1]
					m = so3.moment(R)
					nodemap[n] = self.Gw.number_of_nodes()
					self.Gw.add_node(self.Gw.number_of_nodes(),params = t + m)
				for i,j in CpGraph.edges():
					self.Gw.add_edge(nodemap[i],nodemap[j])
				#redo nearest neighbors structure
				self.buildNearestNeighbors()
				
		elif method == "staggered_grid":
			add_so3 = False
			if len(bmin) == 6:
				add_so3 = True
				bmin = bmin[:3]
				bmax = bmax[:3]
			
			vol = 1 if not add_so3 else pow(math.pi,3)
			for (a,b) in zip(bmin,bmax):
				if b > a:
					vol *= (b-a)
			cellvol = 2*float(vol)/N
			celldim = math.pow(cellvol,1.0/dims)
			assert celldim > 0
			divs = [max(int(math.ceil((b-a)/celldim)),1) for (a,b) in zip(bmin,bmax)]
			active = [(1 if b > a else 0) for (a,b) in zip(bmin,bmax)]
			centers = []
			for cell in itertools.product(*[range(n) for n in divs]):
				params = [float(i)/float(div) for i,div in zip(cell,divs)]
				x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
				self.addWorkspaceNode(x)
				centers.append(False)
				if all(i+a < div for (i,div,a) in zip(cell,divs,active)):
					params = [(float(i)+0.5*a)/float(div) for i,div,a in zip(cell,divs,active)]
					x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
					self.addWorkspaceNode(x)
					centers.append(True)
			print "Grid discretization added",self.Gw.number_of_nodes(),"workspace nodes"
			if k == 'auto':
				print [(b-a)/div for (a,b,div) in zip(bmin,bmax,divs)]
				celldim = max((b-a)/div for (a,b,div) in zip(bmin,bmax,divs))
				print "Connecting each workspace grid node within radius",celldim*math.sqrt(dims)*1.01
				for i in self.Gw.nodes_iter():
					if centers[i]:
						if sum(active)==2:
							self.connectWorkspaceNeighbors(i,radius=celldim*math.sqrt(dims)*0.5*1.01)
						else:
							self.connectWorkspaceNeighbors(i,radius=celldim*1.01)
					else:
						self.connectWorkspaceNeighbors(i,radius=celldim*1.01)
			else:
				for i,d in self.Gw.nodes_iter(data=True):
					self.connectWorkspaceNeighbors(i,k=k+1)

			if add_so3:
				#create rotation graph
				print "Translation graph has",self.Gw.number_of_nodes(),"nodes"
				print "Creating staggered rotation graph with N=",int(math.ceil(1.0/celldim))
				Rgraph = so3_staggered_grid(int(math.ceil(1.0/celldim)))
				print "  ",Rgraph.number_of_nodes(),"nodes created"
				#form cross product of self.Gw and Rgraph
				CpGraph = nx.cartesian_product(self.Gw,Rgraph)
				print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
				#renumber, convert params from pair to 6D [translation,rotation moment] vectors
				self.Gw = nx.Graph()
				nodemap = dict()
				for n,d in CpGraph.nodes(data=True):
					t = d['params'][0]
					R = d['params'][1]
					m = so3.moment(R)
					nodemap[n] = self.Gw.number_of_nodes()
					self.Gw.add_node(self.Gw.number_of_nodes(),params = t + m)
				for i,j in CpGraph.edges():
					self.Gw.add_edge(nodemap[i],nodemap[j])
					
				#redo nearest neighbors structure
				self.buildNearestNeighbors()
		else:
			raise ValueError("Unknown method "+method)
			

	def autoSetWorkspaceRadius(self):
		self.wRadius = 0
		for i,d in self.Gw.nodes_iter(data=True):
			for j in self.Gw.neighbors(i):
				dist = workspace_distance(d['params'],self.Gw.node[j]['params'])
				self.wRadius = max(dist,self.wRadius)
		self.wRadius *= 1.01

	def getResolutionGraph(self):
		"""Returns the graph of resolved workspace nodes and edges, along with configurations"""
		Gres = nx.Graph()
		for i,d in self.Gw.nodes_iter(data=True):
			if d.get('config',None) is not None:
				Gres.add_node(i,**d)
		for i,j,d in self.Gw.edges_iter(data=True):
			if d.get('connected',False):
				Gres.add_edge(i,j)
		return Gres

	def interpolate(self,x,r='auto',multi=False):
		"""Uses inverse distance weighted interpolation at the workspace nodes to
		derive a reasonable guess for the IK solution at workspace setting x.

		If multi=True, returns a list of possible configurations.
		Otherwise, returns a list containing a single "best" configuration.
		"""
		q0 = self.robot.getConfig()
		self.setIKProblem(x)
		s = ik.solver(self.ikTemplate.objectives)
		if not self.hasResolution():
			return [q0]
		if r == 'auto':
			if self.wRadius == None:
				self.autoSetWorkspaceRadius()
			r = self.wRadius
			#if math.log(len(self.Gw))/len(self.Gw)*20*math.e
		nw = self.nnw.neighbors(x,r)
		if len(nw) == 0:
			res = self.nnw.nearest(x)
			if res is None:
				print "Weird: can't evaluate nearest neighbors at",x,"?"
				return [q0]
			xn,n = res
			nw = [(xn,n)]
		if len(nw) == 1:
			#radius too small, look at neighbors
			nn = self.Gw.neighbors(nw[0][1])
			nw += [(self.Gw.node[n]['params'],n) for n in nn]
			r = max(workspace_distance(x,self.Gw.node[n]['params']) for (xw,n) in nw)
		#print "Point distances", [workspace_distance(self.Gw.node[n[1]]['params'],x) for n in nw]
		#print [self.Gw.node[n[1]] for n in nw]

		closest = None
		maxd = float('inf')
		for (xw,iw) in nw:
			d = workspace_distance(x,self.Gw.node[iw]['params'])
			if d < maxd:
				maxd = d
				closest = iw

		#if exactly on the workspace graph, just interpolate as normal
		edgeTol = 1e-3
		Wsubgraph = self.Gw.subgraph([n[1] for n in nw])
		for (iw,jw) in Wsubgraph.edges_iter():
			xi = self.Gw.node[iw]['params']
			xj = self.Gw.node[jw]['params']
			assert len(xi) == len(x)
			(d,u) = segment_point_distance(xi,xj,x)
			#print "Edge distance", segment_point_distance(xi,xj,x)
			if self.Gw.edge[iw][jw].get('connected',False):
				if d < maxd:
					maxd = d
					closest = self.interpolateEdgeLinear(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],xi,xj,u)
		if maxd < edgeTol and hasattr(closest,'__iter__'):
			return [closest]

		#find the largest connected subgraph amongst the config nodes that
		#map to the k workspace nodes
		Sq = nx.Graph()
		for (xw,iw) in nw:
			iq = self.Gw.node[iw].get('config',None)
			if iq is not None:
				Sq.add_node(iw)
				#print iq,'=',iw,
				for jw in self.Gw.neighbors(iw):
					if Sq.has_node(jw) and self.isResolvedEdge(iw,jw):
						Sq.add_edge(iw,jw)
		#print
		ccs = list(nx.connected_components(Sq))

		#find the cc containing the closest
		if not multi:
			if closest is not None:
				for cc in ccs:
					if closest in cc:
						ccs = [cc]
						break
		"""
		ccs = []
		used = dict((s,False) for s in subset)
		for s in subset:
			if used[s]: continue
			used[s] = True
			ccs.append([s])
			for t in Sq.neighbors(s):
				used[t] = True
				ccs[-1].append(t)
		"""
		#print ccs
		qs = []
		for cc0 in ccs:
			distances = []
			for iw in cc0:
				distances.append(workspace_distance(x,self.Gw.node[iw]['params']))
			mindist = min(distances)
			weights = {}
			if mindist < 1e-5:
				#very close to existing node, harmonic interpolation will give
				#numerical difficulty
				for (d,iw) in zip(distances,cc0):
					if d < 1e-5:
						weights[iw] = 1
			else:
				for (d,iw) in zip(distances,cc0):
					weights[iw] = ((r - d)/(r*d))**2
			if len(weights) == 1:
				iq = weights.keys()[0]
				q = self.Gw.node[iw]['config']
			else:
				sumw = 0
				q = [0.0]*self.robot.numLinks()
				for (w,iw) in sorted((-w,iw) for (iw,w) in weights.iteritems()):
					w = -w
					#print "weight",w,"config",self.Gw.node[iw]['config']
					if sumw == 0:
						q = self.Gw.node[iw]['config']
					else:
						#print q[1],
						q = self.robot.interpolate(q,self.Gw.node[iw]['config'],w/(w+sumw))
						#print "to",iq,"frac",w/(w+sumw),"res",q[1]
					sumw += w
				#print
				#print [str(iq)+":"+"%.2f"%(w/sumw,) for iq,w in weights.iteritems()]
			qs.append(q)
		robot = self.robot
		if hasattr(self,'lastqs') and len(self.lastqs)==len(qs):
			lipschitzConstant = 10
			for i,(a,b,cc) in enumerate(zip(qs,self.lastqs,ccs)):
				xstart = self.getIKParameters(a)
				if self.robot.distance(a,b) > lipschitzConstant*workspace_distance(xstart,x):
					print "Big change in configuration",i,"/",len(qs),"magnitude:",self.robot.distance(a,b),"relative to workspace distance",workspace_distance(xstart,x)
					#print zip(a,b)
					#print cc
		self.lastqs = qs
		return qs

	def eval(self,x,qOrig=None):
		"""Solves the IK problem specified by the workspace setting x"""
		seeds = self.interpolate(x,multi=True)
		self.setIKProblem(x)
		if len(seeds) == 0:
			self.ikSolverParams.startRandom = False
			if qOrig is not None:
				self.robot.setConfig(qOrig)
			res = self.ikTemplate.solve(self.robot,self.ikSolverParams)
			self.ikSolverParams.startRandom = True
			return res
		if qOrig is None:
			qOrig = self.robot.getConfig()
		best = None
		dbest = float('inf')
		self.ikSolverParams.startRandom = False
		#pick the closest configuration to qOrig
		for q in seeds:
			self.robot.setConfig(q)
			res = self.ikTemplate.solve(self.robot,self.ikSolverParams)
			if res:
				d = self.robot.distance(res,qOrig)
				if d < dbest:
					best = res
					dbest = d
		return best

	def walkPath(self,x1,x2=None):
		"""Returns a pair x,q of paths, where x=[x[0],...,x[n-1]] is a workspace path
		connecting xstart and xgoal, and q=[q[0],...,q[n-1]] is the configuration space
		path associated with the workspace path. It is guaranteed that every segment 
		of the workspace path from x[1] to x[n-2] is on the workspace graph, is resolved,
		and q[1],...,q[2] are the resolved configurations.

		The arguments x1 and x2 are interpreted as follows.  If only x1 is given, then it
		is considered to be xgoal, and the workspace parameters of the robot's current
		configuration is used as xstart.  If both are given, then x1 is xstart and x2 is xgoal.

		If for some reason the resolution is empty, then None,None is returned.
		"""
		if x2 == None:
			xsrc = self.getIKParameters(self.robot.getConfig())
			xtgt = x1
		else:
			xsrc = x1
			xtgt = x2
		if isinstance(xsrc,int):
			xsrc = self.Gw.node[xsrc]['params']
		if isinstance(xtgt,int):
			xtgt = self.Gw.node[xtgt]['params']

		#get graph of connected nodes in workspace
		Gconnected = nx.Graph()
		for i,d in self.Gw.nodes_iter(data=True):
			assert not Gconnected.has_node(i)
			if d.get('config',None) is not None:
				Gconnected.add_node(i,params=d['params'],config=d['config'])
		for i,j,d in self.Gw.edges_iter(data=True):
			if d.get('connected',False):
				#cost = linf_distance(Gconnected.node[i]['config'],Gconnected.node[j]['config'],self.spinJoints)
				cost = vectorops.distance(Gconnected.node[i]['params'],Gconnected.node[j]['params'])
				Gconnected.add_edge(i,j,cost=cost)
				Gconnected.add_edge(j,i,cost=cost)
		xpath = []
		qpath = []

		#find closest start workspace node
		res = self.nnw.nearest(xsrc)
		if res is None:
			print "No workspace node available at start point"
			return None,None
		nsrc = res[1]
		if xsrc != self.Gw.node[nsrc]['params']:
			if not Gconnected.has_node(nsrc):
				print "Start workspace node not in resolution"
				return None,None
			qpath.append(self.robot.getConfig())
			xpath.append(xsrc[:])

		#find closest target workspace node
		res = self.nnw.nearest(xtgt)
		if res is None:
			print "No workspace node available at target point"
			return None,None
		ntgt = res[1]

		#print "Src node index",nsrc
		#print "Tgt node index",ntgt
		pathindices = nx.shortest_path(Gconnected,nsrc,ntgt,'cost')
		#print "Path is",pathindices
		for idx in pathindices:
			xpath.append(Gconnected.node[idx]['params'])
			qpath.append(Gconnected.node[idx]['config'])

		if xtgt != self.Gw.node[ntgt]['params']:
			self.robot.setConfig(Gconnected.node[xtgt]['config'])
			for q in self.interpolate(xtgt):
				qsolve = self.solve(q,x)
				if qsolve != None:
					qpath.append(qsolve)
					xpath.append(xtgt)
					break
		return xpath,qpath

	def computeDiscontinuities(self,useboundary=False,orientation=None):
		"""Computes a 2D segment mesh or 3D triangulated mesh illustrating the discontinuities
		in the resolution.

		If useboundary=True, the exterior of the reachable workspace is added to the discontinuities.
		Otherwise (default) only the interior discontinuities are calculated.

		If orientation is given, this assumes the workspace is 6D translation / rotation space, and
		only the slice closest to the given orientation is visualized.

		Return value is a pair (vlist,facelist) containing a vertex list vlist and a face list facelist.
		Each entry in facelist is either a pair (2d segment) or a triple (3d triangle)
		"""
		print "Computing discontinuity map..."
		if self.wRadius == None:
			self.autoSetWorkspaceRadius()
		#active dimensions
		bmin,bmax = self.domain
		active = [i for i,(a,b) in enumerate(zip(bmin,bmax)[:3]) if b!=a]
		G = self.Gw
		if orientation is not None: print "Extracting for orientation",orientation
		if orientation != None:
			#find the closest orientation in the graph
			mbest = None
			dbest = float('inf')
			nnodes = 0
			for i,d in self.Gw.nodes(data=True):
				m = d['params'][3:]
				R = so3.from_moment(m)
				dist = so3.distance(orientation,R) 
				if dist < dbest:
					dbest = dist
					mbest = m
					nnodes = 1
				elif dist == dbest:
					nnodes += 1
			assert nnodes > 1,"Only 1 node is closest to the target orientation... is this a random 6D workspace graph?"
			#extract out the subgraph corresponding to those nodes matching the target orientation
			subgraph = []
			for i,d in self.Gw.nodes(data=True):
				m = d['params'][3:]
				if m == mbest:
					subgraph.append(i)
			print "Extracting rotation subgraph of size",len(subgraph)
			graph = nx.subgraph(self.Gw,subgraph)
		else:
			assert len(bmax) <= 3,"Can't compute discontinuity boundaries for 6D workspace graph"

		def getpoint(i):
			return G.node[i]['params'][:3]

		pts = []
		#discontinuity points: indices into pts
		epts = set()
		#maps each discontinuity pt on an edge to its indices (i,j)
		es = dict()
		#maps a discontinuity point on an edge to a range [a,b] of interpolation parameters on the edge 
		eranges = dict()
		#maps vertex points to the vertex index
		vs = dict()
		#maps any other point to the combo of vertices associated with it
		others = dict()

		vfailures = set()
		for i,d in G.nodes_iter(data=True):
			if not self.isResolvedNode(i):
				#external node
				if any(self.isResolvedNode(j) for j in G.neighbors(i)):
					if useboundary:
						vindex = len(pts)
						#epts.add(vindex)
						vs[vindex] = i
						pts.append(getpoint(i))
					else:
						#add a point so that the triangulation doesnt jump concavities
						vs[len(pts)] = i
						pts.append(getpoint(i))
						pass
				continue
			#add a vertex point
			vs[len(pts)] = i
			pts.append(getpoint(i))
		for i,j in G.edges_iter():
			if not useboundary:
				if not self.isResolvedNode(i) or not self.isResolvedNode(j): continue
			else:
				if not self.isResolvedNode(i) and not self.isResolvedNode(j):
					continue
			if self.isResolvedEdge(i,j):
				pts.append(vectorops.interpolate(getpoint(i),getpoint(j),0.5))
			else:
				vfailures.add(i)
				vfailures.add(j)
				#add midpoints of disconnections
				vindex = len(pts)
				epts.add(vindex)
				es[vindex] = (i,j)
				eranges[vindex] = (0,1)
				pts.append(vectorops.interpolate(getpoint(i),getpoint(j),0.5))

				#determine valid ranges by testing resolution
				a,b = G.node[i]['params'],G.node[j]['params']
				resolution = 1e-2
				if self.isResolvedNode(i) != self.isResolvedNode(j):
					#external discontinuity: interpolate from i toward j (or j toward i)
					flip = self.isResolvedNode(j)
					if flip:
						q0 = G.node[j]['config']
						a,b = b,a
					else:
						q0 = G.node[i]['config']
					l,u = 0.0,1.0
					while u > l+resolution:
						m = 0.5*(l+u)
						x = workspace_interpolate(a,b,m)
						if self.solve(q0,x) is not None:
							l = m
						else:
							u = m

					#pts[-1] = vectorops.interpolate(getpoint(i),getpoint(j),l)
					if flip:
						l = 1.0-l
					eranges[vindex] = (l,l)
				else:
					#internal discontinuity: determine maxima from either side
					q0 = G.node[i]['config']
					q1 = G.node[j]['config']
					max0 = 0.0
					min1 = 1.0
					if self.solve(q0,b):
						max0 = 1.0
					else:
						l,u = 0.0,1.0
						while u > l+resolution:
							m = 0.5*(l+u)
							x = workspace_interpolate(a,b,m)
							if self.solve(q0,x) is not None:
								l = m
							else:
								u = m
						max0 = l
					if self.solve(q1,a):
						min1 = 0.0
					else:
						l,u = 0.0,1.0
						while u > l+resolution:
							m = 0.5*(l+u)
							x = workspace_interpolate(a,b,m)
							if self.solve(q1,x) is not None:
								u = m
							else:
								l = m
						min1 = u
					if min1 < max0:
						eranges[vindex] = (min1,max0)
					else:
						#print "Split edge",es[vindex],"range",min1,max0
						eranges[vindex] = (0.5*(min1+max0),0.5*(min1+max0))
					#pts[-1] = vectorops.interpolate(getpoint(i),getpoint(j),(min1+max0)*0.5)
		if len(epts)==0: return [],[]

		#is the graph already simplicial?
		simplicial_already = True
		simplicial_already = (len(active)==2)
		if simplicial_already:
			if len(active)==2:
				simplices = triangles(G)
			else:
				simplices = tetrahedra(G)
			for s in simplices:
				if not any(v in vfailures for v in s):
					continue
				verts = [getpoint(v) for v in s]
				centroid = vectorops.div(vectorops.add(*verts),len(verts))
				if any(self.isResolvedNode(v) for v in s):
					others[len(pts)] = s
				pts.append(centroid)
		else:
			#compute edge-constrained simplicial decomposition of roadmap and add centroids
			#TODO: simplices returned by Delaunay are not edge constrained!
			hpts = [[getpoint(i)[a] for a in active] for i in G.nodes()]
			try:
				delaunay = scipy.spatial.Delaunay(hpts)
			except Exception as e:
				print "Error computing Delaunay triangulation"
				print e
				return None
			for simplex in delaunay.simplices:
				verts = [hpts[v] for v in simplex]
				centroid = vectorops.div(vectorops.add(*verts),len(verts))
				x = [0.0]*3
				for (a,z) in zip(active,centroid):
					x[a] = z
				others[len(pts)] = simplex
				pts.append(x)

		hpts = [[x[a] for a in active] for x in pts]

		try:
			delaunay = scipy.spatial.Delaunay(hpts)
		except Exception as e:
			print "Error computing Delaunay triangulation"
			print e
			return None
		faces = set()
		for simplex in delaunay.simplices:
			output = any(v in epts for v in simplex)
			if not output:
				continue

			if len(active)==3:
				for i in range(len(simplex)):
					tri = [v for j,v in enumerate(simplex) if i!=j]
					if not any([a in epts for a in tri]):
						continue
					"""
					#if it's a conflicting edge along the workspace graph, don't draw it
					conflictEdge = False
					vtri = [vs.get(v,-1) for v in tri]
					for a in tri:
						if a in epts and (a not in es or es[a][0] in vtri or es[a][1] in vtri ):
							conflictEdge = True
							break
					if conflictEdge:
						continue
					if sum(1 for v in vtri if v >= 0) >= 3:
						continue
					"""
					if sum([1 for a in tri if a in others]) != 2:
						continue
					faces.add(tuple(sorted(tri)))
			else:
				for i in range(len(simplex)):
					n = (i+1)%len(simplex)
					a,b = simplex[i],simplex[n]
					if not (a in epts or b in epts):
						continue
					#vseg = (vs.get(a,-1),vs.get(b,-1))
					if not (a in others or b in others):
						continue
					faces.add(tuple(sorted((a,b))))
		faces = list(faces)
		mesh = MeshTopology(hpts,faces)

		#optimize the coordinates of epts so that they 1) lie in their specified range and 2) minimize curvature
		#erase faces on centroids if they are bordered on all sides
		if len(active) == 2:
			for v in others:
				if v not in mesh.incidentFaces: continue
				if len(mesh.incidentFaces[v]) == 2:
					assert len(mesh.incidentVerts[v]) == 2
					#perform edge collapse: v->v1
					mesh.mergeVertex(v,mesh.incidentVerts[v][0])
				elif len(mesh.incidentFaces[v]) == 3:
					#completely surrounded
					pass
		elif len(active)==3:
			collapses = dict()
			for v in others:
				#continue
				if v not in mesh.incidentFaces: continue
				adjacentConflicts = [w for w in set(mesh.incidentVerts[v]) if w in epts]
				collapseto = None
				if len(adjacentConflicts) == 2:
					#collapse down along an edge to a vertex that shares 1 conflict but not the other
					thirdvs = set([w for w in mesh.incidentVerts[v] if w in others])
					for w in thirdvs:
						#share one conflict vertex plus a face of the simplex
						if sum(1 for u in adjacentConflicts if u in mesh.incidentVerts[w])==1 and len(set(others[v]) & set(others[w]))==3:
							if vectorops.distance(pts[v],pts[w]) < self.wRadius*0.5:
								collapseto = w
								break
				elif len(adjacentConflicts) >= 3:
					#collapse to any conflict node?
					#thirdvs = [w for w in incidentVerts[v] if w not in adjacentConflicts]
					#for w in thirdvs:
					#	if sum(1 for u in adjacentConflicts if u in incidentVerts[w])==1:
					#		collapseto = w
					#		break
					collapseto = None
					thirdvs = set([w for w in mesh.incidentVerts[v] if w in others])
					if len(thirdvs) == 3:
						mind = self.wRadius*0.5
						for w in adjacentConflicts:
							d = vectorops.distance(pts[v],pts[w])
							if d < mind:
								mind = d
								collapseto = w
				if collapseto is not None:
					collapses[v] = collapseto
			vertmap = dict()
			for v,w in collapses.iteritems():
				#print "Collapsing",v,"to",w
				numiters = 0
				if w in vertmap:
					continue
				#TEMP: test shifts
				#hpts[v] = vectorops.interpolate(hpts[v],hpts[w],0.5)
				#continue
				mesh.mergeVertex(v,w)
				vertmap[v] = w

		#optimize the boundary
		numPasses = 20
		attached = dict()
		for v in others:
			if v not in mesh.incidentVerts: continue
			#attach centroid points to single conflict points
			adjacentConflicts = [w for w in set(mesh.incidentVerts[v]) if w in epts]
			if len(adjacentConflicts) == 1:
				if adjacentConflicts[0] not in attached:
					attached[adjacentConflicts[0]] = []
				attached[adjacentConflicts[0]].append(v)

		for sweepCnt in xrange(numPasses):
			moveAmt = 0
			numMoved = 0
			numFixed = 0
			for v in es:
				if not v in mesh.incidentFaces: continue
				e = es[v]
				if e[0] < 0:
					numFixed += 1
					continue
				a,b = eranges[v]
				a = max(a,0.1)
				b = min(b,0.9)
				p = getpoint(e[0])
				q = getpoint(e[1])
				newpt = None
				if a==b:
					x = vectorops.interpolate(p,q,a)
					newpt = x
					numFixed+= 1
				else:
					numMoved += 1
					neighbors = [hpts[w] for w in mesh.incidentVerts[v] if w not in attached.get(v,[])]
					if len(neighbors) == 0:
						midpt2 = vectorops.interpolate(p,q,(a+b)*0.5)
					else:
						#optimize between p and q toward the average
						midpt = vectorops.div(vectorops.add(*neighbors),len(neighbors))
						midpt2 = p[:]
						for i,x in zip(active,midpt):
							midpt2[i] = x
					#find the point closest on the line pq to midpt
					d = vectorops.sub(q,p)
					try:
						u = vectorops.dot(d,vectorops.sub(midpt2,p))/vectorops.dot(d,d)
						u = min(max(u,a),b)
						pt = vectorops.interpolate(p,q,u)

						for i,x in zip(active,hpts[v]):
							midpt2[i] = x
						uorig = vectorops.dot(d,vectorops.sub(midpt2,p))/vectorops.dot(d,d)
						moveAmt += abs(uorig-u)*vectorops.norm(d)*(1+len(attached.get(v,[])))

						newpt = pt
					except ZeroDivisionError:
						print "Zero division error?",d
				if newpt is not None:
					oldpt = hpts[v]
					hpts[v] = [newpt[i] for i in active]
					for w in attached.get(v,[]):
						hpts[w] = vectorops.add(hpts[w],vectorops.sub(hpts[v],oldpt))
						numMoved += 1


			#print "Moved",numMoved,"points amount",moveAmt
			if moveAmt < 1e-3: break
		#now output only the hpts that are in the faces
		verts = []
		fmapped = []
		vmap = dict()
		for f in mesh.faces:
			if f == None: continue
			fmap = []
			for v in f:
				if v in vmap:
					fmap.append(vmap[v])
				else:
					vmap[v] = len(verts)
					fmap.append(len(verts))
					verts.append(hpts[v])
			if len(fmap) == 3:
				n = vectorops.cross(vectorops.sub(verts[fmap[1]],verts[fmap[0]]),vectorops.sub(verts[fmap[2]],verts[fmap[0]]))
				if n[2] < 0:
					fmap[1],fmap[2] = fmap[2],fmap[1]
			fmapped.append(fmap)
		return verts,fmapped

