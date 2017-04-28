from klampt import *
from klampt.math import vectorops,so3,se3
from utils.nearestneighbors import *
from utils.metric import *
from utils.disjointset import *
from collections import defaultdict
from ikdb.utils import mkdir_p
from utils.csp import *
from graph import RedundancyResolutionGraph
import networkx as nx
from utils.cgraph import CGraph
from utils.utils import *
import multiprocessing
import random
import math
import os
import copy
import time
import json




def nearest_point(G,nodes,weight=None,p=1):
	"""Finds the node closest to the given nodes on the graph, where the
	node-node distances are combined using an L_p metric.  If no node is
	reachable by all the nodes, this returns the one that is nearest to
	the maximum number of nodes.
	"""
	acc_lengths = defaultdict(float)
	acc_reachable = defaultdict(int)
	for n in nodes:
		L = nx.shortest_path_length(G,source=n,weight=weight)
		for m,v in L.iteritems():
			acc_reachable[m] += 1
			if p == 1:
				acc_lengths[m] += v
			elif p == float('inf'):
				acc_lengths[m] = max(acc_lengths[m],v)
			else:
				acc_lengths[m] += pow(v,p)
	allreachable = False
	numreachable = 0
	best = None
	for (n,r) in acc_reachable.iteritems():
		if r == len(nodes):
			numreachable = r
			allreachable = True
			break
		if r > numreachable:
			numreachable = r
			best = n
	if allreachable:
		#find the optimal node that is reachable by all specified nodes
		minlen = float('inf')
		best = None
		for (n,r) in acc_reachable.iteritems():
			if r != numreachable: 
				continue
			if acc_lengths[n] < minlen:
				minlen = acc_lengths[n]
				best = n
		assert best is not None
		return best
	else:
		raise NotImplementedError("TODO: nearest_point when graph is not connected")
		
class StatCollector:
	def __init__(self, problem = None):
		self.problem = problem
		self.Nw = 0
		self.method = None
		self.Nq = 0
		self.threads = 0
		self.cacheSize = 0
		self.visK = 0
		self.sumConfigs = 0      
		self.sumNumQccs = 0
		self.sumMaxSize = 0
		self.sumAvgSize = 0.0
		self.sumMinSize = 0
		self.configuredW = 0
		self.sampling_time = -1
		self.total_workspace_edges = 0
		self.nondegenerate_workspace_edges = 0
		self.edge_time = -1
		self.min_conflicts = -1
		self.csp_method = None
		self.csp_solve_time = -1
		self.collapse_time = -1
		
	def add(self, qccs):
		if len(qccs) > 0:
			self.configuredW += 1
			self.sumNumQccs += len(qccs)
			maxSize = 0
			sumSizes = 0
			minSize = float('inf')
			for rep in qccs:
				size = len(rep)
				if size > maxSize:
					maxSize = size
				sumSizes += size
				if size < minSize:
					minSize = size
			self.sumConfigs += sumSizes
			self.sumMaxSize += maxSize
			self.sumAvgSize += float(sumSizes)/float(len(qccs))
			self.sumMinSize += minSize
			
	def inc(self):
		self.nondegenerate_workspace_edges += 1
		
	def setProblem(self, problem):
		self.problem = problem
		
	def setNw(self, Nw):
		self.Nw = Nw
		
	def setMethod(self, method):
		self.method= method
		
	def setNq(self, Nq):
		self.Nq = Nq
		
	def setThreads(self, threads):
		self.threads = threads
	
	def setCacheSize(self, size):
		self.cacheSize = size
		
	def setVisK(self, visK):
		self.visK = visK
		
	def setSelfMotionCcsStats(self, stats):
		self.sumConfigs = stats['sumConfigs']
		self.sumNumQccs = stats['sumNumQccs']
		self.sumMaxSize = stats['sumMaxSize']
		self.sumAvgSize = stats['sumAvgSize']
		self.sumMinSize = stats['sumMinSize']
		self.configuredW = stats['configuredW']
		
	def setTotalWorkspaceEdges(self, edges):
		self.total_workspace_edges = edges
		
	def setSamplingTime(self, time):
		self.sampling_time = time
		
	def setEdgeTime(self, time):
		self.edge_time = time

	def setMinConflicts(self, min_conflicts):
		self.min_conflicts = min_conflicts

	def setCSPMethod(self, csp_method):
		self.csp_method = csp_method

	def setCSPSolveTime(self, csp_solve_time):
		self.csp_solve_time = csp_solve_time

	def setCollapseTime(self, collapse_time):
		self.collapse_time = collapse_time
			
	def getReport(self):
		report = {}
		report['Nw'] = self.Nw
		report['method'] = self.method
		report['Nq'] = self.Nq
		report['threads'] = self.threads
		report['cacheSize'] = self.cacheSize
		report['visK'] = self.visK
		report['configuredW'] = self.configuredW
		report['sampling_time'] = self.sampling_time
		report['total_workspace_edges'] = self.total_workspace_edges
		report['nondegenerate_workspace_edges'] = self.nondegenerate_workspace_edges
		report['edge_time'] = self.edge_time
		report['min_conflicts'] = self.min_conflicts
		report['csp_method'] = self.csp_method
		report['csp_solve_time'] = self.csp_solve_time
		report['collapse_time'] = self.collapse_time
		if self.configuredW == 0:
			report['avgNumConfigs'] = float('inf')
			report['avgNumQccs'] = float('inf')
			report['avgMaxSize'] = float('inf')
			report['avgAvgSize'] = float('inf')
			report['avgMinSize'] = float('inf')
		else:
			report['avgNumConfigs'] = float(self.sumConfigs)/float(self.configuredW)
			report['avgNumQccs'] = float(self.sumNumQccs)/float(self.configuredW)
			report['avgMaxSize'] = float(self.sumMaxSize)/float(self.configuredW)
			report['avgAvgSize'] = float(self.sumAvgSize)/float(self.configuredW)
			report['avgMinSize'] = float(self.sumMinSize)/float(self.configuredW)
		return report
		
	def writeReport(self):
		if self.problem == None:
			assert False, "stat collector problem must be initialized"
		out_file = "experiments/" + self.problem + "_k_experiments.csv"
		if not os.path.isdir("experiments"):
			os.makedirs("experiments")
		if not os.path.isfile(out_file):
			with open(out_file, 'a') as f:
				f.write('num_w,num_q_per_w,threads,cache_size,k,workspace_sampling_method,nondegenerate_workspace_nodes,nondegenerate_workspace_edges,total_workspace_edges,avg_num_configs_per_wnode,avg_num_self_motion_ccs,avg_max_size_self_motion_ccs,avg_avg_size_self_motion_ccs,avg_min_size_self_motion_ccs,q_sampling_m_construction_time,edge_tests_time,min_csp_conflicts,csp_method,csp_solve_time,collapse_time,visibility')
		
		report = self.getReport()
		with open(out_file, 'a') as f:
			f.write('\n'+str(report['Nw'])+","+str(report['Nq'])+","+str(report['threads'])+","+str(report['cacheSize'])+","+str(report['visK'])+","+str(report['method'])+","+str(report['configuredW'])+","+str(report['nondegenerate_workspace_edges'])+","+str(report['total_workspace_edges'])+","+str(report['avgNumConfigs'])+","+str(report['avgNumQccs'])+","+str(report['avgMaxSize'])+","+str(report['avgAvgSize'])+","+str(report['avgMinSize'])+","+str(report['sampling_time'])+","+str(report['edge_time'])+","+str(report['min_conflicts'])+","+str(report['csp_method'])+","+str(report['csp_solve_time'])+","+str(report['collapse_time']))


class RedundancyResolver:
	"""Performs the redundancy resolution procedure described in:

	K. Hauser, Continuous pseudoinversion of a multivariate function: application
	to global redundancy resolution, Workshop on the Algorithmic Foundations of Robotics, 2016.

	Attributes:
	- resolution: a RedundancyResolutionGraph structure defining the problem
	- robot: the robot to be resolved (must be same as resolution.robot)
	- Gw: a CGraph of the workspace graph, with an equivalent structure to the one in
	  resolution.Gw but with node attributes:
	  * 'params': the point in the workspace (same as in resolution.Gw)
	  * 'qlist': a list of configuration space node id's associated with this node
	  * 'qccs': (optional) a list-of-lists organizing the configurations in qlist, if self-motion
	     manifold construction is enabled.  Each entry is a list of indices into qlist that are
	     connected together with a self-motion manifold path.
	  * 'qcc_parents': (optional) a topological sort of the self-motion manifold connected
	    components.  A list of length len(qlist).  An entry of -1 indicates no parent.
	  and edge attributes:
	  * 'qlist': a list of feasible configuration space nodes associated with this edge.  These
	    should be interpreted as indexes accessing into the qlist of the incident nodes
	- Qsubgraph: (internal) the subgraph of Gq that should be exported to self.resolution.  This is a map 
	  windex => qindex where windex is the workspace index and qindex is the index into it's qlist.
	- QccSubgraph: (internal) the subgraph of the connected components that should be collapsed to
	  determine self.resolution.  Used if self-motion manifolds are being used.  This is a map
	  windex => qccindex where windex is the workspace index and qindex is the index into it's qcclist.

	Progress can be saved/loaded to a folder via save/load.
	"""
	def __init__(self,resolution,cache=20000,cacheloc=None):
		"""Arguments:
		- cache: the cache size.  Set cache = 0 to not use a cache.
		- cacheloc: if provided, the cache shelf location.  Useful if running multiple
		  programs simultaneously to avoid disk access conflicts.  Default is the
		  'shelf/' folder.
		"""
		self.resolution = resolution
		self.robot = resolution.robot
		self.cachesize = cache
		self.cacheloc = cacheloc
		self.Qsubgraph = dict()
		self.QccSubgraph = dict()
		self.initWorkspaceGraph()

	def initWorkspaceGraph(self, clear=False):
		"""Given a structure in resolution.Gw, creates the values of the workspace graph.  Gw and Gq are cleared.
		If clear=True, the resolution is cleared too.
		"""	
		if self.cachesize == 0:
			self.Gw = nx.Graph()
		elif self.cacheloc != None:
			self.Gw = CGraph(cachesize=self.cachesize,db_loc=self.cacheloc)
			self.Gw.db.printFrequency = 10
		else:
			self.Gw = CGraph(cachesize=self.cachesize)
			self.Gw.db.printFrequency = 10
		#copy Gw structure from resolution
		for i,d in self.resolution.Gw.nodes_iter(data=True):
			self.Gw.add_node(i,params=d['params'],qlist=[])
		for i,j,d in self.resolution.Gw.edges_iter(data=True):
			self.Gw.add_edge(i,j,qlist=[])
		if clear:
			self.resolution.clearResolution()

	def save(self,folder,close=False):
		mkdir_p(folder)
		self.Gw.save(os.path.join(folder,"Gw_cached.pickle"),close=close)
		self.resolution.save(folder)
		if len(self.Qsubgraph) > 0:
			f = open(os.path.join(folder,"Qsubgraph_cached.txt"),"w")
			for k,v in sorted(self.Qsubgraph.items()):
				f.write(str(k)+" "+str(v)+'\n')
			f.close()
		if len(self.QccSubgraph) > 0:
			f = open(os.path.join(folder,"QccSubgraph_cached.txt"),"w")
			for k,v in sorted(self.QccSubgraph.items()):
				f.write(str(k)+" "+str(v)+'\n')
			f.close()
		return

	def load(self,folder):
		try:
			self.Gw.load(os.path.join(folder,"Gw_cached.pickle"))
		except:
			self.resolution.load(folder)
			#TODO: decide whether we want to import the resolution into the graph
			self.fromResolutionToGraph(clear=True)
			for i,d in self.resolution.Gw.nodes_iter(data=True):
				d['qlist'] = []
			for i,j,d in self.resolution.Gw.edges_iter(data=True):
				d['qlist'] = []
			return
		self.resolution.load(folder)
		print "Loaded",self.numConfigurations(),"configs for",self.Gw.number_of_nodes(),'workspace nodes'
		try:
			f = open(os.path.join(folder,"Qsubgraph_cached.txt"),"r")
			self.Qsubgraph = dict()
			for line in f.readlines():
				key,value = line.split()
				key = int(key)
				value = int(value)
				assert self.Gw.has_node(key)
				assert value < len(self.Gw.node[key]['qlist']),"Invalid Qsubgraph entry %d => %d, # of configs at this node is %d"%(key,value,len(self.Gw.node[key]['qlist']))
				self.Qsubgraph[key] = value
			f.close()
		except IOError:
			import traceback
			traceback.print_exc()
			print "Could not load Qsubgraph_cached.txt"
			self.Qsubgraph = dict()
		try:
			f = open(os.path.join(folder,"QccSubgraph_cached.txt"),"r")
			self.QccSubgraph = dict()
			for line in f.readlines():
				key,value = line.split()
				key = int(key)
				value = int(value)
				assert self.Gw.has_node(key)
				assert value < len(self.Gw.node[key]['qlist'])
				self.QccSubgraph[key] = value
			f.close()
		except IOError:
			import traceback
			traceback.print_exc()
			print "Could not load QccSubgraph_cached.txt"
			self.QccSubgraph = dict()
		
		self.sanityCheck()
		"""
		#setup domain
		bmin = self.Gw.node[0]['params'][:]
		bmax = self.Gw.node[0]['params'][:]
		for i,d in self.Gw.nodes_iter(data=True):
			x = d['params']
			for j,v in enumerate(x):
				bmin[j] = min(bmin[j],v)
				bmax[j] = max(bmax[j],v)
		self.domain = (bmin,bmax)
		#setup nearest neighbors
		self.nnw = NearestNeighbors(L2Metric,'kdtree')
		for n,d in self.Gw.nodes_iter(data=True):
			self.nnw.add(d['params'],n)
		"""

	def sanityCheck(self):
		#sanity check
		for i,j,d in self.Gw.edges_iter(data=True):
			for iq,jq in d['qlist']:
				assert iq >= 0 and iq < len(self.Gw.node[i]['qlist'])
				assert jq >= 0 and jq < len(self.Gw.node[j]['qlist'])
		#check consistency with resolution
		for i,d in self.Gw.nodes_iter(data=True):
			assert self.resolution.Gw.has_node(i)
			assert self.resolution.Gw.node[i]['params'] == d['params']
			assert (i in self.Qsubgraph) == self.resolution.isResolvedNode(i),"Mismatch between resolution and Qsubgraph"
		for i,j in self.resolution.Gw.edges_iter():
			eqlist = self.Gw.edge[i][j]['qlist']
			insubgraph = self.edgeInResolution(i,j)
			assert self.resolution.isResolvedEdge(i,j) == insubgraph,"Mismatch between resolution edges and Qsubgraph edges: %d vs %d"%(int(self.resolution.isResolvedEdge(i,j)),int(insubgraph))
		for iw,iq in self.Qsubgraph.iteritems():
			for jw in self.Gw.neighbors(iw):
				assert self.edgeInResolution(iw,jw) == self.resolution.isResolvedEdge(iw,jw),"Mismatch between resolution edges and Qsubgraph edges"

	def edgeInResolution(self,iw,jw):
		if iw > jw:
			return self.edgeInResolution(jw,iw)
		return iw in self.Qsubgraph and jw in self.Qsubgraph and (self.Qsubgraph[iw],self.Qsubgraph[jw]) in self.Gw.edge[iw][jw]['qlist']

	def printStats(self):
		#cache friendlier version
		numNodes = 0
		numConfigs = 0
		numCSpaceEdges = 0
		numPotentialEdges = 0
		numEdgesWithConfigurations = 0
		numEdgesInResolution = 0
		sumdistances = 0
		sumwsdistances = 0
		print "**** REDUNDANCY RESOLUTION STATS *****"
		for i,d in self.Gw.nodes_iter(data=True):
			numConfigs += len(d['qlist'])
			if len(d['qlist']) > 0:
				numNodes += 1
				for j in self.Gw.neighbors(i):
					if i > j: continue
					if len(self.Gw.node[j]['qlist'])==0: continue
					ed = self.Gw.edge[i][j]
					numCSpaceEdges += len(ed['qlist'])
					numPotentialEdges += 1
					if len(ed['qlist'])>0:
						numEdgesWithConfigurations += 1
					if self.resolution.isResolvedEdge(i,j):
						numEdgesInResolution += 1
						sumdistances += self.robot.distance(self.resolution.Gw.node[i]['config'],self.resolution.Gw.node[j]['config'])
						sumwsdistances += workspace_distance(self.resolution.Gw.node[i]['params'],self.resolution.Gw.node[j]['params'])
		print "Roadmap has",self.Gw.number_of_nodes(),"workspace nodes and",self.Gw.number_of_edges(),"edges"
		print "  and",numConfigs,"configuration space nodes and",numCSpaceEdges,"edges"
		print "  of %d nodes, %d have no configuration, and %d / %d are not in redundancy resolution"%(self.Gw.number_of_nodes(),self.Gw.number_of_nodes()-numNodes,numNodes-self.resolution.resolvedCount,numNodes)
		print "  of %d edges, %d  are missing a configuration space edge, and %d / %d are not in redundancy resolution"%(numPotentialEdges,numPotentialEdges-numEdgesWithConfigurations,numPotentialEdges-numEdgesInResolution,numPotentialEdges)
		#if len(self.Qsubgraph) == 0:
		if not self.resolution.hasResolution():
			print "  Resolution is empty."
			return
		if numEdgesInResolution != 0:
			print "  Average configuration space distance:",sumdistances / numEdgesInResolution
			print "  Average configuration space distance / workspace distance:",sumdistances / sumwsdistances

	def numConfigurations(self):
		return sum(len(d['qlist']) for j,d in self.Gw.nodes_iter(data=True))

	def numCSpaceEdges(self):
		return sum(len(d['qlist']) for i,j,d in self.Gw.edges_iter(data=True))

	def addNode(self,iw,q):
		"""Adds a configuration node q corresponding to workspace node index iw.  Returns the index of the
		configuration space node"""
		assert iw < self.Gw.number_of_nodes()
		self.Gw.node[iw]['qlist'].append(q)
		if 'qccs' in self.Gw.node[iw]:
			#add an empty
			print "WARNING: addNode called to a workspace node where qccs attribute is already set up"
			self.Gw.node[iw]['qccs'].append([q])
		return len(self.Gw.node[iw]['qlist'])-1

	def addEdge(self,iw,iq,jw,jq):
		"""Adds a configuration space edge corresponding to workspace/configuration space node index (iw,iq) to
		index (jw,jq).
		"""
		if iw > jw:
			self.addEdge(jw,jq,iw,iq)
			return
		assert iw < self.Gw.number_of_nodes()
		assert jw < self.Gw.number_of_nodes()
		assert iq < len(self.Gw.node[iw]['qlist'])
		assert jq < len(self.Gw.node[jw]['qlist'])
		assert iw != jw,"Can't handle self-motion manifold edges in addEdge"
		self.Gw.edge[iw][jw]['qlist'].append((iq,jq))
	
	def hasEdge(self,iw,iq,jw,jq):
		if iw > jw:
			return self.hasEdge(jw,jq,iw,iq)
		return (iq,jq) in self.Gw.edge[iw][jw]['qlist']

	def removeEdge(self,iw,iq,jw,jq):
		assert iw != jw
		if iw > jw:
			self.removeEdge(jw,jq,iw,iq)
			return
		assert iw < self.Gw.number_of_nodes()
		assert jw < self.Gw.number_of_nodes()
		assert (iq,jq) in self.Gw.edge[iw][jw]['qlist']
		self.Gw.edge[iw][jw]['qlist'].remove((iq,jq))

	def testAndConnectAll(self,iw,jw):
		"""Tries to connect two workspace nodes in configuration space by testing all pairs of configurations."""
		res = False
		for na in xrange(len(self.Gw.node[iw]['qlist'])):
			for nb in xrange(len(self.Gw.node[jw]['qlist'])):
				if self.testAndConnect(iw,na,jw,nb):
					res = True
		return res

	def testAndConnect(self,iw,iq,jw,jq):
		"""Tries to connect two configurations in the configuration map"""
		if iw > jw:
			return self.testAndConnect(jw,jq,iw,iq)
		assert iq >= 0 and iq < len(self.Gw.node[iw]['qlist'])
		assert jq >= 0 and jq < len(self.Gw.node[jw]['qlist'])
		assert self.Gw.has_edge(iw,jw)
		d = self.Gw.edge[iw][jw]
		if (iq,jq) in d['qlist']:
			#already tested
			return True
		a = self.Gw.node[iw]['qlist'][iq]
		b = self.Gw.node[jw]['qlist'][jq]
		if self.resolution.validEdge(a,b,self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
			#add it
			self.addEdge(iw,iq,jw,jq)
			return True
		return False

	def testAndAddEdge(self,iw,iq,jw,q):
		"""Tries to connect an existing configuration to a new one.  If the connection
		is feasible, then the node is added and connected to iw,iq.  The index of
		the new node is returned.  If it's infeasible, None is returned."""
		assert iq < len(self.Gw.node[iw]['qlist'])
		assert jw < self.Gw.number_of_nodes()
		a = self.Gw.node[iw]['qlist'][iw]
		if self.resolution.validEdge(a,q,self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
			#add it
			jq = self.addNode(jw,q)
			self.addEdge(iw,iq,jw,jq)
			return jq
		return None

	def testAndConnectQccs(self,iw,jw,maxTestsPerCC = 10):
		"""Tries to link all configuration-space connected components of a workspace edge. 
		Returns number of connected CC's.

		Arguments:
		- (iw,jw): the workspace edge
		- maxTestsPerCC: the max number of tests run per CC run in increasing distance order.
		  An equivalent number is run with random configurations.
		"""
		assert iw != jw
		numConnected = 0
		ni = self.Gw.node[iw]
		nj = self.Gw.node[jw]
		for qicc in ni['qccs']:
			for qjcc in nj['qccs']:
				connected = False
				#sort by distance
				distances = [(self.robot.distance(ni['qlist'][a],nj['qlist'][b]),a,b) for a in qicc for b in qjcc]
				distances = sorted(distances)
				numtests = 0
				cnt = 0
				for dist,i,j in distances:
					cnt+=1
					if self.testAndConnect(iw,i,jw,j):
						connected=True
						break
					if cnt >= maxTestsPerCC:
						break
				if connected:
					numConnected += 1
					break
				if len(distances) > maxTestsPerCC:
					#random selection
					for (dist,i,j) in random.sample(distances[maxTestsPerCC:],min(maxTestsPerCC,len(distances)-maxTestsPerCC)):
						cnt+=1
						if self.testAndConnect(iw,i,jw,j):
							print "CC connection success with random connection"
							connected = True
							break
				if connected:
					numConnected += 1
				else:
					pass
					#print "CC connection failed after",cnt,"tests"
		return numConnected


	def getSelfMotionCcsStats(self):
		"""Returns statistics describing the connected components in the self-motion graph
		of each workspace node."""
		stats = {}
		sumConfigs = 0        
		sumNumQccs = 0
		sumMaxSize = 0
		sumAvgSize = 0.0
		sumMinSize = 0
		configuredW = 0
		
		for w,d in self.Gw.nodes_iter(data=True):
			if len(d['qlist']) > 0:
				configuredW += 1
			try:
				if len(d['qccs']) > 0:
					sumNumQccs += len(d['qccs'])
					maxSize = 0
					sumSizes = 0
					minSize = float('inf')
					for cc in d['qccs']:
						size = len(cc)
						if size > maxSize:
							maxSize = size
						sumSizes += size
						if size < minSize:
							minSize = size
					sumConfigs += sumSizes
					sumMaxSize += maxSize
					sumAvgSize += float(sumSizes)/float(len(d['qccs']))
					sumMinSize += minSize
			except KeyError:
				if len(d['qlist']) > 0:
					sumNumQccs += len(d['qlist'])
					sumConfigs += len(d['qlist'])
					sumMaxSize += 1
					sumAvgSize += 1
					sumMinSize += 1
		
		stats['sumConfigs'] = sumConfigs    
		stats['sumNumQccs'] = sumNumQccs
		stats['sumMaxSize'] = sumMaxSize
		stats['sumAvgSize'] = sumAvgSize
		stats['sumMinSize'] = sumMinSize
		stats['configuredW'] = configuredW
		if configuredW == 0:
			stats['avgConfigs'] = float('inf')
			stats['avgNumQccs'] = float('inf')
			stats['avgMaxSize'] = float('inf')
			stats['avgAvgSize'] = float('inf')
			stats['avgMinSize'] = float('inf')
		else:
			stats['avgConfigs'] = float(sumConfigs)/float(configuredW)
			stats['avgNumQccs'] = float(sumNumQccs)/float(configuredW)
			stats['avgMaxSize'] = float(sumMaxSize)/float(configuredW)
			stats['avgAvgSize'] = float(sumAvgSize)/float(configuredW)
			stats['avgMinSize'] = float(sumMinSize)/float(configuredW)
		return stats

	def constructSelfMotionManifold(self,iw,k=None):
		"""Constructs self motion manifold for a single workspace node"""
		d = self.Gw.node[iw]
		#construct visibility prm
		wqnodes = len(d['qlist'])
		ccs = DisjointSet()
		Gccs = nx.Graph()
		for i in range(wqnodes):
			ccs.add(i)
			Gccs.add_node(i)
		qis = d['qlist']
		xi = d['params']
		if k is not None:
			nnq = NearestNeighbors(L2Metric,'kdtree')
			for i,qi in enumerate(qis):
				nnq.add(qi,i)
			for i,qi in enumerate(qis):
				knn = nnq.knearest(qi,k)
				for qj,j in knn:
					if not ccs.same(i,j):
						if self.resolution.validEdge(qi,qj,xi,xi):
							ccs.merge(i,j)
							Gccs.add_edge(i,j)
		else:
			#visibility prm with all-pairs connections 
			for i,qi in enumerate(qis):
				distances = []
				for j in xrange(i):
					if not ccs.same(i,j):
						distances.append((self.robot.distance(qi,qis[j]),j))
				for dists,j in sorted(distances):
					if not ccs.same(i,j):
						if self.resolution.validEdge(qi,qis[j],xi,xi):
							ccs.merge(i,j)
							Gccs.add_edge(i,j)

		#now create 'qccs' and 'qcc_parents'
		components = dict()
		for rep in ccs.getReps():
			components[rep] = [v for v in ccs.getSetRef(rep)]
		d['qccs'] = components.values()

		parentList = [-1]*wqnodes
		for e in nx.dfs_edges(Gccs):
			parentList[e[1]]=e[0]
		d['qcc_parents'] = parentList
		return len(d['qccs'])
		
	def constructSelfMotionManifolds(self,k=None):
		"""Create self-motion C-space edges to minimize number of C-space connected
		components for all workspace points.

		If k is provided, does k-nearest neighbor connections.  Otherwise, tries all connections
		"""
		print "Constructing self-motion connected components"
		c = ProgressUpdater(self.Gw.number_of_nodes(),5)
		start = time.time()
		for w in self.Gw.nodes_iter():
			c.update()
			self.constructSelfMotionManifold(w,k)
		end = time.time()
		stats = self.getSelfMotionCcsStats()
		print "CC statistics: %f configs, %f components, size max %f avg %f min %f"%(stats['avgNumConfigs'], stats['avgNumQccs'], stats['avgMaxSize'], stats['avgAvgSize'], stats['avgMinSize'])
		return end - start

	def connectConfigurationSpaceOnEdge(self,iw,jw,istart=0,jstart=0):
		#standard method
		numConnected = 0
		for na,qa in enumerate(self.Gw.node[iw]['qlist']):
			for nb,qb in enumerate(self.Gw.node[jw]['qlist']):
				if na < istart and nb < jstart: continue
				if self.testAndConnect(iw,na,jw,nb):
					numConnected += 1
		return numConnected

		
	def connectConfigurationSpace(self,visibility=False,k=None,counts=None):
		self_motion_time = self.constructSelfMotionManifolds(k=k) if visibility else 0
				
		#create C-space edges for each workspace edge as long as the edge is feasible and
		#the movement monotonically decreases the distance
		print "Testing configuration space edges"
		numNondegenerateEdges = 0
		numPotentialEdges = 0
		minVisibilityEdgesChecked = 0
		numActualizedEdges = 0
		if counts is None:
			counts = [0]*self.Gw.number_of_nodes()
		else:
			assert visibility == False,"Can't have counts (incremental construction) plus visibility graph construction"
		start = time.time()
		for (iw,jw) in self.Gw.edges_iter():
			numPotentialEdges += (len(self.Gw.node[iw]['qlist'])*len(self.Gw.node[jw]['qlist']))-(counts[iw]*counts[jw])
			if visibility:
				minVisibilityEdgesChecked += len(self.Gw.node[iw]['qccs'])*len(self.Gw.node[jw]['qccs'])
		print "Potentially " + str(numPotentialEdges) + " edges to check"
		if visibility:
			print "Visibility method will check at least" + u"\u2248" + str(minVisibilityEdgesChecked) + " edges"
		c = ProgressUpdater(self.Gw.number_of_edges(),1)
		for (iw,jw) in self.Gw.edges_iter():
			c.update()
			if visibility:
				numConnections = self.testAndConnectQccs(iw,jw)
			else:
				numConnections = self.connectConfigurationSpaceOnEdge(iw,jw,counts[iw],counts[jw])
			if numConnections > 0:
				numNondegenerateEdges += 1
			numActualizedEdges += numConnections
		end = time.time()
		inter_wnode_time = end - start
		c.done()
		#check for isolated nodes
		print numNondegenerateEdges,"of",self.Gw.number_of_edges(),"workspace edges are nondegenerate"
		if visibility:
			print "The self-motion graph of each workspace node has an average number of",round(self.getSelfMotionCcsStats()['avgNumQccs'],2),"connected components"
			print numActualizedEdges,"edges out of possible",minVisibilityEdgesChecked
		else:
			"""
			#TODO: print stats with new cached Gw structure
			ccs = list(nx.connected_components(self.Gq))
			print "Configuration space graph has",self.Gq.number_of_nodes(),"nodes and",len(ccs),"connected components"
			print "CC sizes:"
			for i in range(len(ccs)):
				if len(ccs[i]) == 1:
					print 
					print len(ccs)-i,"isolated nodes"
					break
				print "  ",len(ccs[i]),
			print
			"""
			print numActualizedEdges,"edges out of",numPotentialEdges,"checked are feasible"
		return numNondegenerateEdges, self.Gw.number_of_edges(), self_motion_time, inter_wnode_time

	def sampleConfigurationSpace(self,NConfigsPerPoint=100,connect=True,biased=False,seedFromAdjacent=True,visibility=False,k=None):
		print "Sampling",NConfigsPerPoint,"configurations at each workspace point"
		start_node_counts = [len(d['qlist']) for j,d, in self.Gw.nodes_iter(data=True)]
		start_node_count = sum(start_node_counts)
		if start_node_count > 0:
			print "   Already have",start_node_count,"configurations in graph"
		configuration_count = start_node_count
		number_of_nodes_with_configurations = sum(1 for c in start_node_counts if c > 0)
		numConnectedEdges = 0
		totalAdjacentSamples = 0
		totalRandomSamples = 0
		successfulAdjacentSamples = 0
		successfulRandomSamples = 0
		c = ProgressUpdater(self.Gw.number_of_nodes(),5)
		t_connect = 0
		start = time.time()
		for i,d in self.Gw.nodes_iter(data=True):
			c.update()
			#adjacentWithConfigs = [j for j in self.Gw.neighbors(i) if j < i and start_node_counts[j] > 0]
			adjacentWithConfigs = [j for j in self.Gw.neighbors(i) if j < i and len(self.Gw.node[j]['qlist']) > 0]
			fracAdjacent = float(len(adjacentWithConfigs)) / float(len(self.Gw.neighbors(i)))
			numAdjacentSamples = int(fracAdjacent*NConfigsPerPoint)
			numRandomSamples = NConfigsPerPoint - numAdjacentSamples 
			if biased and start_node_count != 0:
				#numconfigs * average # of configs per workspace point / # of configs for this point
				if len(d['qlist']) == 0:
					numRandomSamples = 0
				else:
					numRandomSamples = NConfigsPerPoint * configuration_count / (number_of_nodes_with_configurations * len(d['qlist']))
			#print fracAdjacent,"adjacent"

			self.resolution.setIKProblem(d['params'])

			#try starting from existing neighboring configurations
			self.resolution.ikSolverParams.startRandom = False
			numFailed = 0
			if seedFromAdjacent:
				for it in xrange(numAdjacentSamples):
					totalAdjacentSamples += 1
					j = random.choice(adjacentWithConfigs)
					qstart = random.choice(self.Gw.node[j]['qlist'])
					self.robot.setConfig(qstart)
					q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
					if q is not None:
						successfulAdjacentSamples += 1
						#this checking is useful for reducing # of edge checks esp for nonredundant robots
						same = False
						for iq in d['qlist']:
							if self.robot.distance(q,iq) < 1e-2:
								same = True
								break
						if not same:
							self.addNode(i,q)
							configuration_count += 1
						else:
							#print "Same config"
							pass
					else:
						numFailed += 1

			self.resolution.ikSolverParams.startRandom = True
			for it in xrange(numRandomSamples+numFailed):
				#solver will sample randomly
				totalRandomSamples += 1
				q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
				if q is not None:
					successfulRandomSamples += 1
					#this checking is useful for reducing # of edge checks esp for nonredundant robots
					same = False
					for iq in d['qlist']:
						if self.robot.distance(q,iq) < 1e-2:
							same = True
							break
					if not same:
						self.addNode(i,q)
						configuration_count += 1
					else:
						#print "Same config"
						pass

			if visibility:
				self.constructSelfMotionManifold(i,k)
			
			if connect:
				cstart = time.time()
				if len(d['qlist']) == 0:
					continue
				for j in self.Gw.neighbors(i):
					if j > i or len(self.Gw.node[j]['qlist']) == 0:
						continue
					if visibility:
						numConnections = self.testAndConnectQccs(i,j)
					else:
						numConnections = self.connectConfigurationSpaceOnEdge(i,j,start_node_counts[i],start_node_counts[j])
					numConnectedEdges += numConnections
				cend = time.time()
				t_connect += cend - cstart

		end = time.time()
		t_sample = end - start - t_connect
		c.done()
		print "Sampled",configuration_count-start_node_count,"of possible",self.Gw.number_of_nodes()*NConfigsPerPoint
		print successfulAdjacentSamples,"/",totalAdjacentSamples,"succeeded from adjacent configs"
		print successfulRandomSamples,"/",totalRandomSamples,"succeeded from random configs"
		if visibility:
			print "The self-motion graph of each workspace node has an average number of",round(self.getSelfMotionCcsStats()['avgNumQccs'],2),"connected components"

		#Batch method: less cache efficient than checking in the above loop
		#if connect:
		#	self.connectConfigurationSpace(counts=start_node_counts)
		if connect:
			numPotentialEdges = 0
			for (iw,jw) in self.Gw.edges_iter():
				numPotentialEdges += (len(self.Gw.node[iw]['qlist'])*len(self.Gw.node[jw]['qlist']))-(start_node_counts[iw]*start_node_counts[jw])
			print "Connected",numConnectedEdges,"out of possible",numPotentialEdges
		
		return t_sample, t_connect
		
	def parallelQSamplingAtP(self, params, NConfigs=100, selfmotion=True, k=None):
		configs = []
		ccs = DisjointSet()
		
		#sample configurations
		self.resolution.setIKProblem(params)
		self.resolution.ikSolverParams.startRandom = True
		for it in xrange(NConfigs):
			#solver will sample randomly
			q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
			if q is not None:
				#this checking is useful for reducing # of edge checks esp for nonredundant robots
				same = False
				for qj in configs:
					if self.robot.distance(q,qj) < 1e-2:
						same = True
						break
				if not same:
					ccs.add(len(configs))
					configs.append(q)
				else:
					#print "Same config"
					pass
				
		#construct self motion manifolds
		wqnodes = len(configs)
		Gccs = nx.Graph()
		for i in range(wqnodes):
			Gccs.add_node(i)
		if selfmotion and k is not None:
			nnq = NearestNeighbors(L2Metric,'kdtree')
			for iq in range(wqnodes):
				nnq.add(configs[iq],iq)
			for iq in range(wqnodes):
				knn = nnq.knearest(configs[iq],k)
				for q,jq in knn:
					if not ccs.same(iq,jq):
						if self.resolution.validEdge(configs[iq],configs[jq],params,params):
							#keep track of qedges?
							ccs.merge(iq,jq)
							Gccs.add_edge(iq,jq)
		elif selfmotion:
			#visibility prm with all-pairs connections 
			for i,qi in enumerate(configs):
				distances = []
				for j in xrange(i):
					if not ccs.same(i,j):
						distances.append((self.robot.distance(qi,configs[j]),j))
				for dist,j in sorted(distances):
					if not ccs.same(i,j):
						if self.resolution.validEdge(qi,configs[j],params,params):
							ccs.merge(i,j)
							Gccs.add_edge(i,j)
			"""
			for iq in range(wqnodes):
				for jq in xrange(iq):
					if not ccs.same(iq,jq):
						if self.resolution.validEdge(configs[iq],configs[jq],params,params):
							#keep track of qedges?
							ccs.merge(iq,jq)
							Gccs.add_edge(iq,jq)
			"""			
		components = dict()
		for rep in ccs.getReps():
			components[rep] = [v for v in ccs.getSetRef(rep)]
		qccs = components.values()

		parentList = [-1]*wqnodes
		for e in nx.dfs_edges(Gccs):
			parentList[e[1]]=e[0]
		qcc_parents = parentList
		
		return configs, qccs, qcc_parents
		
	def parallelEdgeConnection(self, iconfigs, iqccs, iparams, jconfigs, jqccs, jparams, maxTestsPerCC=20):
		"""This method assumes that the workspace index of node i is less than that of node j."""
		qlist = []
		
		for qicc in iqccs:
			for qjcc in jqccs:
				connected = False
				#sort by distance
				distances = [(self.robot.distance(iconfigs[a],jconfigs[b]),a,b) for a in qicc for b in qjcc]
				distances = sorted(distances)
				numtests = 0
				cnt = 0
				for dist,i,j in distances:
					cnt+=1
					if self.resolution.validEdge(iconfigs[i],jconfigs[j],iparams,jparams):
						#add it
						qlist.append((i,j))
						connected = True
						break
					if cnt >= maxTestsPerCC:
						break
				if connected:
					continue
				if len(distances) > maxTestsPerCC:
					#random selection
					for (dist,i,j) in random.sample(distances[maxTestsPerCC:],min(maxTestsPerCC,len(distances)-maxTestsPerCC)):
						if self.resolution.validEdge(iconfigs[i],jconfigs[j],iparams,jparams):
							#add it
							qlist.append((i,j))
							print "CC connection success with random connection"
							connected = True
							break
						
		return qlist

	def makeCSP(self,visibility=False):
		"""Creates a CSP whereby each workspace node is a variable and configs are values.
		Each neighboring node must satisfy the constraint that its workspace edge corresponds to
		a feasible configuration space edge.
		"""
		csp = CSP()
		if visibility:
			for n,d in self.Gw.nodes_iter(data=True):
				if len(d['qlist']) > 0:
					csp.addVariable(range(len(d['qccs'])),name="qcc_"+str(n))
			for a,b,d in self.Gw.edges_iter(data=True):
				assert a < b,"Edges aren't ordered in networkx data structure?"
				if len(d['qlist']) == 0: continue
				validValues = set()
				accmap = dict()
				bccmap = dict()
				for i,cc in enumerate(self.Gw.node[a]['qccs']):
					for j in cc:
						accmap[j] = i
				for i,cc in enumerate(self.Gw.node[b]['qccs']):
					for j in cc:
						bccmap[j] = i
				for (iq,jq) in d['qlist']:
					assert iq >= 0 and iq < len(self.Gw.node[a]['qlist'])
					assert jq >= 0 and jq < len(self.Gw.node[b]['qlist'])
					assert iq in accmap
					assert jq in bccmap
					validValues.add((accmap[iq],bccmap[jq]))
				csp.addConstraint(validValues,("qcc_"+str(a),"qcc_"+str(b)))
		else:
			for n,d in self.Gw.nodes_iter(data=True):
				if len(d['qlist']) > 0:
					csp.addVariable(range(len(d['qlist'])),name="q_"+str(n))
			for a,b,d in self.Gw.edges_iter(data=True):
				assert a < b,"Edges aren't ordered in networkx data structure?"
				if len(d['qlist']) == 0: continue
				validValues = set(d['qlist'])
				csp.addConstraint(validValues,("q_"+str(a),"q_"+str(b)))
		return csp

	def collapseQccSubgraphNode(self,iw,fixed=None):
		"""Tries a simple node-centric collapse so that the representative configuration in the associated
		connected component in QccSubgraph is connected to all anchors in neighboring entries."""
		assert iw in self.QccSubgraph
		if fixed is None:
			fixed = dict()
		ccindex = self.QccSubgraph[iw]
		#the set of configs in same self-motion manifold (SMM)
		manifold = self.Gw.node[iw]['qccs'][ccindex]
		#the dict mapping configs in SMM to lists of their connections at neighboring workspace points
		anchors = dict((iq,[]) for iq in manifold)
		#the set of all attached configs at neighboring workspace points
		neighbors = set()
		#the dictionary mapping workspace neighbors to the configuration space edge connecting them
		wneighbors = dict()
		fixedNeighbors = dict()
		#print 
		#print "**** Beginning collapse of node",iw,"****"
		smanifold = set(manifold)
		for jw in self.Gw.neighbors(iw):
			if jw not in self.QccSubgraph: continue
			jmanifold = self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]]
			jsmanifold = set(jmanifold)
			added = False
			if iw < jw:
				for (iq,jq) in self.Gw.edge[iw][jw]['qlist']:
					if iq in smanifold and jq in jsmanifold:
						if jw in fixed and fixed[jw] != jq:
							continue
						if jw in wneighbors:
							#print "Workspace edge (%d,%d) duplicated, C-space edges (%d,%d) and (%d,%d) exist"%(iw,jw,iq,jq,wneighbors[jw][0],wneighbors[jw][1])
							continue
						wneighbors[jw] = (iq,jq)
						anchors[iq].append((jw,jq))
						neighbors.add((jw,jq))
			else:
				for (jq,iq) in self.Gw.edge[jw][iw]['qlist']:
					if iq in smanifold and jq in jsmanifold:
						if jw in fixed and fixed[jw] != jq:
							continue
						if jw in wneighbors:
							#print "Workspace edge (%d,%d) duplicated, C-space edges (%d,%d) and (%d,%d) exist"%(iw,jw,iq,jq,wneighbors[jw][0],wneighbors[jw][1])
							continue
						wneighbors[jw] = (iq,jq)
						anchors[iq].append((jw,jq))
						neighbors.add((jw,jq))
			if jw in fixed:
				#print "Neighbor",jw,"fixed at config",fixed[jw]
				#print "  My neighbors",neighbors
				if (jw,fixed[jw]) not in neighbors:
					#print "  There doesn't seem to be a C-space edge to it?"
					#print "  ",self.Gw.edge[iw][jw]['qlist']
					#if jw in wneighbors:
					#	print "But there is an edge in the C-space list!",wneighbors[jw]
					#	raw_input()
					pass
				else:
					fixedNeighbors[jw] = fixed[jw]
		
		#if we add new edges, we'll take them out before this is complete so that only edges between anchors are kept.
		newEdges = set()

		if len(anchors) == 0: #TODO should we ever reach this case?
			return manifold[0],"No anchors found"
		elif len(anchors)  == 1:
			return anchors.keys()[0],"Single anchor found to "+' '.join(str(x) for x in wneighbors.keys())
		#try to find a hub connecting edges to all adjacent workspace nodes, if possible
		#print "Anchors for",iw,":",anchors.keys()
		for iq in anchors:
			isHub = True
			for jw,jq in neighbors:
				if not self.hasEdge(iw,iq,jw,jq):
					if self.testAndConnect(iw, iq, jw, jq):
						newEdges.add((iw,iq,jw,jq))
						anchors[iq].append((jw,jq))
					else:
						isHub = False							
						break
			if isHub:
				#quick exit: there exists a hub, just use that
				#print "Found a hub for workspace node",iw,":",iq
				#remove all non-hub edges
				"""
				toPrune = [iq2 for iq2 in anchors if iq2 != iq]
				for iq2 in toPrune:
					#remove jq's inter-wnode edges
					cut = []
					for (jw,jq) in neighbors:
						if self.hasEdge(iw,iq2,jw,jq):
							cut.append((jw,jq))
					for kw,kq in cut:
						self.removeEdge(iw, iq2, kw, kq)
				"""
				return iq,"Found a single hub configuration connecting to "+" ".join(str(x) for x in neighbors)
	
		print "Failed to collapse self-motion manifold of node",iw,":",self.Gw.node[iw]['params'],"using hub heuristic"
		#collapse using a set-cover algorithm
		result = set_cover(anchors,neighbors)
		print "  Set cover reduced from",len(anchors),"anchors to",len(result)
		anchors = result

		"""
		#test centroid connections?
		midpt = robot_average(self.robot,[self.Gw.node[iw]['qlist'][iq] for iq in anchors])
		ix = self.Gw.node[iw]['params']
		midpt = self.resolution.solve(midpt,ix)
		if midpt != None:
			print "Trying midpoint connections..."
			numConnections = 0
			for (jw,jq) in neighbors:
				if self.resolution.validEdge(midpt,self.Gw.node[jw]['qlist'][jq],ix,self.Gw.node[jw]['params']):
					numConnections += 1
			print numConnections,"of",len(neighbors),"connections succeeded"
			if numConnections == len(neighbors):
				print "TODO: connect midpoint to neighbors"
		else:
			print "Midpoint configuration can't be solved"
		"""

		why = "Unknown reason"
		#now try to shrink down the distance between anchors via a spanning tree algorithm
		try:
			Gcc = nx.Graph()
			for n in manifold:
				Gcc.add_node(n)
			toposort = self.Gw.node[iw]['qcc_parents']
			for n in manifold:
				if toposort[n] >= 0:
					assert Gcc.has_node(toposort[n])
					Gcc.add_edge(n,toposort[n])
			components = 0
			for component in nx.connected_components(Gcc):
				components += 1
			assert components == 1,"Somehow the connected components computed by DisjointSet don't give a connected graph?"
			nroot = nearest_point(Gcc,anchors.keys())
			P = nx.predecessor(Gcc,nroot)
			#make the anchors graph
			Ganchors = nx.DiGraph()
			for n in anchors:
				Ganchors.add_node(n)
				p = P[n]
				while len(p) != 0 and p[0] not in anchors:
					if Ganchors.has_node(p[0]):
						break
					Ganchors.add_node(p[0])
					Ganchors.add_edge(p[0],n)
					p = P[p[0]]
			print Ganchors.number_of_nodes(),"nodes intervening between anchors"
			#collapse along edges in Ganchors
			for iq in reversed(nx.topological_sort(Ganchors)):
				#try attaching all anchors for iq to P[iq]
				if iq not in anchors:
					continue
				if len(P[iq]) == 0:
					#hit the root
					continue
				pq = P[iq][0]
				connected = True
				for jw,jq in anchors[iq]:
					if not self.hasEdge(iw,pq,jw,jq):
						if self.testAndConnect(iw, pq, jw, jq):
							newEdges.add((iw,pq,jw,jq))
						else:
							connected = False
				if connected:
					#print "Moving anchors up from",iq,"to",pq
					if pq not in anchors:
						anchors[pq] = set()
					anchors[pq] |= anchors[iq]
					del anchors[iq] 
			print "Collapsed further to",len(anchors),"anchors using spanning tree"
			if len(anchors) == 1:
				why = "Spanning tree algorithm collapsed successfully"
		except Exception as e:
			import traceback
			traceback.print_exc()
			print "Manifold:",manifold
			toposort = self.Gw.node[iw]['qcc_parents']
			print "Topological sort:",[toposort[v] for v in manifold]
			print "Exception occurred during spanning tree algorithm, continuing"

		bestq,bestanchored = sorted([(len(anchored),(iq,anchored)) for iq,anchored in anchors.iteritems()])[-1][1]
		#try neighbors that are not anchored
		disconnectedWNeighbors = set()
		for iq,anchored in anchors.iteritems():
			if iq != bestq:
				for (jw,jq) in anchored:
					disconnectedWNeighbors.add(jw)	
		if len(disconnectedWNeighbors) > 0:
			print "Node",iw,"best anchor",bestq,"has",len(bestanchored),"anchors"
			print "  and",len(disconnectedWNeighbors),"neighbors are disconnected from best"
			numConnections = 0
			for jw in disconnectedWNeighbors:
				jmanifold = self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]]
				jsmanifold = set(jmanifold)
				connected = False
				for jq in jmanifold:
					if jw in fixed and fixed[jw] != jq:
						continue
					if self.hasEdge(iw,bestq,jw,jq):
						#print "Hmmm. already have an edge from best config to neighbor",jw,"config",jq
						bestanchored.add((jw,jq))
						neighbors.add((jw,jq))
						numConnections += 1
						connected = True
						break
					if self.testAndConnect(iw, bestq, jw, jq):
						#print "Successfully connected from",iw,"anchor to",jw,"config",jq
						bestanchored.add((jw,jq))
						neighbors.add((jw,jq))
						numConnections += 1
						connected = True
						break
				if not connected:
					#print "Could not connect to neighbor",jw,"from config",bestq
					pass
			print "Connected",numConnections,"of",len(disconnectedWNeighbors),"greedily"
			if numConnections == len(disconnectedWNeighbors):
				anchors = {bestq:bestanchored}
				why = "Greedy hub identification algorithm collapsed successfully"

		#do we keep around the edges?
		"""
		cspaceEdgeList = dict()
		for jw in self.Gw.neighbors(iw):
			cspaceEdgeList[jw] = []
			if iw < jw:
				for (iq,jq) in self.Gw.edge[iw][jw]['qlist']:
					if (jw,jq) in neighbors:
						cspaceEdgeList[jw].append((iq,jq))
			else:
				for (jq,iq) in self.Gw.edge[jw][iw]['qlist']:
					if (jw,jq) in neighbors:
						cspaceEdgeList[jw].append((iq,jq))
		for jw in cspaceEdgeList:
			if jw in fixed:	
				print "Edges to fixed neighbor",jw,":",cspaceEdgeList[jw]
		print "*** My anchors ****"
		for jw in anchors:
			print "  ",jw,":",sorted(list(anchors[jw]))
		for jw,edges in cspaceEdgeList.iteritems():
			if len(edges) > 1:
				anchored = []
				nonanchored = []
				for (iq,jq) in edges:
					if iq in anchors and (jw,jq) in anchors[iq]:
						anchored.append((iq,jq))
					else:
						nonanchored.append((iq,jq))
				edges = anchored + nonanchored
				if jw in fixedNeighbors:
					print "Have fixed edge",iw,jw," to anchor",fixedNeighbors[jw]
					print "keeping",edges[0]
					print "Removing",edges[1:]
				for (iq,jq) in edges[1:]:
					self.removeEdge(iw, iq, jw, jq)
		"""
		if len(anchors) > 1:
			#can't resolve this CC, fix it using full resolution afterward
			return None,"Cannot be resolved, "+str(len(anchors))+" anchors left"
		iq = anchors.keys()[0] 
		for jw,jq in fixedNeighbors.iteritems():
			if iw < jw:
				if (iq,jq) not in self.Gw.edge[iw][jw]['qlist']:
					print "Edge deleted from anchor to fixed neighbor anchor"
					print "Workspace edge",iw,jw
					print "Picked C-space edge",(iq,jq)
					print "Original edge",wneighbors[jw]
					print "Qlist",self.Gw.edge[iw][jw]['qlist']
					raw_input()
			else:
				if (jq,anchors.keys()[0]) not in self.Gw.edge[jw][iw]['qlist']:
					print "Edge deleted from anchor to fixed neighbor anchor"
					print "Workspace edge",iw,jw
					print "Picked C-space edge",(iq,jq)
					print "Original edge",wneighbors[jw]
					print "Qlist",self.Gw.edge[iw][jw]['qlist']
					raw_input()

		for jw in wneighbors:
			jmanifold = self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]]
			jsmanifold = set(jmanifold)
			connected = False
			if iw < jw:
				for (iq2,jq) in self.Gw.edge[iw][jw]['qlist']:
					if iq == iq2 and jq in jsmanifold:
						connected=True
			else:
				for (jq,iq2) in self.Gw.edge[jw][iw]['qlist']:
					if iq == iq2 and jq in jsmanifold:
						connected=True
			if not connected:
				print "Not connected to neighbor",jw,"?"
				print "Chosen anchor",iq
				print "Original edge",wneighbors[jw]
				print self.Gw.edge[iw][jw]['qlist']
				raw_input()

		return anchors.keys()[0],why


	def collapseQccSubgraph(self, parallel=False, threads=None, tasks=None, results=None):
		start = time.time()
		self.Qsubgraph = dict()
		totalCollapsed = 0
		triedToCollapse = 0
		for iw,d in self.Gw.nodes_iter(data=True):
			if len(d['qlist']) > 0:
				assert iw in self.QccSubgraph,"Workspace node %d not covered by QccSubgraph?"%(iw,)

		assignedWEdges = set()
		resolvedAnchors = dict()
		resolutionFailures = set()
		order = sorted(self.QccSubgraph.keys())
		reasons = dict()
		#build topology of resolved connected components
		QccEdges = set()
		for iw in order:
			manifold = self.Gw.node[iw]['qccs'][self.QccSubgraph[iw]]
			smanifold = set(manifold)
			for jw in self.Gw.neighbors(iw):
				if jw not in self.QccSubgraph: continue
				jmanifold = self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]]
				jsmanifold = set(jmanifold)
				added = False
				if iw < jw:
					for (iq,jq) in self.Gw.edge[iw][jw]['qlist']:
						if iq in smanifold and jq in jsmanifold:
							QccEdges.add((iw,jw))
							break
				else:
					for (jq,iq) in self.Gw.edge[jw][iw]['qlist']:
						if iq in smanifold and jq in jsmanifold:
							QccEdges.add((jw,iw))
							break
		#try collapsing individual nodes
		for iw in order:
			if iw not in self.QccSubgraph: continue
			collapsed,why = self.collapseQccSubgraphNode(iw,resolvedAnchors)
			reasons[iw] = why
			triedToCollapse += 1
			if collapsed is not None:
				resolvedAnchors[iw] = collapsed
				self.Qsubgraph[iw] = collapsed
				for jw in self.Gw.neighbors(iw):
					if jw in resolvedAnchors:
						if iw < jw:
							if (iw,jw) in QccEdges:
								if (resolvedAnchors[iw],resolvedAnchors[jw]) not in self.Gw.edge[iw][jw]['qlist']:
									print "Resolution of edge",(iw,jw),"is configs",(resolvedAnchors[iw],resolvedAnchors[jw])
									print "Solution",why,"has a problem"
									print "Edge not here:",self.Gw.edge[iw][jw]['qlist']
								assert (resolvedAnchors[iw],resolvedAnchors[jw]) in self.Gw.edge[iw][jw]['qlist']
						else:
							if (jw,iw) in QccEdges:
								if (resolvedAnchors[jw],resolvedAnchors[iw]) not in self.Gw.edge[jw][iw]['qlist']:
									print "Resolution of edge",(jw,iw),"is configs",(resolvedAnchors[jw],resolvedAnchors[iw])
									print "Solution",why,"has a problem"
									print "Edge not here:",self.Gw.edge[jw][iw]['qlist']
								assert (resolvedAnchors[jw],resolvedAnchors[iw]) in self.Gw.edge[jw][iw]['qlist']
				totalCollapsed += 1
			else:
				resolutionFailures.add(iw)

		#some resolutions failed -- try doing a CSP only on these and their neighbors
		if len(resolutionFailures) > 0:
			print "FALLING BACK TO ALL_PAIRS RESOLUTION TECHNIQUE"
			failureWNodes = list(resolutionFailures)
			failureManifolds = dict([(iw,self.Gw.node[iw]['qccs'][self.QccSubgraph[iw]]) for iw in resolutionFailures])
			failureWNodes = set(failureWNodes)
			neighbors = ring(self.Gw,failureWNodes)
			neighborAnchors = ring(self.Gw,failureWNodes,2)
			neighborHubs = dict()
			print len(failureWNodes),"nodes could not be resolved"
			print "1-ring size",len(neighbors)
			print "2-ring size",len(neighborAnchors)
			print "Determining 1-ring hubs..."

			assignedWEdges = set()
			for iw,jw,d in self.Gw.edges_iter(data=True):
				assert iw < jw
				if not iw in self.QccSubgraph or not jw in self.QccSubgraph:
					continue
				imanifold = set(self.Gw.node[iw]['qccs'][self.QccSubgraph[iw]])
				jmanifold = set(self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]])
				for (iq,jq) in self.Gw.edge[iw][jw]['qlist']:
					if iq in imanifold and jq in jmanifold:
						assignedWEdges.add((iw,jw))
						assignedWEdges.add((jw,iw))
						break

			numManifoldConfigs = sum(len(v) for k,v in failureManifolds.iteritems())
			for jw in neighbors:
				if jw not in self.QccSubgraph:
					continue
				#might be missing a C-space edge from iw to jw?
				if jw not in resolvedAnchors:
					#jw might be isolated...
					continue
				jmanifold = self.Gw.node[jw]['qccs'][self.QccSubgraph[jw]]
				kwlist = [kw for kw in self.Gw.neighbors(jw) if kw in neighborAnchors and (jw,kw) in assignedWEdges]
				kwanchors = [resolvedAnchors.get(kw,None) for kw in kwlist]
				#determine potential hub configurations
				hubs = []
				for jq in jmanifold:
					valid = True
					for kw,anch in zip(kwlist,kwanchors):
						if anch != None:
							if not self.testAndConnect(jw,jq,kw,anch):
								valid = False
								break
					if valid:
						hubs.append(jq)
				if len(hubs) == 0:
					print "There appears to be a node in the 1 ring (%d) that has no hub vertex?"%(jw,)
					print "  Assignment is configuration",resolvedAnchors[jw]
					print "  What happened:",reasons[jw]
					print "  This configuration CANNOT connect to",[kw for kw,anch in zip(kwlist,kwanchors) if anch != None and not self.testAndConnect(jw,jq,kw,anch)]
					assert resolvedAnchors[jw] in jmanifold
					print "  Neighboring wnodes",kwlist
					print "  Neighboring anchors",kwanchors
					print "  Reasons for neighbors..."
					for kw in kwlist:
						print "   ",kw,reasons[kw]
					raise ValueError("Collapse function is mssed up")

				neighborHubs[jw] = hubs
				del self.Qsubgraph[jw]
			numHubConfigs = sum(len(hub) for hub in neighborHubs.itervalues())
			print "# of failed manifold configs total",numManifoldConfigs,"# of neighboring hub configs total",numHubConfigs
			#now do exhaustive connections between manifold configs in failureWNodes and hub configs in neighbors
			#add these options to the CSP
			nameToWnode = dict()
			csp = CSP()
			for iw,imanifold in failureManifolds.iteritems():
				csp.addVariable(list(imanifold),'q_'+str(iw))
				nameToWnode['q_'+str(iw)] = iw
			for jw,jhub in neighborHubs.iteritems():
				csp.addVariable(jhub,'qhub_'+str(jw))
				nameToWnode['qhub_'+str(jw)] = jw

			print "Determining all-pair connections..."
			numValidEdges = 0
			numTestedEdges = 0
			stepcount = 0
			for iw,imanifold in failureManifolds.iteritems():
				stepcount += len(self.Gw.neighbors(iw))
			c = ProgressUpdater(stepcount,1)
			if parallel:
				launched = 0
				validValues = {}
				for iw,imanifold in failureManifolds.iteritems():
					for jw in self.Gw.neighbors(iw):
						c.update()
						if jw in failureWNodes:
							for iq in imanifold:
								for jq in failureManifolds[jw]:
									numTestedEdges += 1
									tasks.put(('test',iw,iq,self.Gw.node[iw]['qlist'][iq],self.Gw.node[iw]['params'],jw,jq,self.Gw.node[jw]['qlist'][jq],self.Gw.node[jw]['params']))
									if launched < threads:
										launched += 1
										continue
									else:
										res = results.get()
										if res[0] == True:
											self.addEdge(res[1],res[2],res[3],res[4])
											try:
												validValues[(res[1],res[3])].append((res[2],res[4]))
											except KeyError:
												validValues[(res[1],res[3])] = [(res[2],res[4])]
											numValidEdges += 1
						elif jw in neighborHubs:
							for iq in imanifold:
								for jq in neighborHubs[jw]:
									numTestedEdges += 1
									tasks.put(('test',iw,iq,self.Gw.node[iw]['qlist'][iq],self.Gw.node[iw]['params'],jw,jq,self.Gw.node[jw]['qlist'][jq],self.Gw.node[jw]['params']))
									if launched < threads:
										launched += 1
										continue
									else:
										res = results.get()
										if res[0] == True:
											self.addEdge(res[1],res[2],res[3],res[4])
											try:
												validValues[(res[1],res[3])].append((res[2],res[4]))
											except KeyError:
												validValues[(res[1],res[3])] = [(res[2],res[4])]
											numValidEdges += 1
				for i in range(threads):
					res = results.get()
					if res[0] == True:
						self.addEdge(res[1],res[2],res[3],res[4])
						try:
							validValues[(res[1],res[3])].append((res[2],res[4]))
						except KeyError:
							validValues[(res[1],res[3])] = [(res[2],res[4])]
						numValidEdges += 1
				for key, values in validValues.iteritems():
					iw, jw = key
					if jw in failureWNodes:
						csp.addConstraint(values,("q_"+str(iw),"q_"+str(jw)))
					elif jw in neighborHubs:
						csp.addConstraint(values,("q_"+str(iw),"qhub_"+str(jw)))
			else:
				for iw,imanifold in failureManifolds.iteritems():
					for jw in self.Gw.neighbors(iw):
						c.update()
						if jw in failureWNodes:
							validValues = []
							for iq in imanifold:
								for jq in failureManifolds[jw]:
									numTestedEdges += 1
									if self.testAndConnect(iw,iq,jw,jq):
										validValues.append((iq,jq))
										numValidEdges += 1
							csp.addConstraint(validValues,("q_"+str(iw),"q_"+str(jw)))
						elif jw in neighborHubs:
							validValues = []
							for iq in imanifold:
								for jq in neighborHubs[jw]:
									numTestedEdges += 1
									if self.testAndConnect(iw,iq,jw,jq):
										validValues.append((iq,jq))
										numValidEdges += 1
							csp.addConstraint(validValues,("q_"+str(iw),"qhub_"+str(jw)))
			c.done()
			print "# of valid edges",numValidEdges,"# total",numTestedEdges
			#raw_input("Press enter to begin CSP assignment")
			print("Beginning CSP assignment")
			assignment = csp.heuristicMaxAssignment()
			assignment = csp.randomDescentAssignment(assignment)
			print "# of conflicts",csp.numConflicts(assignment)
			varsinconflict = []
			failuresinconflict = []
			neighborsinconflict = []
			for c in xrange(len(csp.constraints)):
				if not csp.testConstraint(c,assignment):
					for v in csp.conDomains[c]:
						assert isinstance(v,int)
						vname = csp.variables[v]
						assert isinstance(vname,str)
						varsinconflict.append(vname)
						if vname.startswith('qhub_'):
							neighborsinconflict.append(int(vname[5:]))
						else:
							failuresinconflict.append(int(vname[2:]))
			print "# of centers in conflict",len(failuresinconflict)
			print "# of neighbors in conflict",len(neighborsinconflict)
			for v,vname in enumerate(csp.variables):
				self.Qsubgraph[nameToWnode[vname]] = assignment[v]
			#sanity check
			for (k,v) in self.Qsubgraph.iteritems():
				assert v < len(self.Gw.node[k]['qlist'])
		else:
			print "NO NEED FOR FURTHER RESOLUTION"
		
		end = time.time()
		self.calculateResolutionFromSubgraph()	
		if len(self.QccSubgraph) > 0:
			return 0,0,totalCollapsed,triedToCollapse,end-start
		else:
			return -1,-1,-1,-1,end-start
			
	def decodeCSPAssignment(self,csp_assignment,visibility=False,parallel=False,threads=None,tasks=None,results=None):
		"""For an assignment produced by makeCSP, calculates self.Qsubgraph and the resolution."""
		if visibility:		
			self.QccSubgraph = dict()
		else:
			self.Qsubgraph = dict()
		vcount = 0
		for n,d in self.Gw.nodes_iter(data=True):
			if len(d['qlist']) > 0:
				assert vcount < len(csp_assignment)
				if csp_assignment[vcount] != None:
					if visibility:
						self.QccSubgraph[n] = csp_assignment[vcount]
					else:
						self.Qsubgraph[n] = csp_assignment[vcount]
				vcount += 1
		if visibility:
			print "Now collapsing QccSubgraph"
			return self.collapseQccSubgraph(parallel=parallel,threads=threads,tasks=tasks,results=results)
		else:
			self.calculateResolutionFromSubgraph()
			return None,None,None,None,0

	def encodeCSPAssignment(self):
		"""Given self.Qsubgraph, returns an assignment for the CSP produced by makeCSP."""
		if len(self.Qsubgraph) == 0:
			return None
		res = []
		for n,d in self.Gw.nodes_iter(data=True):
			if len(d['qlist'])==0: continue
			res.append(self.Qsubgraph.get(n,None))
		return res
		
	def saveNoncollapsibleManifolds(self, path, create=True):
		"""Call this in experimenting.py with path = 'experiments'.
		The saved .dot files can be visualized with GraphViz, i.e. the command
		'dot -Tpng manifold186.dot -o manifold186.png'"""
		raise NotImplementedError("Noncollapsible manifold save")
		if create and not os.path.isdir(path):
			os.makedirs(path)
		#loop over and identify noncollapsible representatives
		noncollapsible = set()
		manifolds = set()
		for iq in self.Qsubgraph:
			iw = self.Gq.node[iq]['windex']
			rep = self.Gw.node[iw]['qccs'].getSetRep(iq)
			if rep in manifolds:
				noncollapsible.add(rep)
			else:
				manifolds.add(rep)
		#for each noncollapsible representative, save subgraph as .dot file
		for irep in noncollapsible:
			windices = set()
			iw = self.Gq.node[irep]['windex']
			windices.add(iw)
			manifold = self.Gw.node[iw]['qccs'].getSetRef(irep)
			adjacentManifolds = set()
			for iq in manifold:
				assert self.Gq.node[iq]['windex'] == iw
				touched = [jq for jq in self.Gq.neighbors(iq) if self.Gq.node[jq]['windex'] != iw]
				for jq in touched:
					jw = self.Gq.node[jq]['windex']
					adjacentManifolds.add(self.Gw.node[jw]['qccs'].getSetRep(jq))
			nbunch = set()
			nbunch.update(manifold)
			for jrep in adjacentManifolds:
				jw = self.Gq.node[jrep]['windex']
				windices.add(jw)
				nbunch.update(self.Gw.node[jw]['qccs'].getSetRef(jrep))
			subgraph = nx.Graph(self.Gq.subgraph(nbunch))
			windices = list(windices)
			for n,d in subgraph.nodes_iter(data=True):
				d['style'] = 'filled'
				d['fillcolor'] = '/greys7/' + str((windices.index(d['windex']) % 7) + 1)
				if d['windex'] == iw: d['shape'] = 'diamond'
			nx.drawing.nx_pydot.write_dot(subgraph,path+'/manifold'+str(irep)+'.dot')
			
	def uniqueAssignment(self):
		"""Tries to solve the optimization max{S} cover(S) such that each point in Gw is assigned to
		at most one element of Gq, and cover(S) measures the number of edges of Gw covered.

		Each edge in Gw is assigned to at most one edge of Gq.  This is a constraint satisfaction
		optimization problem, where each edge is a variable, and they must meet the
		constraint that an adjacent edge's endpoint must match the endpoint of the edge.

		Specifically, for edge (u,v), need to meet the constraint that for all t in Gw.neighbors(u),
			Val(t) in Gq.neighbors(Val(u))
		and for all W in Gw.neighbors(V),
			Val(w) in Gq.neighbors(Val(v))

		I've done some experimenting with heuristics.
		A typical CSP heuristic is to assign the edge with the fewest remaining values
		(most-constrained-value) and breaks ties with largest degree (most-constraining-variable). 
		The value to which it is assigned is the least-constraining-value heuristic. The only issue
		is that this tries to find a contradiction along a branch as quickly as possible.
		If we want to maximize the number of edges assigned, then perhaps a different ordering would be
		better.  It performs OK.

		Least-constraining variable might be a good heuristic but it's relatively expensive
		least-constraining-value works well for a value assignment.
		"""
		print "Trying CSP solver assignment..."
		t0 = time.time()
		csp = self.makeCSP()
		print "Time to make CSP:",time.time()-t0
		assignment = csp.heuristicMaxAssignment()
		if assignment != None:
			self.decodeCSPAssignment(assignment)
		else:
			print "CSP solver failed"

	def randomDescentAssignment(self,randomize=True):
		"""Tries to solve the optimization max{S} cover(S) such that each point in Gw is assigned to
		at most one element of Gq, and cover(S) measures the number of edges of Gw covered.

		Each node in Gw is assigned to at most one node of Gq.  This is a constraint satisfaction
		optimization problem, where each node is a variable, and they must meet the
		constraint that the assignments to two adjacent nodes in Gw correspond to an edge in Gq

		Specifically, for node i, need to meet the constraint that for all j in Gw.neighbors(i),
			The edge (Val(i),Val(j)) exists in in Gq

		This heuristic samples a random assignment for each node (if randomize=True).  Then, for all nodes
		that conflict with their neighbors, it attempts to switch the assignment to a less-
		conflicting one.  This repeats until no progress can be made.
		"""
		assignment = (self.encodeCSPAssignment() if not randomize else None)
		t0 = time.time()
		csp = self.makeCSP()
		print "Time to make CSP:",time.time()-t0
		assignment = csp.randomDescentAssignment(assignment,perturb=True)
		if assignment != None:
			self.decodeCSPAssignment(assignment)
		else:
			print "CSP solver failed"
		return

	def pointwiseAssignment(self,NConfigsPerPoint=10):
		"""Moves pointwise through the workspace grid and generates a configuration by
		local optimization from the previous grid point.  If a configuration cannot be
		generated, then the edge is marked as infeasible and a random configuration is
		generated.
		"""
		print "Generating pointwise assignment"
		start_node_count = self.numConfigurations()
		qmap = dict()
		emap = list()
		efailures = 0
		c = ProgressUpdater(self.Gw.number_of_nodes(),5)
		for i,d in self.Gw.nodes_iter(data=True):
			c.update()
			adjacentWithConfigs = [j for j in self.Gw.neighbors(i) if j < i and j in qmap]
			solved = False
			if len(adjacentWithConfigs) > 0:
				self.resolution.setIKProblem(d['params'])
				self.resolution.ikSolverParams.startRandom = False
				qstart = robot_average(self.robot,[qmap[j] for j in adjacentWithConfigs])
				self.robot.setConfig(qstart)
				q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
				if q is not None:
					for j in adjacentWithConfigs:
						if self.resolution.validEdge(qmap[j],q,self.Gw.node[j]['params'],d['params']):
							emap.append((i,j))
							qmap[i] = q
							solved = True
						else:
							efailures += 1
			if not solved:
				efailures += len([j for j in self.Gw.neighbors(i) if j < i])
				#no edge between neighbors and this node
				self.resolution.setIKProblem(d['params'])
				self.resolution.ikSolverParams.startRandom = True
				for sample in xrange(NConfigsPerPoint):
					#solver will sample randomly
					q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
					if q is not None:
						qmap[i] = q
						break
		c.done()
		print "Sampled",len(qmap),"of possible",self.Gw.number_of_nodes()
		print "Connected",len(emap),"of possible",self.Gw.number_of_edges(),"failed to connect",efailures
		nodemap = dict()
		for (i,q) in qmap.iteritems():
			iq = self.addNode(i,q)
			nodemap[i] = iq
		for (i,j) in emap:
			iq,jq = nodemap[i],nodemap[j]
			self.addEdge(i,iq,j,jq)
		self.Qsubgraph = nodemap
		self.calculateResolutionFromSubgraph()

	def optimize(self,numIters=10):
		"""Optimizes the configurations of the current resolution to minimize
		joint-space path lengths using coordinate descent."""
		Gresolution = self.resolution.getResolutionGraph()
		for iters in xrange(numIters):
			print "Iteration",iters
			sumdistances = 0
			c = ProgressUpdater(Gresolution.number_of_nodes(),10)
			for v in Gresolution.nodes_iter():
				c.update()
				if len(Gresolution.neighbors(v)) == 0: 
					continue
				x = self.resolution.Gw.node[v]['params']
				q0 = self.resolution.Gw.node[v]['config']
				wneighbors = [w for w in Gresolution.neighbors(v)]
				qneighbors = [self.resolution.Gw.node[w]['config'] for w in Gresolution.neighbors(v)]
				xneighbors = [self.resolution.Gw.node[i]['params'] for i in wneighbors]
				d0 = sum(self.robot.distance(q0,qw) for qw in qneighbors)
				qavg = robot_average(self.robot,qneighbors)
				#try to move toward average
				self.resolution.setIKProblem(x)
				self.resolution.ikSolverParams.startRandom = False
				maxTries = 10
				moved = False
				for tries in xrange(maxTries):
					self.robot.setConfig(qavg)
					q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
					if q is not None:
						#check all edges from neighbors
						d = sum(self.robot.distance(q,qw) for qw in qneighbors)
						valid = True
						if d > d0:
							valid = False
						else:
							for qw,xw in zip(qneighbors,xneighbors):
								if not self.resolution.validEdge(q,qw,x,xw):
									valid = False
									break
						if valid:
							assert q is not None
							self.resolution.Gw.node[v]['config'] = q
							sumdistances += d
							moved = True
							break
					#not moved, try subdividing
					qavg = self.robot.interpolate(q0,qavg,0.5)
				if not moved:
					#print "Unable to move configuration",v
					sumdistances += d0
			c.done()
			print "Changed average path length to",sumdistances/Gresolution.number_of_edges()
		#now add new q's and edges between to the graph
		self.fromResolutionToGraph()

	def sampleConflicts(self,NConfigsPerPoint=10):
		conflicts = set()
		for i,j,d in self.Gw.edges_iter(data=True):
			if self.resolution.isResolvedEdge(i,j):
				continue
			if not self.resolution.isResolvedNode(i) or not self.resolution.isResolvedNode(j):
				continue
			conflicts.add(i)
			conflicts.add(j)
		added = []
		for iw in conflicts:
			d = self.Gw.node[iw]
			self.resolution.setIKProblem(d['params'])
			self.resolution.ikSolverParams.startRandom = True
			for sample in xrange(NConfigsPerPoint):
				#solver will sample randomly
				q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
				if q is not None:
					iq = self.addNode(iw,q)
					added.append((iw,iq))
			#try local solving from neighbors too
			self.resolution.ikSolverParams.startRandom = True
			Qseed = sum([self.Gw.node[jw]['qlist'] for jw in self.Gw.neighbors(iw) if jw in conflicts],[])
			if len(Qseed) > 0:
				for sample in xrange(NConfigsPerPoint):
					qseed = random.choice(Qseed)
					self.robot.setConfig(qseed)
					q = self.resolution.ikTemplate.solve(self.robot,self.resolution.ikSolverParams)
					if q is not None:
						iq = self.addNode(iw,q)
						added.append((iw,iq))
		print len(added),"configs added for",len(conflicts),"workspace nodes"
		numEdgesAttempted = 0
		numEdgesAdded = 0
		c = ProgressUpdater(len(added),5)
		print "Connecting edges..."
		for (iw,iq) in added:
			c.update()
			for jw in self.Gw.neighbors(iw):
				for jq in xrange(len(self.Gw.node[jw]['qlist'])):
					if (jw,jq) in added and jw > iw: continue
					numEdgesAttempted += 1
					if self.testAndConnect(iw,iq,jw,jq):
						numEdgesAdded += 1
		c.done()
		print "Attempted",numEdgesAttempted,"edge connections and added",numEdgesAdded

	def solveCSP(self,visibility=False,parallel=False,numRandomDescents=2,numInitialGuesses=10,threads=None,tasks=None,results=None):
		csp = self.makeCSP(visibility=visibility)
		start = time.time()
		#assignment = csp.gecodeAssignment(maximize=True)
		#assignment = csp.sugarMaxAssignment()
		#assignment = csp.backtrackingAssignment()
		#if assignment == None:						
		assignment = csp.heuristicMaxAssignment()
		for i in xrange(numRandomDescents):
			assignment = csp.randomDescentAssignment(assignment)
		bestMethod = 'heuristic max'
		bestAssignment = assignment
		bestConflicts = csp.numConflicts(assignment)
		for guess in xrange(numInitialGuesses-1):
			assignment = csp.randomDescentAssignment()
			for i in xrange(numRandomDescents):
				assignment = csp.randomDescentAssignment(assignment)
			conflicts = csp.numConflicts(assignment)
			if conflicts < bestConflicts:
				bestMethod = 'random'
				bestAssignment = assignment
				bestConflicts = conflicts
		end = time.time()
		solve_time = end - start
		print "CSP solve time",solve_time,"minimum # of conflicts",bestConflicts,"using method",bestMethod
		collapse_time = self.decodeCSPAssignment(bestAssignment,visibility=visibility,parallel=parallel,threads=threads,tasks=tasks,results=results)[4]
		if visibility:
			print "QCC CSP resolve time",collapse_time
		return bestConflicts,bestMethod,solve_time,collapse_time

	def select(self,x,nodeRadius,edgeRadius,cc_filter=True,use_resolution='auto'):
		if use_resolution == 'auto':
			if len(self.Qsubgraph) == 0:
				use_resolution = False
		demin = float('inf')
		emin = None
		for aw,bw,data in self.Gw.edges_iter(data=True):
			if len(data['qlist']) == 0: 
				continue
			a = self.Gw.node[aw]['params']
			b = self.Gw.node[bw]['params']
			(d,u) = segment_point_distance(a,b,x)
			if d < demin:
				demin = d
				emin = (aw,bw,u)
		if demin > edgeRadius:
			#print "Closest edge:",demin
			emin = None
		if emin:
			(aw,bw,u) = emin
			if use_resolution:
				ares = [aq for aq in range(len(self.Gw.node[aw]['qlist'])) if aq in self.Qsubgraph]
				bres = [bq for bq in range(len(self.Gw.node[bw]['qlist'])) if bq in self.Qsubgraph]
				res = []
				for aq in ares:
					for bq in bres:
						if self.Gq.has_edge(aq,bq) or self.Gq.has_edge(bq,aq):
							res.append(self.robot.interpolate(self.Gq.node[aq]['config'],self.Gq.node[bq]['config'],u))
						else:
							#whatever, do it anyway
							#res.append(self.robot.interpolate(self.Gq.node[aq]['config'],self.Gq.node[bq]['config'],u))
							pass
				if len(res) == 0:
					print "Nothing in Qsubgraph???"
				return res
			elif cc_filter:
				#get connected component
				configs = set(aq for (aq,bq) in self.Gw.edge[aw][bw]['qlist']) | set(bq for (aq,bq) in self.Gw.edge[aw][bw]['qlist'])
				Gsub = self.Gq.subgraph(configs)
				ccs = list(nx.connected_components(Gsub))
				touched = set()
				res = []
				for (aq,bq) in self.Gw.edge[aw][bw]['qlist']:
					aindex = [i for i,cc in enumerate(ccs) if aq in cc]
					bindex = [i for i,cc in enumerate(ccs) if bq in cc]
					assert len(aindex)==1 and len(bindex)==1
					assert aindex[0] == bindex[0]
					if aindex[0] in touched: continue
					touched.add(aindex[0])
					res.append(self.robot.interpolate(self.Gq.node[aq]['config'],self.Gq.node[bq]['config'],u))
				return res
			else:
				#return all
				res = []
				for (aq,bq) in self.Gw.edge[aw][bw]['qlist']:
					res.append(self.robot.interpolate(self.Gq.node[aq]['config'],self.Gq.node[bq]['config'],u))
				return res
		dnmin = float('inf')
		nmin = None
		for (iw,data) in self.Gw.nodes_iter(data=True):
			if len(data['qlist']) == 0: 
				continue
			d = workspace_distance(x,data['params']) 
			if d < dnmin:
				nmin = iw
				dnmin = d
		if dnmin > nodeRadius:
			#print "Node distance",dnmin
			nmin = None
		if nmin:
			if use_resolution:
				qs = [qi for iq,qi in enumerate(self.Gw.node[nmin]['qlist']) if iq in self.Qsubgraph]
				return qs
			elif cc_filter:
				#get connected component
				configs = set(range(len(self.Gw.node[nmin]['qlist'])))
				for n in self.Gw.neighbors(nmin):
					configs |= set(range(len(self.Gw.node[n]['qlist'])))
				Gsub = self.Gq.subgraph(configs)
				ccs = list(nx.connected_components(Gsub))
				res = []
				touched = set()
				for v in range(len(self.Gw.node[nmin]['qlist'])):
					vindex = [i for i,cc in enumerate(ccs) if v in cc]
					assert len(vindex) == 1
					if vindex[0] in touched: continue
					touched.add(vindex[0])
					res.append(self.Gq.node[v]['config'])
				return res
			else:
				return self.Gw.node[nmin]['qlist']
		return None

	def calculateResolutionFromSubgraph(self):
		"""Computes the resolution member from the Qsubgraph set."""
		self.resolution.clearResolution()
		for iw,iq in self.Qsubgraph.iteritems():
			self.resolution.setConfig(iw, self.Gw.node[iw]['qlist'][iq])
			for jw in self.Gw.neighbors(iw):
				if jw > iw: continue
				if self.resolution.isResolvedNode(jw) and self.hasEdge(iw,iq,jw,self.Qsubgraph[jw]):
					self.resolution.markConnected(iw,jw)
		#self.sanityCheck()

	def fromResolutionToGraph(self,clear=False):
		"""Adds all resolved configs/edges in the resolution member to Gq, and
		updates Qsubgraph.  If clear = True, the graph is cleared first."""
		if clear:
			self.initWorkspaceGraph()
		iqmap = dict()
		for i,d in self.resolution.Gw.nodes_iter(data=True):
			if d.get('config',None) is not None:
				iqmap[i] = self.addNode(i,d['config'])
		for i,j,d in self.resolution.Gw.edges_iter(data=True):
			if d.get('connected',False):
				iq = iqmap[i]
				jq = iqmap[j]
				self.addEdge(i,iq,j,jq)
		self.Qsubgraph = iqmap
		#self.sanityCheck()
		
class SequentialResolver:
	def __init__(self, problem, rr, Nw, method, visibility):
		self.stats = StatCollector(problem=problem)
		self.rr = rr
		try:
			self.stats.setCacheSize(self.rr.Gw.db.getSize())
		except:
			pass
		self.Nw = Nw
		self.method = method	
		self.visibility = visibility
		
	def sampleAndConnect(self, Nq=100, k=None):
		self.stats.setNw(self.Nw)
		self.stats.setMethod(self.method)
		self.stats.setTotalWorkspaceEdges(self.rr.Gw.number_of_edges())
		self.stats.setNq(Nq)
		try:
			self.stats.setCacheSize(self.rr.Gw.db.getSize())
		except:
			pass
		self.stats.setVisK(k)
		t_sample, t_connect = self.rr.sampleConfigurationSpace(NConfigsPerPoint=Nq,connect=True,visibility=self.visibility,k=k)
		self.stats.setSelfMotionCcsStats(self.rr.getSelfMotionCcsStats())
		self.stats.setSamplingTime(t_sample)
		self.stats.setEdgeTime(t_connect)
		
	def solveCSP(self):
		min_conflicts, csp_method, solve_time, collapse_time = self.rr.solveCSP(visibility=self.visibility, parallel=False)
		self.stats.setMinConflicts(min_conflicts)
		self.stats.setCSPMethod(csp_method)
		self.stats.setCSPSolveTime(solve_time)
		self.stats.setCollapseTime(collapse_time)
		
	def writeReport(self):
		self.stats.writeReport()
		
	def run(self, Nq, visK, folder, solveCSP=True):
		self.sampleAndConnect(Nq, visK)
		self.rr.save(folder,close=not solveCSP)
		if solveCSP:
			self.solveCSP()
			self.rr.save(folder,close=True)
		self.writeReport()

class ParallelResolver:
	def __init__(self, rr, makeParams=None, threads=4):
		try:
			self.stats = StatCollector(problem=makeParams['problem'])
		except:
			self.stats = StatCollector()
		self.rr = rr
		try:
			self.stats.setCacheSize(rr.Gw.db.getSize())
		except:
			pass
		self.makeParams = makeParams
		self.setThreads(threads)
		self.tasks = multiprocessing.Queue()
		self.results = multiprocessing.Queue()
		
	def setThreads(self, threads):
		self.threads = threads
		self.stats.setThreads(threads)
		
	def getStats(self):
		return self.stats
		
	def spawnWorkers(self):
		"""Assumes makeParams != None"""
		for i in range(self.threads):
			multiprocessing.Process(target=worker, args=(i, self.makeParams, self.tasks, self.results)).start()
		
	def superviseQSampling(self, Nq=100, k=None):		
		#sample configs and construct self-motion manifolds
		print "Sampling configurations and constructing self-motion connected components"
		self.stats.setNw(self.makeParams['Nw'])
		self.stats.setMethod(self.makeParams['workspace_graph'])
		self.stats.setTotalWorkspaceEdges(self.rr.Gw.number_of_edges())
		self.stats.setNq(Nq)
		try:
			self.stats.setCacheSize(self.rr.Gw.db.getSize())
		except:
			pass
		self.stats.setVisK(k)
		c = ProgressUpdater(self.rr.Gw.number_of_nodes(),5)
		launched = 0
		start = time.time()
		for i,d in self.rr.Gw.nodes(data=True):
			self.tasks.put(('sample',d['params'],Nq,k,i))
			if launched < self.threads:
				launched += 1
				continue
			else:
				configs, qccs, qcc_parents, uid = self.results.get()
				self.rr.Gw.node[uid]['qlist'] = configs
				self.rr.Gw.node[uid]['qccs'] = qccs
				self.rr.Gw.node[uid]['qcc_parents'] = qcc_parents
				self.stats.add(qccs)
				c.update()
		for i in range(self.threads):
			configs, qccs, qcc_parents, uid = self.results.get()
			self.rr.Gw.node[uid]['qlist'] = configs
			self.rr.Gw.node[uid]['qccs'] = qccs
			self.rr.Gw.node[uid]['qcc_parents'] = qcc_parents
			self.stats.add(qccs)
			c.update()
		end = time.time()
		self.stats.setSamplingTime(end - start)
		c.done()
		
	def superviseEdgeTests(self):
		#connect edges
		print "Testing configuration space edges"
		self.stats.setTotalWorkspaceEdges(self.rr.Gw.number_of_edges())
		try:
			self.stats.setCacheSize(self.rr.Gw.db.getSize())
		except:
			pass
		c = ProgressUpdater(self.rr.Gw.number_of_edges(),5)
		launched = 0
		start = time.time()
		for (iw,jw) in self.rr.Gw.edges_iter():
			#need iw < jw for proper bookkeeping
			if iw > jw:
				zw = iw
				iw = jw
				jw = zw
			di = self.rr.Gw.node[iw]
			dj = self.rr.Gw.node[jw]
			self.tasks.put(('connect',iw,di['qlist'],di['qccs'],di['params'],jw,dj['qlist'],dj['qccs'],dj['params']))
			if launched < self.threads:
				launched += 1
				continue
			else:
				result = self.results.get()
				self.rr.Gw.edge[result[1]][result[2]]['qlist'] = result[0]
				if len(result[0]) > 0:
					self.stats.inc()
				c.update()
		for i in range(self.threads):
			result = self.results.get()
			self.rr.Gw.edge[result[1]][result[2]]['qlist'] = result[0]
			if len(result[0]) > 0:
				self.stats.inc()
			c.update()
		end = time.time()
		self.stats.setEdgeTime(end - start)
		c.done()
		
	def solveCSP(self):
		min_conflicts, csp_method, solve_time, collapse_time = self.rr.solveCSP(visibility=True, parallel=True, threads=self.threads, tasks=self.tasks, results=self.results)
		self.stats.setMinConflicts(min_conflicts)
		self.stats.setCSPMethod(csp_method)
		self.stats.setCSPSolveTime(solve_time)
		self.stats.setCollapseTime(collapse_time)

	def writeReport(self):
		self.stats.setProblem(self.makeParams['problem'])
		self.stats.writeReport()
		
	def killWorkers(self):
		#terminate processes
		for p in multiprocessing.active_children():
			p.terminate()
			
	def run(self, Nq, visK, folder, solveCSP=True):
		self.spawnWorkers()
		self.superviseQSampling(Nq,visK)
		self.rr.save(folder)
		self.superviseEdgeTests()
		self.rr.save(folder,close=not solveCSP)
		if solveCSP:
			self.solveCSP()
			self.rr.save(folder,close=True)
		self.writeReport()
		self.killWorkers()

				
def worker(index, pdef, tasks, results):
	pdef['worker'] = index
	pdef['clean'] = True
	pdef['Nw'] = 0
	filename = pdef['filename']
	
	world = WorldModel()
	world.readFile(filename)
	resolution = RedundancyResolutionGraph(world)
	rr = RedundancyResolver(resolution,cache=0)
	resolution.readJsonSetup(pdef)

	while True:
		task = tasks.get()
		if task[0] == 'sample':
			configs, qccs, qcc_parents = rr.parallelQSamplingAtP(task[1], NConfigs=task[2], k=task[3])
			results.put((configs, qccs, qcc_parents, task[4]))
		elif task[0] == 'connect':
			qlist = rr.parallelEdgeConnection(task[2], task[3], task[4], task[6], task[7], task[8])
			results.put((qlist, task[1], task[5]))
		elif task[0] == 'test':
			validity = resolution.validEdge(task[3],task[7],task[4],task[8])
			results.put((validity, task[1], task[2], task[5], task[6]))
		elif task[0] == 'kill':
			return

