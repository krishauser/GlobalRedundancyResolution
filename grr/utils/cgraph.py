import networkx as nx
from cache import CachedDB

class CGraph:
	"""A networkx.Graph-like interface where the data associated with graph nodes and edges are
	stored in a persistent database.  The graph is assumed undirected.

	All node/edge data must be pickle-able.

	One small difference between networkx and this is that edges (i,j) are sorted in order
	of increasing id.
	"""
	class DBCaller:
		def __init__(self, db, key1, twokeys=True):
			self.db = db
			self.key1 = key1
			self.twokeys = twokeys
		def __getitem__(self, key2):
			if self.key1 == None:
				return self.db[key2] if not self.twokeys else CGraph.DBCaller(self.db, key2)
			else:
				key = (self.key1, key2) if self.key1 < key2 else (key2, self.key1)
				return self.db[key]
			
	def __init__(self, cachesize=10000, db_loc="shelf", flush=True):
		"""Initializer. 

		Arguments:
		- cachesize: the maximum number of entries to store in-memory from the database.
		- db_loc: the location of the database's backing store.
		- flush: if True (default), clears the database first.  Otherwise, the old data is
		  loaded from disk.  Behavior is undefined when flush = True?
		"""
		self.G = nx.Graph(qnodes=0)
		self.db = CachedDB(maxsize=cachesize, db_loc=db_loc, flush=flush)
		self.node = CGraph.DBCaller(self.db, None, twokeys=False)
		self.edge = CGraph.DBCaller(self.db, None, twokeys=True)

	def load(self,fn):
		import os
		import shelve
		self.G = nx.read_gpickle(fn)
		fn = os.path.splitext(fn)[0]+'.shelf'
		self.db.table = {}
		self.db.qhead = None
		self.db.qtail = None
		self.db.numDBGets = 0
		self.db.numDBPuts = 0
		self.db.fn = fn
		self.db.db.close()
		self.db.db = shelve.open(fn,'w')
		print "Shelf load opened",len(self.db.db),"entries",

	def save(self,fn,close=False):
		import os
		import shutil
		import shelve
		nx.write_gpickle(self.G,fn)
		fn = os.path.splitext(fn)[0]+'.shelf'
		self.db.dump()
		self.db.db.close()
		if fn != self.db.fn:
			print "Copying shelf file to",fn,"..."
			shutil.copy(self.db.fn,fn)
		if not close:
			self.db.db = shelve.open(self.db.fn)
		
	def add_node(self, n, **attr):
		"""Adds node to graph if does not exist. Updates its attributes otherwise."""
		self.G.add_node(n)
		self.db[n] = attr
		
	def add_edge(self, u, v, **attr):
		"""Add edge between u and v. Add u and v to graph if necessary."""
		if u not in self.G.node:
			self.add_node(u)
		if v not in self.G.node:
			self.add_node(v)
		self.G.add_edge(u, v)
		uid = (u,v) if u < v else (v,u)
		self.db[uid] = attr
		
	def nodes(self, data=False):
		if data:
			return [(i, self.db[i]) for i in self.G.nodes()]
		else:
			return self.G.nodes()

	def nodes_iter(self, data=False):
		if data:
			for i in self.G.nodes_iter():
				yield (i, self.db[i])
		else:
			for i in self.G.nodes_iter():
				yield i
	
	def edges(self, data=False):
		if data:
			return [(i, j, self.db[(i, j)]) if i < j else (j, i, self.db[(j, i)]) for i, j in self.G.edges()]
		else:
			return [(i, j) if i < j else (j, i) for i, j in self.G.edges()]
				
	def edges_iter(self, data=False):
		if data:
			for i,j in self.G.edges():
				if i < j:
					yield i, j, self.db[(i, j)]
				else:
					yield j, i, self.db[(j, i)]
		else:
			for i,j in self.G.edges_iter():
				if i < j:
					yield i, j
				else:
					yield j, i
					
	def neighbors(self, n):
		return self.G.neighbors(n)

	def neighbors_iter(self, n):
		return self.G.neighbors_iter(n)
		
	def has_node(self,n):
		return self.G.has_node(n)

	def has_edge(self,i,j):
		return self.G.has_edge(i,j)

	def number_of_nodes(self):
		return self.G.number_of_nodes()
		
	def number_of_edges(self):
		return self.G.number_of_edges()
