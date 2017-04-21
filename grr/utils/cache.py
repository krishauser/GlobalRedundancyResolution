import os
from shutil import rmtree
import time
import shelve
import random
import weakref

class CachedDB:
	"""A database with in-memory cache and backing storage on disk.

	Retrieve items from the database using the standard bracket []
	operator (item = db[key]).  The items' properties are not changed
	until the item is assigned back to the database (db[key] = item).
	For compatibility with the shelve module, keys are cast to str's
	first (so "1" and 1 refer to the same key).

	Uses a write-back strategy (i.e., items are written to disk
	whenever they are updated and need to be dumped).  To force a dump
	to disk, call dump() on this object.

	Attributes:
	- maxsize: maximum size of the database
	- printFrequency: if you wish to have some printed DB updates, set this
	  to some positive integer.  Every printFrequency DB operations, a
	  printout will occur.
	- table: resident elements
	- db: backing store, a "shelf"
	- qhead/qtail: the LRU queue.
	- numDBGets/numDBPuts: counts the number of gets/puts to the database
	"""
	def __init__(self, maxsize=10000, db_loc="shelf", flush=True):
		"""Initializes the database.

		Arguments:
		- maxsize: maximum number of items to store resident in memory.
		- db_loc: the location of the backing store.
		- flush: if True, creates an empty database.  If False, loads
		  the previous database from disk.
		"""
		assert maxsize > 1
		self.maxsize = maxsize
		self.verbose = False
		self.printFrequency = None
		self.table = {}
		self.qhead = None
		self.qtail = None
		self.numDBGets = 0
		self.numDBPuts = 0
		
		try:
			if flush and os.path.isdir(db_loc):
				rmtree(db_loc)
			if not os.path.isdir(db_loc):
				os.mkdir(db_loc)
			self.fn = os.path.join(db_loc,"shelf")
			self.db = shelve.open(self.fn)
		except Exception as e:
			import traceback
			traceback.print_exc()
			raise ValueError("Error opening shelf file "+self.fn,", is it closed for editing?")

		
	def getSize(self):
		"""Returns maximum size of cache."""
		return self.maxsize

	def dump(self):
		"""Dumps all changes to disk."""
		print "Dumping Cache to",self.fn
		qnode = self.qtail
		numDumped = 0
		numIgnored = 0
		while qnode:
			if qnode.updated():
				numDumped += 1
				self.db[qnode.getID()] = qnode.getProps()
				qnode.recorded()
			else:
				numIgnored += 1
			qnode = qnode.getNext()
		#may need to do this to flush
		self.db.close()
		print "Dumped",numDumped,"entries, ignored",numIgnored
		self.db = shelve.open(self.fn,'w')

	def _strconv(self, uid):
		try:
			return str(uid)
		except:
			assert False, "uid must be convertible to string"
		
	def _enqueue(self, uid, props={}):
		uid = self._strconv(uid)
		#TODO assert uid is not in database?
		if uid in self.table:
			assert False, "item already exists"
			
		#Dequeue
		if len(self.table) >= self.maxsize:
			if self.qhead.updated():
				self.numDBPuts += 1
				if self.verbose and (self.numDBGets + self.numDBPuts)%self.printFrequency == 0:
					print "Cache: Putting %s to database %s. Overall %d gets, %d puts"%(self.qhead.getID(),self.fn,self.numDBGets,self.numDBPuts)
				self.db[self.qhead.getID()] = self.qhead.getProps()
				self.qhead.recorded()
			del self.table[self.qhead.getID()]
			self.qhead = self.qhead.getPrev()
			self.qhead.getNext().setPrev(None)
			self.qhead.setNext(None)
		
		#Enqueue
		entry = CacheEntry(uid, props)
		entry.update()
		#TODO weak ref?
		self.table[entry.getID()] = entry
		if self.qtail != None:
			self.qtail.setPrev(entry)
		entry.setNext(self.qtail)
		if self.qhead == None:
			self.qhead = entry
		self.qtail = entry
		return entry
	
	def _bump(self, uid):
		uid = self._strconv(uid)
		assert uid in self.table
		entry = self.table[uid]
		if entry.getPrev() == None:
			return entry
		elif entry.getNext() == None:
			self.qhead = entry.getPrev()
		else:
			entry.getNext().setPrev(entry.getPrev())
		entry.getPrev().setNext(entry.getNext())
		entry.setPrev(None)
		self.qtail.setPrev(entry)
		entry.setNext(self.qtail)
		self.qtail = entry
		return entry
		
	def __contains__(self, uid):
		uid = self._strconv(uid)
		return uid in self.table or uid in self.db
		
	def __getitem__(self, uid):
		uid = self._strconv(uid)
		
		#if in cache, bump and return it
		if uid in self.table:
			return self._bump(uid)
			
		#elif in database, cache and return it
		self.numDBGets += 1
		if self.verbose and (self.numDBGets + self.numDBPuts)%self.printFrequency == 0:
			print "Cache: Getting %s from database %s. Overall %d gets, %d puts"%(uid,self.fn,self.numDBGets,self.numDBPuts)
		props = self.db[uid] #raises KeyError
		return self._enqueue(uid, props)
		
	def __setitem__(self, uid, props):
		"""Updates props if uid exists. Creates new entry otherwise."""
		uid = self._strconv(uid)
			
		#in cache
		if uid in self.table:
			entry = self._bump(uid)
			entry.setProps(props)
			
		#in db or nonexistent
		else:
			self.db[uid] = props
			self._enqueue(uid, props)
		
class CacheEntry:
	def __init__(self, uid, props={}):
		self.props = props
		self.uid = uid
		self.modified = False
		self.prev = None
		self.next = None
		
	def __setitem__(self, key, item):
		self.update()
		self.props[key] = item
		
	def __getitem__(self, key):
		item = self.props[key]
		if type(item) == list or type(item) == dict or type(item) == set:
			#in case user mutates the property
			self.update()
		return item
		
	def __contains__(self, key):
		return key in self.props
		
	def setProps(self, props):
		self.update()
		self.props = props
		
	def setPrev(self, prev):
		self.prev = (weakref.proxy(prev) if (prev is not None and not isinstance(prev,weakref.ProxyType)) else prev)
		
	def setNext(self, nxt):
		self.next = nxt
		
	def getPrev(self):
		return self.prev
		
	def getNext(self):
		return self.next
		
	def getProps(self):
		return self.props
		
	def getID(self):
		return self.uid
		
	def recorded(self):
		self.modified = False
		
	def update(self):
		self.modified = True
		
	def updated(self):
		return self.modified

if __name__ == "__main__":
	#Testing
	try:
		N = 100000
		cachesize = 10000
		cache = CachedDB(maxsize=cachesize)
		store_start = time.clock()
		for i in range(N):
			cache[(i,i)] = {'joints' : [float(i)]*30}
		store_end = time.clock()
		indices = [random.randint(0,N-1) for i in range(cachesize)]
		get_from_db_start = time.clock()
		count = 0
		for i in indices:
			count += 1
			entry = cache[(i,i)]
			if count % 1000 == 0:
				print(entry.getID())
				print(entry['joints'])
		get_from_db_end = time.clock()
		get_from_cache_start = time.clock()
		count = 0
		for i in indices:
			count += 1
			entry = cache[(i,i)]
			if count % 1000 == 0:
				print(entry.getID())
				print(entry['joints'])
		get_from_cache_end = time.clock()
		print("time to store: " + str(store_end - store_start))
		print("time to get from db: " + str(get_from_db_end - get_from_db_start))
		print("time to get from cache: " + str(get_from_cache_end - get_from_cache_start))
	finally:
		cache.db.close()
