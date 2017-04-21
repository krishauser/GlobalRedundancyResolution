class RepSet:
        """Auxiliary class used by DisjointSet class.
        Lets one associate a representative element with a set object."""
	def __init__(self, rep):
		self.rep = rep
		self.s = set()

	def __len__(self):
		return len(self.s)

	def add(self,elem):
		self.s.add(elem)

	def discard(self,elem):
		self.s.discard(elem)

	def remove(self,elem):
		self.s.remove(elem)

	def getRep(self):
		return self.rep

	def getSet(self):
		return self.s

	def setRep(self,rep):
		self.rep = rep

class DisjointSet:
        """DisjointSet implementation that lets one get the internal subsets
        associated with any member of the larger DisjointSet object in addition
        to querying in constant time whether or not two elements of the
        DisjointSet object are in the same internal subset."""
	def __init__(self):
		self.map = {}
		self.reps = set()
		
	def __len__(self):
                """Return the number of internal subsets."""
		return len(self.reps)

	def __iter__(self):
                """Iterate over the representative of each of the internal
                subsets."""
		return iter(self.reps)
		
	def add(self,elem,join=None): #O(1)
                """Add an element to the DisjointSet object. Optionally add it
                to the internal subset of which the 'join' parameter is an
                element."""
		if elem in self.map:
			return
		elif join == None:
			self.map[elem] = RepSet(elem)
			self.map[elem].add(elem)
			self.reps.add(elem)
		else:
			try:
				#this creates reference equality
				self.map[elem] = self.map[join]
				self.map[join].add(elem)
			except KeyError:
				print("Error: join must be in a set")
				assert False

	def remove(self,elem): #O(1)
                """Remove an element from the DisjointSet object."""
		try:
			s = self.map[elem]
			s.remove(elem)
			self.map.pop(elem)
			if s.getRep() == elem:
				self.reps.remove(elem)
				for e in s.getSet():
					self.reps.add(e)
					s.setRep(e)
					break
		except KeyError:
			assert False,"Error: elem must be in a set"

	def getReps(self):
                """Get a set containing the representative of each of the
                internal subsets."""
		return self.reps

	def getSetRef(self,elem):
                """Get a reference to the set of which elem is an element."""
		return self.map[elem].getSet()

	def getSetCopy(self,elem):
                """Get a copy of the set of which elem is an element."""
		return self.map[elem].getSet().copy()

	def getSetRep(self,elem):
                """Get the represenative of the set of which elem is an
                element."""
		return self.map[elem].getRep()
		
	def length(self,elem):
                """Return the cardinality of the set of which elem is an
                element."""
		return len(self.getSetRef(elem))

	def iterate(self,elem):
                """Iterate over the elements of the set of which elem is an
                element."""
		return (i for i in self.getSetRef(elem))

	def same(self,elemA,elemB): #O(1)
                """Returns True if elemA is in the same internal subset as
                elemB. Otherwise returns False."""
		return self.map[elemA] is self.map[elemB]

	def merge(self,elemA,elemB): #O(min(|A|,|B|))
                """Merges the internal subsets of which elemA and elemB are
                elements."""
		if not self.same(elemA, elemB):
			setA = self.map[elemA]
			setB = self.map[elemB]
			shorter = setA if len(setA) < len(setB) else setB
			longer = setB if len(setA) < len(setB) else setA
			for elem in shorter.getSet():
				longer.add(elem)
				self.map[elem] = longer
				self.reps.discard(elem)
#			# Testing
#			for a in setA.getSet():
#				for b in setB.getSet():
#					assert self.same(a,b)
#					assert a in self.getSetRef(b)
#					assert b in self.getSetRef(a)

## Testing
#sets = DisjointSet()
#sets.add(2)
#assert sets.getSetRef(2) == set([2])
#sets.add(3)
#assert sets.getSetRef(3) == set([3])
#assert not sets.same(2,3)
#sets.merge(2,3)
#assert sets.getSetRef(2) == set([2,3])
#assert sets.getSetRef(3) == set([2,3])
#container = set()
#for i in sets.iterate(3):
#	container.add(i)
#assert container == set([2,3])
#assert sets.same(2,3)
#sets.add(4)
#assert sets.getSetRef(4) == set([4])
#sets.add(5,2)
#sets.getSetRef(2) == set([2,3,5])
#sets.getSetRef(3) == set([2,3,5])
#sets.getSetRef(5) == set([2,3,5])
#container = set()
#for i in sets.iterate(3):
#	container.add(i)
#assert container == set([2,3,5])
#assert sets.same(2,3)
#assert sets.same(2,5)
#assert sets.same(3,3)
## Removing like this is safe
#sets.getSetCopy(5).remove(2)
#sets.getSetRef(2) == set([2,3,5])
#sets.getSetRef(3) == set([2,3,5])
#sets.getSetRef(5) == set([2,3,5])
## Removing like this will cause side effects
#sets.getSetRef(5).remove(2)
#sets.getSetRef(2) == set([3,5])
#sets.getSetRef(3) == set([3,5])
#sets.getSetRef(5) == set([3,5])
#sets.getSetRef(5).add(2)
## Removing and set representatives
#assert sets.getSetRep(5) == 2
#sets.remove(2)
#assert sets.getSetRep(5) == 3 or sets.getSetRep(5) == 5
#print('Passed tests.')
