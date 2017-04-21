from collections import defaultdict

class Trie:
	def __init__(self,items=None,values=None):
		self.value = None
		self.children = {}
		if items != None:
			self.build(items,values)
	def haslabel(self):
		return self.value != None
	def isleaf(self):
		return len(self.children)==0
	def lookup(self,value):
		if len(value)==0:
			return self.value
		if self.isleaf():
			return None
		try:
			return self.children[value[0]].lookup(value[1:])
		except KeyError:
			return None
	def build(self,items,values=None):
		if values == None:
			values = [True]*len(items)
		prefixes = defaultdict(list)
		childvalues = defaultdict(list)
		for it,val in zip(items,values):
			if len(it) == 0:
				if self.value != None:
					print "Trie: warning, more than one value specified for an item"
				self.value = val
			else:
				prefixes[it[0]].append(it[1:])
				childvalues[it[0]].append(val)
		self.children = {}
		for prefix,childitems in prefixes.iteritems():
			self.children[prefix] = Trie(childitems,childvalues[prefix])
