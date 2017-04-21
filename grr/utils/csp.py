import random
from collections import defaultdict
import itertools
from trie import Trie
import time

class Conjunction:
	def __init__(self,pnlist):
		self.pnlist = pnlist
	def __call__(self,*args):
		assert len(args) == len(self.pnlist),"Must have same number of arguments as +/- list"
		for (pn,a) in zip(self.pnlist,args):
			if pn == a:
				return True
		return False

class CSP:
	def __init__(self):
		#variable names
		self.variables = []
		#map from names to indices
		self.varToIndex = {}
		#domains of each variable
		self.varDomains = []
		#constraints: either functions or tables
		self.constraints = []
		#names of each constraint
		self.conNames = []
		#arguments to each constraint
		self.conDomains = []

		#only important messages
		self.verbose = 1

		self.varIncidentConstraints = []
		self.varNeighbors = []
		self.conToIndex = {}
		self.conDomainIncidentConstraints = defaultdict(list)

	def varIndex(self,v):
		if isinstance(v,int):
			return v
		return self.varToIndex[v]
	def conIndex(self,v):
		if isinstance(v,int):
			return v
		return self.conToIndex[v]
	def addVariable(self,domain,name=None):
		if name != None:
			assert name not in self.varToIndex,"Duplicate variable "+str(name)
			assert not isinstance(name,int),"Variable cannot be named an integer"
			self.varToIndex[name] = len(self.variables)
		self.variables.append(name)
		self.varDomains.append(domain)
		self.varIncidentConstraints.append([])
		self.varNeighbors.append(set())
	def addConstraint(self,constraint,vars,name=None):
		"""constraint must either be a function on the variables in vars or a 
		set of valid tuples of the variables in vars."""
		conindex = len(self.constraints)
		if name != None:
			assert name not in self.conToIndex,"Duplicate constraint "+str(name)
			self.conToIndex[name] = conindex
		self.conNames.append(name)
		self.constraints.append(constraint)
		self.conDomains.append([])
		for v in vars:
			v = self.varIndex(v)
			self.conDomains[-1].append(v)
			self.varIncidentConstraints[v].append(conindex)
			for u in vars:
				u = self.varIndex(u)
				if u != v:
					self.varNeighbors[v].add(u)
		self.conDomainIncidentConstraints[tuple(self.conDomains[-1])].append(conindex)
	def assignmentCost(self,var,value,partialAssignment):
		"""Can be overloaded to assign a cost to assigning var=value given the current partial assignment"""
		return 0
	def numConflicts(self,assignment):
		"""Returns the number of conflicts in the given assignment"""
		count = 0
		for c in xrange(len(self.constraints)):
			if not self.testConstraint(c,assignment):
				count += 1
		return count
	def testConstraint_raw(self,con,args):
		if callable(self.constraints[con]):
			return self.constraints[con](*args)
		else:
			return args in self.constraints[con]
	def testConstraint(self,con,assignment):
		"""Tests whether the given assignment satisfies constraint con.  assignment can be a partial
		assignment or full assignment."""
		con = self.conIndex(con)
		args = tuple([assignment[v] for v in self.conDomains[con]])
		return self.testConstraint_raw(con,args)
	def reduceUnaryDomains(self):
		"""For all unary constraints, reduces the domains of the variables and removes the
		constraint from further consideration"""
		for var in xrange(len(self.variables)):
			if (var,) in self.conDomainIncidentConstraints:
				cons = self.conDomainIncidentConstraints[(var,)]
				valid = []
				for value in self.varDomains[var]:
					assignment = {var:value}
					for i in cons:
						if self.testConstraint(i,assignment):
							valid.append(value)
							break
				self.varDomains[var] = valid
				for i in cons:
					self.varIncidentConstraints[var].erase(i)
					self.conDomains[i] = []
					self.constraints[i] = lambda : True
				del self.conDomainIncidentConstraints[(var,)]
	def isCompatible(self,var,value,partialAssignment):
		"""
		- var: the variable to test
		- value: the value to set
		- partialAssignment: a list of assigned values, or None if the indexed variable is unassigned"""
		assert isinstance(var,int)
		assert len(partialAssignment) == len(self.variables)
		prevValue = partialAssignment[var]
		partialAssignment[var] = value
		for c in self.varIncidentConstraints[var]:
			cdomain = self.conDomains[c]
			if all(partialAssignment[v] != None for v in cdomain):
				if not self.testConstraint_raw(c,tuple([partialAssignment[v] for v in cdomain])):
					partialAssignment[var] = prevValue
					return False
		partialAssignment[var] = prevValue
		return True
	def numIncompatible(self,var,value,partialAssignment):
		"""
		- var: the variable to test
		- value: the value to set
		- partialAssignment: a list of assigned values, or None if the indexed variable is unassigned
		"""
		assert isinstance(var,int)
		assert len(partialAssignment) == len(self.variables)
		prevvalue = partialAssignment[var]
		numFailed = 0
		partialAssignment[var] = value
		for c in self.varIncidentConstraints[var]:
			cdomain = self.conDomains[c]
			if all(partialAssignment[v] != None for v in cdomain):
				if not self.testConstraint_raw(c,tuple([partialAssignment[v] for v in cdomain])):
					#print dict(zip([self.variables[v] for v in cdomain],[partialAssignment[v] for v in cdomain])),"Incompatible with constraint",c
					numFailed += 1
		partialAssignment[var] = prevvalue
		return numFailed
	def isCompatible2(self,var1,value1,var2,value2):
		"""Returns true if the variables var1=value1, var2=value2 are pairwise compatible."""
		assert isinstance(var1,int) and isinstance(var2,int)
		if (var1,var2) in self.conDomainIncidentConstraints:
			cons = self.conDomainIncidentConstraints[(var1,var2)]
			for c in cons:
				if self.testConstraint_raw(c,(value1,value2)):
					return True
		if (var2,var1) in self.conDomainIncidentConstraints:
			cons = self.conDomainIncidentConstraints[(var2,var1)]
			for c in cons:
				if self.testConstraint_raw(c,(value2,value1)):
					return True
		return True
	def compatible2(self,var1,value1,var2,domain2):
		"""Returns a list of all values in var2's domain domain2 that are compatible with var1=value1."""
		assert isinstance(var1,int) and isinstance(var2,int)
		assert var2 in self.varNeighbors[var1] and var1 in self.varNeighbors[var2]
		invalid = domain2[:]
		valid = []
		if (var1,var2) in self.conDomainIncidentConstraints:
			cons = self.conDomainIncidentConstraints[(var1,var2)]
			for value2 in invalid:
				compat = False
				for c in cons:
					if self.testConstraint_raw(c,(value1,value2)):
						compat = True
						break
				if compat:
					valid.append(value2)
			for value2 in valid:
				invalid.remove(value2)
		if (var2,var1) in self.conDomainIncidentConstraints:
			cons = self.conDomainIncidentConstraints[(var2,var1)]
			for value2 in invalid:
				compat = False
				for c in cons:
					if self.testConstraint_raw(c,(value2,value1)):
						compat = True
						break
				if compat:
					valid.append(value2)
		return valid
	def compatible(self,var,domain,partialAssignment,varNew=None):
		"""Returns a list of all values in var's domain domain that are compatible with the given partial assignment.
		If varNew is provided, only constraints in which both var and varNew appear will be checked."""
		assert isinstance(var,int)
		if varNew != None:
			assert isinstance(varNew,int)
			assert varNew in self.varNeighbors(var)
		if len(domain)==0:
			return domain
		partialAssignment[var] = domain[0]
		ctest = []
		for c in self.varIncidentConstraints[var]:
			cdomain = self.conDomains[c]
			if varNew == None or varNew in cdomain:
				if all(partialAssignment[v] != None for v in cdomain):
					ctest.append(c)
		newdomain = []
		for value in domain:
			partialAssignment[var] = value
			compat = True
			for c in ctest:
				cdomain = self.conDomains[c]
				if not self.testConstraint_raw(c,tuple([partialAssignment[v] for v in cdomain])):
					compat = False
					break
			if compat:
				newdomain.append(value)
		partialAssignment[var] = None
		return newdomain

	def loadDIMACS(self,fn):
		"""Assumes an empty CSP"""
		nvar = 0
		ncons = 0
		with  open(fn,'r') as f:
			for line in f:
				if line[0] == 'c':
					#comment
					continue
				elif line[0] == 'p':
					args = line.split()
					if len(args) != 4:
						raise IOError("Invalid DIMACS CNF file")
					if args[1] != 'cnf':
						raise IOError("Invalid DIMACS CNF file")
					nvar = int(args[2])
					ncons = int(args[3])
					for v in xrange(nvar):
						self.addVariable(domain=[0,1],name=str(v+1))
				else:
					args = line.split()
					allcons = []
					cons = []
					for a in args:
						if a == '0':
							allcons.append(cons)
							cons = []
						else:
							cons.append(a)
					if cons != []:
						allcons.append(cons)
					for cons in allcons:
						domain = []
						pn = []
						for c in cons:
							if c[0] == '-':
								pn.append(0)
								domain.append(c[1:])
							else:
								pn.append(1)
								domain.append(c)
						self.addConstraint(Conjunction(pn),domain)

	def saveCSPSugar(self,fn):
		"""Saves in .csp file format used by Sugar"""
		f = open(fn,"w")
		f.write("; output from Python CSP class\n")
		mynames = self.variables[:]
		for i,var in enumerate(self.variables):
			if not var[0].isalpha():
				mynames[i] = 'x_'+var
		domainmaps = [None]*len(self.variables)
		for i,var in enumerate(mynames):
			domain = self.varDomains[i]
			vmin = min(domain)
			vmax = max(domain)
			if vmax - vmin >= len(domain):
				domainmaps[i] = dict((v,i) for i,v in enumerate(domain))
				#f.write("(int %s (%s))\n"%(var," ".join(str(v) for v in domain)))
				f.write("(int %s 0 %d)\n"%(var,len(domain)-1))
			else:
				f.write("(int %s %d %d)\n"%(var,vmin,vmax))
		for c,con in enumerate(self.constraints):
			#TODO: specialized output for Conjunction, AllDiff constraints
			cdomain = self.conDomains[c]
			if self.conNames[c] != None:
				f.write(";Constraint %s\n",self.conNames[c])
			for args in itertools.product(*[self.varDomains[v] for v in cdomain]):
				if not self.testConstraint_raw(c,args):
					#add a disjunction
					f.write("(or ")
					for v,val in zip(cdomain,args):
						if domainmaps[v] != None:
							val = domainmaps[v][val]
						f.write("(!= %s %d) "%(mynames[v],val))
					f.write(")\n")
		f.write(";END")
		f.close()


	def sugarMaxAssignment(self,SUGAR_PATH="/home/motion/sugar-v2-2-1/bin/sugar"):
		"""Uses the Sugar solver... must be installed at the given SUGAR_PATH."""
		import subprocess
		self.saveCSPSugar("temp.csp")
		proc = subprocess.Popen([SUGAR_PATH, "-max", "temp.csp"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		(out, err) = proc.communicate()
		if "ERROR" in out or "ERROR" in err:
			if self.verbose > 0:
				print "Error running Sugar CSP solver, program stdout:"
				print out
				print "stderr"
				print err
			raise IOError("Error running Sugar")
		if "OPTIMUM FOUND" not in out:
			if self.verbose > 0:
				print "Error running Sugar CSP solver, program stdout:"
				print out
			raise IOError("Error running Sugar")
		lines = out.split("\n")
		#find the last line of the form "o #"
		oindex = 0
		objective_value = None
		for i,l in enumerate(lines):
			if len(l) > 1 and l[0] == 'o':
				oindex = i
			if l.startswith("s OPTIMUM FOUND"):
				objective_value = l.split()[-1]
		assert(lines[oindex].split()[-1] == objective_value)
		aindex = oindex+1
		variable_settings = {}
		while aindex < len(lines):
			if lines[aindex][0] == 'a':
				if len(lines[aindex].split()) < 3:
					break
				a,var,value=lines[aindex].split()
				variable_settings[var] = int(value)
			else:
				break
			aindex += 1
		mynames = self.variables[:]
		domainmaps = [None]*len(self.variables)
		for i,var in enumerate(self.variables):
			if not var[0].isalpha():
				mynames[i] = 'x_'+var
			domain = self.varDomains[i]
			vmin = min(domain)
			vmax = max(domain)
			if vmax - vmin >= len(domain):
				domainmaps[i] = dict((v,i) for i,v in enumerate(domain))
		namemap = dict((v,i) for (i,v) in enumerate(mynames))

		for var in mynames:
			if var not in variable_settings:
				raise RuntimeError("Unable to find variable "+var+" in Sugar solver output")
		assignment = [None]*len(self.variables)
		for var,value in variable_settings.iteritems():
			if var not in namemap:
				raise RuntimeError("Invalid variable "+var+" in Sugar solver output")
			if domainmaps[namemap[var]] != None:
				value = domainmaps[namemap[var]][value]
			assignment[namemap[var]] = value
		if self.verbose > 0:
			print "Note: temp.csp stores the CSP used by Sugar solver..."
		return assignment


	def gecodeAssignment(self,maximize=False):
		"""Uses the gecode python solver."""
		import gecode
		import gecode.fsa
		from bab_example import bab

		t0 = time.time()
		mynames = self.variables[:]
		domainmaps = [None]*len(self.variables)
		domainranges = [None]*len(self.variables)
		for i,var in enumerate(self.variables):
			if not var[0].isalpha():
				mynames[i] = 'x_'+var
			domain = self.varDomains[i]
			vmin = min(domain)
			vmax = max(domain)
			domainmaps[i] = dict((v,i) for i,v in enumerate(domain))
			domainranges[i] = (0,len(domain)-1)
			#	domainmaps[i] = dict((v,v) for v in domain)
			#	domainranges[i] = (vmin,vmax)
		namemap = dict((v,i) for (i,v) in enumerate(mynames))

		space = gecode.space()
		gecode_vars = []
		gecode_constraint_vars = []
		for (i,var) in enumerate(self.variables):
			gecode_vars.append(space.intvar(domainranges[i][0],domainranges[i][1]))
		for c,con in enumerate(self.constraints):
			#TODO: specialized output for Conjunction, AllDiff constraints
			validtuples = []
			cdomain = self.conDomains[c]
			if not callable(self.constraints[c]):
				validtuples = list(self.constraints[c])
			else:
				for args in itertools.product(*[self.varDomains[v] for v in cdomain]):
					if self.testConstraint_raw(c,args):
						validtuples.append(args)
			assert len(validtuples) > 0
			#map tuples to gecode indices
			for i,args in enumerate(validtuples):
				validtuples[i] = [domainmaps[v][val] for v,val in zip(cdomain,args)]
			if not maximize:
				space.extensional([gecode_vars[v] for v in cdomain],validtuples)
			else:
				#form a trie on the set of valid tuples
				trie = Trie(validtuples)
				#build an FSA on the trie
				vmin = min(0,min(domainranges[v][0] for v in cdomain))
				vmax = max(1,max(domainranges[v][1] for v in cdomain))
				fsa = gecode.fsa.FSA(range(vmin,vmax+1)).add_extremities()
				failnodes = []
				start = fsa.new_node()
				for i,v in enumerate(cdomain):
					failnodes.append(start+i)
				successnode = start+len(cdomain)
				#print "Success-end node connect",successnode,fsa.zend
				fsa.connect(successnode,fsa.zend,0)
				for i in xrange(len(failnodes)-1):
					#print "Fail-fail node connect",failnodes[i],failnodes[i+1]
					for v in xrange(domainranges[cdomain[i]][0],domainranges[cdomain[i]][1]+1):
						fsa.connect(failnodes[i],failnodes[i+1],v)
				#print "Fail-end node connect",failnodes[-1],fsa.zend
				fsa.connect(failnodes[-1],fsa.zend,1)
				assert successnode == len(fsa.states)

				def doconnect(trienode,fsanode,depth=0):
					if trienode.isleaf():
						return
					leaves = []
					fsachildren = []
					for v,tc in trienode.children.iteritems():
						assert v >= vmin and v <= vmax
						if tc.isleaf():
							#print "Success node connect",fsanode,successnode
							fsa.connect(fsanode,successnode,v)
							fsachildren.append(None)
						else:
							#this is slow
							#n = fsa.new_node()
							n = len(fsa.states)+1
							assert n not in fsa.states
							#print "Internal node connect",fsanode,n
							fsa.connect(fsanode,n,v)
							fsachildren.append(n)
					#invalid values go to the corresponding failure node
					for v in xrange(domainranges[cdomain[depth]][0],domainranges[cdomain[depth]][1]+1):
						if v not in trienode.children:
							#print "Fail node connect",fsanode,failnodes[depth]
							fsa.connect(fsanode,failnodes[depth],v)
					for tc,fc in zip(trienode.children.itervalues(),fsachildren):
						if not tc.isleaf():
							doconnect(tc,fc,depth+1)
				doconnect(trie,fsa.zbeg)
				cviolated = space.intvar(0,1)

				"""					
				print "Begins",fsa.gecodify(False)[0]
				print "Ends",fsa.gecodify(False)[1]
				print "Failure track",failnodes
				print "Transitions",fsa.gecodify(False)[2]
				print validtuples
				print "Edges",[v for v in fsa.edges]
				"""
				#debugging
				"""
				for item in validtuples:
					state = fsa.zbeg
					for i,l in enumerate(item):
						next = None
						for s in fsa.succs(state,l):
							assert next==None,"Multiple successors at depth "+str(i)+" in item "+str(item)
							next = s
						assert next != None,"Item "+str(item)+" fails FSA at depth "+str(i)
						state = next
					next = None
					for s in fsa.succs(state,0):
						assert s == fsa.zend
						next = s
					assert next != None,"Item "+str(item)+" fails FSA, last step"
					for s in fsa.succs(state,1):
						raise ValueError("Item "+str(item)+" incorrect in FSA, last step = 0 works")
				"""
				gecode_constraint_vars.append(cviolated)
				space.extensional([gecode_vars[v] for v in cdomain]+[cviolated],fsa)
				if space.status() == 0:
					if self.verbose > 0:
						print "Uh... Gecode stopped on constraint",cdomain
						print "Domains",[domainmaps[v] for v in cdomain]
						print validtuples

					raw_input()
		all_gecode_vars = gecode_vars[:]
		if maximize:
			numconflicts = space.intvar(0,len(self.variables))
			space.linear(gecode_constraint_vars,gecode.IRT_LE,numconflicts)
			space.minimize(numconflicts)
			all_gecode_vars += gecode_constraint_vars
			all_gecode_vars.append(numconflicts)
		if self.verbose > 0:
			print "Time to create Gecode CSP:",time.time()-t0
			print "Solving Gecode CSP with",len(all_gecode_vars),"variables and",sum(len(v) for v in self.varDomains)+3*len(gecode_constraint_vars),"total values"
		space.branch(all_gecode_vars,gecode.INT_VAR_SIZE_MIN,gecode.INT_VAL_MIN)
		best = None
		bestval = len(self.constraints)
		#for sol in space.search():
		for sol in bab(space):
			if maximize:
				if self.verbose > 0:
					#print sol.val(gecode_vars)
					print sol.val(numconflicts),"/",len(gecode_constraint_vars),"constraints unsatisfied"
					#for c in gecode_constraint_vars:
					#	print "  ",sol.val(c)
				if sol.val(numconflicts) < bestval:
					best = sol.val(gecode_vars)
					bestval = sol.val(numconflicts)
			else:
				best = sol.val(gecode_vars)
				break
		if best == None:
			return None
		#map back to original domains
		for i,val in enumerate(best):
			best[i] = self.varDomains[i][val]
		return best

	def backtrackingAssignment(self):
		"""Finds an assignment using backtracking + forward checking / AC3 conflict resolution / heuristics.
		"""
		res = [None]*len(self.variables)
		def ac3(domains,v1,v2):
			#for two values v1 and v2
			d1 = domains[v1]
			d2 = domains[v2]
			numchanged = 0
			changed = True
			while changed:
				changed = False
				newd1 = []
				newd2 = []
				for val1 in d1:
					compatible = False
					for val2 in d2:
						if self.isCompatible2(v1,val1,v2,val2):
							compatible = True
							break
					if not compatible:
						#incompatible
						changed = True
						numchanged += 1
					else:
						newd1.append(val1)
				for val2 in d2:
					compatible = False
					for val1 in newd1:
						if self.isCompatible2(v1,val1,v2,val2):
							compatible = True
							break
					if not compatible:
						#incompatible
						changed = True
						numchanged += 1
					else:
						newd2.append(val2)
				if changed:
					d1 = newd1
					d2 = newd2
			#print "AC3 removed",numchanged,"elements from",self.variables[v1],self.variables[v2]
			return d1,d2
		def run_ac3(domains,unassigned,v1=None,v2=None):
			if v1 == None:
				for v1 in unassigned:
					if not run_ac3(domains,unassigned,v1):
						return False
			elif v2 == None:
				for v2 in self.varNeighbors[v1]:
					if v2 in unassigned:
						if not run_ac3(domains,unassigned,v1,v2):
							return False
			else:
				d1,d2 = ac3(domains,v1,v2)
				domains[v1] = d1
				domains[v2] = d2
				if len(d1) == 0 or len(d2) == 0:
					return False
			return True
		def run_forward_checking(var,domains,assignment):
			for n in self.varNeighbors[var]:
				if assignment[n] != None:
					#already assigned 
					continue
				newdomain = self.compatible(n,domains[n],assignment)
				#if len(newdomain) < len(domains[n]):
				#	print "Forward checking reduced domain of",n,"from",len(domains[n]),"to",len(newdomain)
				if len(newdomain) == 0:
					#contradiction, quit early
					return False
				domains[n] = newdomain
			return True


		#note make sure not to modify the domains variable
		def descend(assignment,unassigned,domains):
			#success!
			if len(unassigned)==0: 
				if self.verbose > 1:
					print "Backtracking success!"
				return True
			#pick a variable according to most constrained + most constraining value heuristics
			bestlist = []
			varorder = []
			#print len(unassigned),"unassigned"
			for var in unassigned:
				domainsize = len(domains[var])
				degree = 0
				for c in self.varIncidentConstraints[var]:
					cdomain = self.conDomains[c]
					if any((m in unassigned and m != var) for m in cdomain):
						degree += 1
				score = (-domainsize,degree)
				varorder.append((score,var))
			varorder = sorted(varorder)[::-1]
			for score,var in varorder[:1]:
				domain = domains[var]
				if len(domain) == 0:
					#backtrack
					if self.verbose > 1:
						print "Backtracking on variable",self.variables[var]
					return False

				#print "Assigning variable",self.variables[var],", branching",-score[0],"and affects",score[1],"neighbors"
				#look at all possible edge assignments and how many neighboring values
				#will remain after assignment
				valueOrder = []
				for value in domain:
					mincompatible = 100000
					numcompatible = 0
					distance = self.assignmentCost(var,value,assignment)
					assignment[var] = value
					for n in self.varNeighbors[var]:
						if n not in unassigned: 
							continue
						ncompatible = 0
						for nval in domains[n]:
							if self.isCompatible(n,nval,assignment):
								ncompatible += 1
								if ncompatible >= mincompatible:
									break
						if ncompatible != 0: numcompatible += 1
						mincompatible = min(ncompatible,mincompatible)
					#print "Number of values compatible with",value,":",mincompatible
					score = (mincompatible,numcompatible,-distance)
					valueOrder.append((score,value))
				valueOrder = sorted(valueOrder)[::-1]
				unassigned.remove(var)
				for (score,value) in valueOrder:
					assignment[var] = value
					do_descent = True
					#TODO: don't copy so many domains
					newdomains = domains[:]
					#now propagate domains so that they are compatible with bestValue (forward checking)
					if not run_forward_checking(var,newdomains,assignment):
						do_descent = False
					if not do_descent:
						#print "Assigning variable",self.variables[var],"=",value,"caused contradiction in forward checking"
						continue

					#run AC3
					for n in self.varNeighbors[var]:
						if n not in unassigned: 
							continue
						if not run_ac3(newdomains,unassigned,n):
							do_descent = False
							break
					if not do_descent:
						#print "Assigning variable",self.variables[var],"=",value,"caused failure in AC3"
						continue

					determined = []
					for n in unassigned:
						if len(newdomains[n]) == 1:
							assignment[n] = newdomains[n][0]
							determined.append(n)
							if not run_forward_checking(n,newdomains,assignment):
								do_descent=False
					if not do_descent:
						for n in determined:
							assignment[n] = None
						continue
					for n in determined:
						unassigned.remove(n)
					#finally do the descent
					if descend(assignment,unassigned,newdomains):
						return True
					#undo determined variables
					for n in determined:
						assignment[n] = None
						unassigned.add(n)
				#failed, restore state for next variable
				unassigned.add(var)
				assignment[var] = None
				#print "Variable",self.variables[var],"didn't work"
				#raw_input()
			if self.verbose > 1:
				print "Tried all values at depth",len(self.variables)-len(unassigned),"/",len(self.variables),", backtracking"
			#raw_input()
			return False

		assignment = [None]*len(self.variables)
		domains = [d[:] for d in self.varDomains]
		unassigned = set(range(len(self.variables)))
		#forward checking on unary domains
		for var in xrange(len(self.variables)):
			if (var,) in self.conDomainIncidentConstraints:
				cons = self.conDomainIncidentConstraints[(var,)]
				valid = []
				for value in self.varDomains[var]:
					assignment = {var:value}
					for i in cons:
						if self.testConstraint(i,assignment):
							valid.append(value)
							break
				domains[var] = valid
		run_ac3(domains,unassigned)
		res = descend(assignment,unassigned,domains)
		if res:
			return assignment
		return None

	def heuristicMaxAssignment(self):
		"""Tries to maximize the number of satisfied constraints using simple heuristics.

		I've done some experimenting with heuristics.
		A typical CSP heuristic is to assign the variable with the fewest remaining values
		(most-constrained-value) and breaks ties with largest degree (most-constraining-variable). 
		The value to which it is assigned is the least-constraining-value heuristic. The only issue
		is that this tries to find a contradiction along a branch as quickly as possible.
		If we want to maximize the number of variables assigned, then perhaps a different ordering would be
		better.  It performs OK.

		Least-constraining variable might be a good heuristic but it's relatively expensive
		least-constraining-value works well for a value assignment.
		"""
		res = [None]*len(self.variables)
		unassigned = set(range(len(self.variables)))
		domains = [d[:] for d in self.varDomains]
		def ac3(v1,v2):
			#for two values v1 and v2
			d1 = [v for v in domains[v1]]
			d2 = [v for v in domains[v2]]
			numchanged = 0
			changed = True
			while changed:
				changed = False
				newd1 = []
				newd2 = []
				for val1 in d1:
					compatible = False
					for val2 in d2:
						if self.isCompatible2(v1,val1,v2,val2):
							compatible = True
							break
					if not compatible:
						#incompatible
						changed = True
						numchanged += 1
					else:
						newd1.append(val1)
				for val2 in d2:
					compatible = False
					for val1 in newd1:
						if self.isCompatible2(v1,val1,v2,val2):
							compatible = True
							break
					if not compatible:
						#incompatible
						changed = True
						numchanged += 1
					else:
						newd2.append(val2)
				if changed:
					d1 = newd1
					d2 = newd2
			if numchanged == 0:
				return domains[v1],domains[v2]
			#print "AC3 removed",numchanged,"elements"
			return d1,d2
		def run_ac3(v1=None,v2=None):
			if v1 == None:
				for v1 in unassigned:
					run_ac3(v1)
			elif v2 == None:
				for v2 in self.varNeighbors[v1]:
					if v2 in unassigned:
						run_ac3(v1,v2)
			else:
				d1,d2 = ac3(v1,v2)
				if len(d1) == 0 or len(d2) == 0:
					#ack, we've eliminated everything
					return
				else:
					domains[v1] = d1
					domains[v2] = d2
		run_ac3()
		def getdegree(n):
			return len([m for m in self.varNeighbors[n] if m in unassigned])
		degrees = dict((n,getdegree(n)) for n in unassigned)
		def getscore(var):
			
			#this is the number of unassigned neighbors of the variable
			degree = degrees[var]
			#most-constrained-value heuristic doesn't work well for minimizing number of conflicts
			#since we're not backtracking
			return (-len(domains[var]),degree)

		scores = dict((v,getscore(v)) for v in unassigned)
		failures = []
		numAssigned = 0
		while len(unassigned) > 0:
			bestlist = []
			bestScore = (-float('inf'),0)
			#find the max scoring node in scores
			for var,score in scores.iteritems():
				assert res[var] == None
				if score > bestScore:
					bestScore = score
					bestlist = [var]
				elif score == bestScore:
					bestlist.append(var)
			#Deterministic choice
			#best = min(bestlist)
			#Random choice among best
			best = random.choice(bestlist)
			domain = domains[best]
			if len(domain) == 0:
				if self.verbose > 1:
					print "Empty domain",self.variables[best],"score",bestScore
				#assign the minimum conflict value from the original domain
				minIncompatible = 10000000
				bestVal = None
				for val in self.varDomains[best]:
					nincompatible = self.numIncompatible(best,val,res)
					if nincompatible < minIncompatible:
						bestVal = val
						minIncompatible = nincompatible
				assert res[best] == None
				assert bestVal != None
				res[best] = bestVal
				numAssigned += 1

				del scores[best]
				del degrees[best]
				unassigned.remove(best)
				failures.append(best)
				continue
			if self.verbose > 1:
				print "Assigning variable",self.variables[best],"score:",bestScore
			#look at all possible edge assignments and how many neighboring values
			#will remain after assignment
			bestValue = None
			bestScore = (-1,0,0)
			for value in domain:
				mincompatible = 100000
				numcompatible = 0
				distance = self.assignmentCost(best,value,res)
				res[best] = value
				for n in self.varNeighbors[best]:
					if n not in unassigned: 
						continue
					ncompatible = 0
					for nval in domains[n]:
						if self.isCompatible(n,nval,res):
							ncompatible += 1
							if ncompatible >= mincompatible:
								break
					if ncompatible != 0: numcompatible += 1
					mincompatible = min(ncompatible,mincompatible)
				#print "Number of values compatible with",value,":",mincompatible
				score = (mincompatible,numcompatible,-distance)
				if score > bestScore:
					bestScore = score
					bestValue = value
			if self.verbose > 1:
				print "  to value",bestValue,"score %d min compatible, %d total compatible, distance %g"%bestScore
			assert bestValue != None
			res[best] = bestValue
			numAssigned += 1
			#if bestScore[0] > 1:
			#	raw_input()
			if bestScore == (0,0):
				if self.verbose > 1:
					print "Warning, not compatible with any neighbors"
					print "Domain:",domain
			#now propagate domains so that they are compatible with bestValue (forward checking)
			for n in self.varNeighbors[best]:
				if n not in unassigned: 
					continue
				newdomain = self.compatible(n,domains[n],res)
				#if len(newdomain) < len(domains[n]):
				#	print "Forward checking reduced domain of",n,"from",len(domains[n]),"to",len(newdomain)
				domains[n] = newdomain
			#run AC3
			for n in self.varNeighbors[best]:
				if n not in unassigned: 
					continue
				run_ac3(n)
			unassigned.remove(best)

			#now update scores and degrees
			del scores[best]
			del degrees[best]
			for n in self.varNeighbors[best]:
				if n in unassigned:
					degrees[n] = getdegree(n)
					scores[n] = getscore(n)
			for e in unassigned:
				assert e in scores
			for e in scores:
				assert e in unassigned
		if len(failures) > 0:
			if self.verbose > 0:
				print "Unique resolution failed, unable to assign",len(failures),"/",len(self.variables),"variables","violating",self.numConflicts(res),"/",len(self.constraints),"constraints"
		return res

	def randomDescentAssignment(self,startingAssignment=None,perturb=False):
		"""Tries to minimize the number of conflicting constraints.

		This heuristic samples a random assignment for each non-assigned value.  Then, for all
		variables that conflict with their neighbors, it attempts to switch the assignment to a less-
		conflicting one.  This repeats until no progress can be made.
		"""
		assignment = [None]*len(self.variables)
		if startingAssignment != None:
			if isinstance(startingAssignment,dict):
				#may be a partial assignment
				for (var,value) in startingAssignment.iteritems():
					assignment[var] = value
			else:
				#assume it's a list
				assert len(startingAssignment) == len(self.variables)
				assignment = startingAssignment[:]
		for i in xrange(len(self.variables)):
			if assignment[i] == None:
				domain = self.varDomains[i]
				if len(domain) == 0: 
					continue
				assignment[i] = random.choice(domain)
			
		def getNumConflicts(i,val=None):
			if val is None:
				val = assignment[i]
			return self.numIncompatible(i,val,assignment)
		def quality(i,val=None):
			if val is None:
				val = assignment[i]
			return self.assignmentCost(i,val,assignment)
		numConflicts = dict((v,getNumConflicts(v)) for v in xrange(len(self.variables)))
		if perturb:
			pRandomizeConflicts = 0.1
			pRandomizeBorder = 0.0
			#randomize conflicting vertices
			for (v,nc) in numConflicts.iteritems():
				if nc > 0 and random.random() < pRandomizeConflicts:
					assignment[v] = random.choice(self.varDomains[v])
				if nc > 0 and random.random() < pRandomizeBorder:
					w = self.varNeighbors[v]
					assignment[w] = random.choice(self.varDomains[w])

		totalConflicted = len([c for c in numConflicts.itervalues() if c != 0])
		active = set([k for k,v in numConflicts.iteritems()])
		while True:
			if len(active) == 0:
				if self.verbose > 1:
					print "Found a conflict-free setting!"
				#greedily build subgraph of conflict-free vertices
				return assignment
			if self.verbose > 1:
				print "Descent pass through active list of size",len(active)
			#print "Active:",[self.variables[v] for v in active]
			changed = False
			changeList = []
			order = list(active)
			random.shuffle(order)
			for i in order:
				iconflicts = numConflicts[i],quality(i)
				bestnewval = None
				bestconflicts = iconflicts
				if iconflicts[0] == 0:
					#only seek to improve quality
					#print "Var",self.variables[i],"=",assignment[i],"quality",iconflicts[1]
					for val in self.varDomains[i]:
						q = quality(i,val)
						if q < bestconflicts[1]:
							bestconflicts = 0,q
							bestnewval = val
				else:
					#print "Var",self.variables[i],"=",assignment[i],"conflicts",iconflicts
					for val in self.varDomains[i]:
						nc = getNumConflicts(i,val),quality(i,val)
						if nc < bestconflicts:
							bestconflicts = nc
							bestnewval = val
				if bestnewval != None:
					changed=True
					changeList.append(i)
					assignment[i] = bestnewval
					#print "  Changed var",self.variables[i],"=",bestnewval,"from",iconflicts,"conflicts to",bestconflicts
					#now update neighbors' conflict counts, wake up
					#dormant neighbors
					for j in self.varNeighbors[i]:
						if assignment[j] != None:
							newconflicts = getNumConflicts(j)
							#if newconflicts < numConflicts[j]:
							#	print "Reduced number of conflicts for",self.variables[j],"too"
							numConflicts[j] = newconflicts
							active.add(j)
					numConflicts[i] = bestconflicts[0]
				else:
					if iconflicts[0] > 0:
						#print "  No change to var",self.variables[i],"still",iconflicts,"conflicts"
						assert iconflicts[0] == numConflicts[i]
				if i in active:
					active.remove(i)
			newConflicted = len([c for c in numConflicts.itervalues() if c != 0])
			if self.verbose > 1:
				print "Beginning # of conflicts",totalConflicted,"now",newConflicted
			totalConflicted = newConflicted
			if self.verbose > 1:
				print "# of changed variables",len(changeList)
			#raw_input()
			if not changed:
				return assignment

if __name__ == '__main__':
	import sys
	csp = CSP()
	fn = 'simple_v3_c2.cnf'
	if len(sys.argv) > 1:
		fn = sys.argv[1]
	csp.loadDIMACS(fn)
	#csp.saveCSPSugar(fn[:-3]+"csp")
	print len(csp.variables),"variables",len(csp.constraints),"constraints"
	#assignment = csp.sugarMaxAssignment()
	assignment = csp.gecodeAssignment(maximize=True)
	#assignment = csp.backtrackingAssignment()
	#assignment = csp.heuristicMaxAssignment()
	#assignment = csp.randomDescentAssignment()
	print dict(zip(csp.variables,assignment)),"has",csp.numConflicts(assignment),"conflicts"
