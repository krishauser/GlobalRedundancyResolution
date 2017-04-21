from klampt import *
from klampt.math import vectorops,so3,se3
from grr.graph import RedundancyResolutionGraph
from grr.visualization import GLRedundancyProgram,parse_args
from grr.solver import RedundancyResolver,ParallelResolver
from grr.utils.utils import *
from OpenGL.GL import *
import os
import time
import json



class RedundancyProgram(GLRedundancyProgram):
	def __init__(self,world,robot,cache=2000):
		GLRedundancyProgram.__init__(self,world,robot)
		self.rr = RedundancyResolver(self.resolution,cache=cache)
		self.pr = ParallelResolver(self.rr)
		self.drawCSpaceRoadmap = False

	def display(self):
		GLRedundancyProgram.display(self)

		if self.drawCSpaceRoadmap:
			#draw workspace - c-space roadmap
			glDisable(GL_LIGHTING)
			glEnable(GL_BLEND)
			glDisable(GL_DEPTH_TEST)
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
			glPointSize(5.0)
			xshift = -4
			glBegin(GL_POINTS)
			for n,d in self.rr.Gq.nodes_iter(data=True):
				x = self.rr.Gw.node[d['windex']]['params']
				q = d['config']
				depth = (-q[2] + 1.0)*0.25
				if n in self.rr.Qsubgraph:
					continue
				else:
					glColor4f(1,1,depth,0.2)
				glVertex3f(x[2]+xshift,q[2],q[1])
			glEnd()
			glBegin(GL_LINES)
			for i,j,d in self.rr.Gq.edges_iter(data=True):
				if i > j: continue
				xi = self.rr.Gw.node[self.rr.Gq.node[i]['windex']]['params']
				xj = self.rr.Gw.node[self.rr.Gq.node[j]['windex']]['params']
				qi = self.rr.Gq.node[i]['config']
				qj = self.rr.Gq.node[j]['config']
				depth = (-qi[2] + 1.0)*0.25
				a = 1.0/(1.0+5.0*self.rr.robot.distance(qi,qj))
				if i in self.rr.Qsubgraph and j in self.rr.Qsubgraph:
					continue
				else:
					glColor4f(1,1,depth,a)
				glVertex3f(xi[2]+xshift,qi[2],qi[1])
				glVertex3f(xj[2]+xshift,qj[2],qj[1])
			glEnd()
			#draw path
			glColor3f(1,0.5,0)
			glBegin(GL_POINTS)
			for n,d in self.rr.Gq.nodes_iter(data=True):
				x = self.rr.Gw.node[d['windex']]['params']
				q = d['config']
				if n in self.rr.Qsubgraph:
					glVertex3f(x[2]+xshift,q[2],q[1])
			glEnd()
			glBegin(GL_LINES)
			for i,j,d in self.rr.Gq.edges_iter(data=True):
				if i > j: continue
				xi = self.rr.Gw.node[self.rr.Gq.node[i]['windex']]['params']
				xj = self.rr.Gw.node[self.rr.Gq.node[j]['windex']]['params']
				qi = self.rr.Gq.node[i]['config']
				qj = self.rr.Gq.node[j]['config']
				if i in self.rr.Qsubgraph and j in self.rr.Qsubgraph:
					glVertex3f(xi[2]+xshift,qi[2],qi[1])
					glVertex3f(xj[2]+xshift,qj[2],qj[1])
			glEnd()
			glEnable(GL_DEPTH_TEST)
			glDisable(GL_BLEND)


	def keyboardfunc(self,c,x,y):
		if c == 'h':
			print "Keyboard help:"
			#print "- [space]: samples more points in the workspace (not very functional)"
			print "INTERACTIVE SOLVER"
			print "- q: samples more points in the configuration space (10 per workspace point)"
			print "- Q: samples more points in the configuration space attempting to add more in poorly-sampled areas"
			print "- c: samples more configurations at conflicting workspace edges"
			print "- u: performs a CSP assignment"
			print "- d: min-conflicts descent of the assignment"
			print "- r: randomizes the assignment and then descends"
			print "- p: performs pointwise redundancy resolution"
			print "- o: 10 iterations of coordinate descent optimization of overall path length"
			print "VISIBILITY GRAPH SOLVER"
			print "- v: build the network using the visibility graph algorithm"
			print "- R: solves a network built using the visibility graph algorithm"
			print "PARALLEL SOLVER"
			print "- t: multithreaded parallel solve"
			print "INSPECTION / VISUALIZATION"
			print "- s: saves roadmap and resolution to disk"
			print "- x: validate the existing resolution"
			print "- i: toggles between interpolation mode and graph inspection mode"
			print "- g: toggles drawing the workspace graph"
			print "- G: toggles drawing the C-space roadmap (only valid for planar problems)"
			print "- w: performs a walk to a random workspace node"
			print "- m: saves a real-time movie"
			print "- M: saves a 360 spin movie"
		elif c == ' ':
			pass
			#print "Sampling workspace..."
			#self.resolution.sampleWorkspace([-3,0,-3],[3,0,3],10)
		elif c == 'q':
			print "Sampling configuration space..."
			self.rr.sampleConfigurationSpace(10)
		elif c == 'Q':
			print "Sampling configuration space..."
			self.rr.sampleConfigurationSpace(10,biased=True)
		elif c == 'c':
			print "Sampling configurations at conflicting edges..."
			self.rr.sampleConflicts()
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 'u':
			print "Performing unique assignment..."
			self.rr.uniqueAssignment()
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 'r':
			print "Doing random descent..."
			self.rr.randomDescentAssignment(True)
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 'd':
			print "Doing descent..."
			self.rr.randomDescentAssignment(False)
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 's':
			print "Saving results..."
			if self.folder == None:
				raise RuntimeError("folder element was not set?")
			self.rr.save(self.folder)
			if self.settings != None:
				f = open(os.path.join(self.folder,'settings.json'),'w')
				json.dump(self.settings,f)
				f.close()
		elif c == 'p':
			print "Pointwise assignment..."
			self.rr.pointwiseAssignment()
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 'o':
			print "Optimizing..."
			self.rr.optimize()
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 'v':
			self.rr.sampleConfigurationSpace(self.Nq,connect=True,visibility=True,k=20)
		elif c == 'R':
			self.rr.solveCSP(visibility=True)
			self.rr.printStats()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		elif c == 't':
			print("There are " + str(multiprocessing.cpu_count()) + " CPUs available.")
			self.pr.setThreads(int(raw_input("Enter the number of threads you would like to spawn: ")))
			self.pr.spawnWorkers()
			self.pr.superviseQSampling(Nq=self.Nq,k=self.visK)
			self.pr.superviseEdgeTests()
			self.pr.writeReport()
			self.pr.solveCSP()
			self.pr.killWorkers()
			self.roadmapDisplayList.markChanged()
			self.disconnectionDisplayList.markChanged()
		else:
			GLRedundancyProgram.keyboardfunc(self,c,x,y)
				

def make_program(pdef):
	"""Loads a file from disk and returns an appropriately configured
	RedundancyProgram.
	"""
	print('made it here')
	makeParams = pdef
	Nw=pdef.get('Nw',1000)
	wgraphtype=pdef.get('workspace_graph','staggered_grid')	
	pdef['workspace_graph'] = wgraphtype
	cache=pdef.get('cache',20000)
	filename = pdef['filename']
	folder = pdef['output']

	world = WorldModel()
	world.readFile(filename)
	robot = world.robot(0)
	program = RedundancyProgram(world,robot,cache=cache)

	program.workspace_graph = wgraphtype
	program.Nq = pdef.get('Nq',100)
	program.visK = pdef.get('vis')
	program.pr.makeParams = makeParams
	program.resolution.readJsonSetup(pdef)
	program.folder = folder

	if not pdef.get('clean',False):
		try:
			print "Loading from",folder,"..."
			program.rr.load(folder)
			#TEMP: DO PADDING FOR BAXTER
			for iw,d in program.rr.Gw.nodes_iter(data=True):
				for i,q in enumerate(d['qlist']):
					while len(q) < robot.numLinks():
						q.append(0.0)
					d['qlist'][i] = q
			for iq,d in program.resolution.Gw.nodes_iter(data=True):
				if 'config' in d and len(d['config']) < robot.numLinks():
					d['config'] = d['config'] + [0.0]*(robot.numLinks()-len(d['config']))
			program.rr.printStats()
		except IOError as e:
			print "Did not successfully load from folder %s, generating from scratch"%(folder,)
			program.resolution.sampleWorkspace(Nw,method=wgraphtype)
			program.rr.initWorkspaceGraph()
			print "Outputs will be saved to folder",folder
	else:
		program.resolution.sampleWorkspace(Nw,method=wgraphtype)
		program.rr.initWorkspaceGraph()
		print "Outputs will be saved to folder",folder

	#TEMP: make workspace grid following exact
	#if problem.startswith('planar'):
	#	for n,data in program.rr.Gw.nodes_iter(data=True):
	#		data['params'][1] = 0

	program.settings = pdef
	return program


if __name__ == "__main__":
	import sys
	if len(sys.argv) == 1:
		print "Usage: redundancy.py problem [json file] [opt1=arg1] [opt2=arg2] ..."

	opts = parse_args(sys.argv)
	program = make_program(opts)
	program.run()
	exit(0)
