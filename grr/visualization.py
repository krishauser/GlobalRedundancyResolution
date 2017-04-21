from klampt import *
from klampt.model import trajectory
from klampt.math import vectorops,so3,se3
from klampt import vis
from klampt.vis.glcommon import GLWidgetPlugin
from klampt import PointPoser,TransformPoser
from klampt.vis import gldraw
from klampt.vis.glcommon import CachedGLObject
from utils.metric import *
from ikdb.utils import mkdir_p
from graph import RedundancyResolutionGraph
from OpenGL.GL import *
from utils.utils import *
import os
import json
import time

class GLRedundancyProgram(GLWidgetPlugin):
	def __init__(self,world,robot):
		GLWidgetPlugin.__init__(self)
		self.world = world
		self.robot = robot
		assert robot.index == 0,"Can only work with the 0'th robot in the world"
		self.resolution = RedundancyResolutionGraph(world)
		self.mode = 'interpolate'
		self.configs = None
		self.folder = None
		self.settings = None
		self.drawWorkspaceRoadmap = False
		self.solveConstraint = True
		self.clippingplanes = (0.1,50)
		self.rotationAsDepth = False
		self.pointWidget = PointPoser()
		self.xformWidget = TransformPoser()
		self.roadmapDisplayList = CachedGLObject()
		self.disconnectionDisplayList = CachedGLObject()
		self.movie_frame = None
		self.movie_rotate = False
		self.walk_path = None
		self.walk_workspace_path = None
		self.walk_progress = None
		self.temp_config = None
		self.useboundary = None

	def initialize(self):
		GLBaseClass.initialize(self)
		self.window.program.clearColor = [1,1,1,1]
		
		assert self.resolution.domain != None
		print "Domain",self.resolution.domain
		if len(self.resolution.domain[0]) == 6:
			#if a 2D problem and want to show depth, turn this to true
			self.rotationAsDepth = True
			link = self.resolution.ikTemplate.objectives[0].link()
			self.xformWidget.enableTranslation(True)
			self.xformWidget.enableRotation(True)
			print "Initial transform",self.robot.link(link).getTransform()
			#self.xformWidget.set(*self.robot.link(link).getTransform())
			self.addWidget(self.xformWidget)
		else:
			self.rotationAsDepth = False
			link = self.resolution.ikTemplate.objectives[0].link()
			local,world = self.resolution.ikTemplate.objectives[0].getPosition()
			print "Initial position",self.resolution.robot.link(link).getWorldPosition(local)
			self.pointWidget.set(self.resolution.robot.link(link).getWorldPosition(local))
			self.addWidget(self.pointWidget)
		return True

	def run(self):
		"""Standardized interface for running program"""
		vp = vis.getViewport()
		#Square screen
		#vp.w,vp.h = 800,800
		#For saving HD quality movies
		vp.w,vp.h = 1024,768
		vp.clippingplanes = self.clippingplanes
		vis.setViewport(vp)
		#vis.run(program)
		vis.setPlugin(self)
		vis.show()
		while vis.shown():
			time.sleep(0.1)
		vis.setPlugin(None)
		vis.kill()

	def workspaceToPoint(self,x):
		if self.rotationAsDepth:
			return (x[0],x[4],x[2])
		else:
			return x[:3]

	def display(self):
		if self.configs == None:
			self.robot.drawGL()
		else:
			for q in self.configs:
				self.robot.setConfig(q)
				self.robot.drawGL()

		if self.walk_workspace_path != None:
			if self.temp_config:
				self.robot.setConfig(self.temp_config)
				glEnable(GL_BLEND)
				#glDisable(GL_DEPTH_TEST)
				for i in xrange(self.robot.numLinks()):
					self.robot.link(i).appearance().setColor(1,0,0,0.5)
				self.robot.drawGL()
				for i in xrange(self.robot.numLinks()):
					self.robot.link(i).appearance().setColor(0.5,0.5,0.5,1)
				#glEnable(GL_DEPTH_TEST)
				glDisable(GL_BLEND)
			glColor3f(1,1,0)
			glLineWidth(5.0)
			glBegin(GL_LINE_STRIP)
			for w in self.walk_workspace_path.milestones:
				if len(w) == 2:
					glVertex2f(w[0],w[1])
				else:
					glVertex3f(w[0],w[1],w[2])
			glEnd()
			glLineWidth(1.0)


		#draw workspace graph
		if self.drawWorkspaceRoadmap:
			def drawWorkspaceRoadmap():
				if not self.resolution.hasResolution():
					glDisable(GL_LIGHTING)
					glPointSize(5.0)
					glColor3f(1,0,0)
					glBegin(GL_POINTS)
					for n,d in self.resolution.Gw.nodes_iter(data=True):
						glVertex3fv(self.workspaceToPoint(d['params']))
					glEnd()
					glColor3f(1,0.5,0)
					glBegin(GL_LINES)
					for (i,j) in self.resolution.Gw.edges_iter():
						glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[i]['params']))
						glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[j]['params']))
					glEnd()
				else:
					#maxn = max(len(d['qlist']) for n,d in self.resolution.Gw.nodes_iter(data=True))
					#maxe = max(len(d['qlist']) for i,j,d in self.resolution.Gw.edges_iter(data=True))
					#maxn = max(maxn,1)
					#maxe = max(maxe,1)
					glDisable(GL_LIGHTING)
					glPointSize(5.0)
					glBegin(GL_POINTS)
					for n,d in self.resolution.Gw.nodes_iter(data=True):
						if not self.resolution.isResolvedNode(n):
							continue
						#u = float(len(d['qlist']))/float(maxn)
						#nsubset = sum(1 for iq in d['qlist'] if iq in self.rr.Qsubgraph)
						#if nsubset > 1:
						#	glColor3f(1,0,1)
						#else:
						#	glColor3f(u,nsubset*0.5,0)
						glColor3f(1,0.5,0)
						glVertex3fv(self.workspaceToPoint(d['params']))
					glEnd()
					glBegin(GL_LINES)
					for (i,j,d) in self.resolution.Gw.edges_iter(data=True):
						if not self.resolution.isResolvedNode(i) or not self.resolution.isResolvedNode(j):
							continue
						#nsubset = sum(1 for (ia,ib) in d['qlist'] if (ia in self.Qsubgraph and ib in self.Qsubgraph))
						#u = float(len(d['qlist']))/float(maxe)
						r,g,b = 1,1,0
						if not self.resolution.isResolvedEdge(i,j):
							r,g,b = 1,0,1
						else:
							qd = self.robot.distance(self.resolution.Gw.node[i]['config'],self.resolution.Gw.node[j]['config'])
							wd = workspace_distance(self.resolution.Gw.node[i]['params'],self.resolution.Gw.node[j]['params'])
							u = 1.0 - math.exp(-max(0.0,2.0*qd/wd-1.0))
							if u > 0.8:
								r,g,b = 1,1.0-u*5.0,0.0
							else:
								r,g,b = u/0.8,1,0.0
						#assert (nsubset >=1) == self.resolution.isResolvedEdge(i,j),"Mismatch between Qsubgraph and resolution?"
						glColor3f(r,g,b)
						glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[i]['params']))
						glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[j]['params']))
					glEnd()
					"""
					glEnable(GL_LIGHTING)
					q0 = self.robot.getConfig()
					for iw,d in self.rr.Gw.nodes_iter(data=True):
						qs = [iq for iq in d['qlist'] if iq in self.rr.Qsubgraph]
						if len(qs) > 1:
							for iq in qs:
								self.robot.setConfig(self.rr.Gq.node[iq]['config'])
								self.robot.drawGL()
					self.robot.setConfig(q0)
					glDisable(GL_LIGHTING)
					"""
			#self.roadmapDisplayList.draw(drawWorkspaceRoadmap)
			drawWorkspaceRoadmap()
		else:
			#render boundaries only
			def drawDisconnections(orientation=None):
				bmin,bmax = self.resolution.domain
				active = [i for i,(a,b) in enumerate(zip(bmin,bmax)[:3]) if b!=a]
				if self.useboundary == None:
					self.useboundary = (len(active)==2)
				verts,faces = self.resolution.computeDiscontinuities(self.useboundary,orientation)
				
				if len(active)==3:
					glEnable(GL_LIGHTING)
					glEnable(GL_BLEND)
					glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
					#if self.walk_workspace_path != None:
					#	gldraw.setcolor(1,0,1,0.25)
					#else:
					#	gldraw.setcolor(1,0,1,0.5)

					for tri in faces:
						tverts = [verts[v] for v in tri]
						centroid = vectorops.div(vectorops.add(*tverts),len(tri))
						params = [(x-a)/(b-a) for (x,a,b) in zip(centroid,self.resolution.domain[0],self.resolution.domain[1])]
						if self.walk_workspace_path != None:
							gldraw.setcolor(params[0],params[1],params[2],0.25)
						else:
							gldraw.setcolor(params[0],params[1],params[2],1.0)
						gldraw.triangle(*tverts)
				else:
					glDisable(GL_LIGHTING)
					glEnable(GL_BLEND)
					glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
					glColor4f(1,0,1,0.5)
				
					glLineWidth(4.0)
					glBegin(GL_LINES)
					for s in faces:
						spts = verts[s[0]],verts[s[1]]
						glVertex3f(spts[0][0],0,spts[0][1])
						glVertex3f(spts[1][0],0,spts[1][1])
					glEnd()
					glLineWidth(1.0)
			#draw stabbing lines
			"""
			glDisable(GL_LIGHTING)
			glColor3f(1,0,0)
			glBegin(GL_LINES)
			for (i,j,d) in self.resolution.Gw.edges_iter(data=True):
				if self.resolution.isResolvedEdge(i,j):
					continue
				if not self.resolution.isResolvedNode(i) or not self.resolution.isResolvedNode(j):
					continue
				glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[i]['params']))
				glVertex3fv(self.workspaceToPoint(self.resolution.Gw.node[j]['params']))
			glEnd()
			"""
			orientation = None
			if len(self.resolution.domain[0]) == 6:
				orientation = self.xformWidget.get()[0]
			self.disconnectionDisplayList.draw(drawDisconnections,args=(orientation,))

			bmin,bmax = self.resolution.domain

			#draw workspace bound
			#glDisable(GL_LIGHTING)
			#glColor3f(1,0.5,0)
			#gldraw.box(bmin[:3],bmax[:3],lighting=False,filled=False)


		GLBaseClass.display(self)

	def mousefunc(self,button,state,x,y):
		GLBaseClass.mousefunc(self,button,state,x,y)
		self.do_click(x,y)

	def do_click(self,x,y):
		if self.pointWidget.hasFocus():
			p = self.pointWidget.get()
			if len(self.resolution.domain[0]) == 6:
				x = p+[0,0,0]
			else:
				x = p
			if self.mode == 'interpolate':
				if self.solveConstraint:
					qOrig = self.configs[0] if self.configs is not None else self.robot.getConfig()
					qx = self.resolution.eval(x,qOrig) 
					self.configs = [qx if qx is not None else qOrig ]
				else:
					self.configs = self.resolution.interpolate(x,multi=True)
			else:
				self.configs = self.rr.select(x,0.05,0.02)
			self.refresh()
		if self.xformWidget.hasFocus():
			R,t = self.xformWidget.get()
			m = so3.moment(R)
			x = t + m
			if self.mode == 'interpolate':
				if self.solveConstraint:
					qOrig = self.configs[0] if self.configs is not None else self.robot.getConfig()
					qx = self.resolution.eval(x,qOrig) 
					self.configs = [qx if qx is not None else qOrig ]
				else:
					self.configs = self.resolution.interpolate(x,multi=True)
			else:
				self.configs = self.rr.select(x,0.05,0.02)
			self.refresh()

	def motionfunc(self,x,y,dx,dy):
		GLBaseClass.motionfunc(self,x,y,dx,dy)
		self.do_click(x,y)
		
	def keyboardfunc(self,c,x,y):
		if c == 'h':
			print "Keyboard help:"
			print "- x: verifies edge checks for existing resolution"
			print "- i: toggles between interpolation mode and graph inspection mode"
			print "- g: toggles drawing the workspace graph"
			print "- w: performs a walk to a random workspace node"
			print "- m: saves a real-time movie"
			print "- M: saves a 360 spin movie"
		elif c == 'x':
			self.resolution.verify()
			self.disconnectionDisplayList.markChanged()
		elif c == 'i':
			if self.mode == 'interpolate':
				self.mode = 'inspect'
			else:
				self.mode = 'interpolate'
			print "Toggled visualization mode to",self.mode
		elif c == 'g':
			self.drawWorkspaceRoadmap = not self.drawWorkspaceRoadmap
		elif c == 'b':
			if self.useboundary: self.useboundary = False
			else: self.useboundary = True
			self.disconnectionDisplayList.markChanged()
		elif c == 'm':
			if self.movie_frame is None:
				self.movie_frame = 0
				self.movie_rotate = False
			else:
				self.movie_frame = None
		elif c == 'M':
			if self.movie_frame is None:
				self.movie_frame = 0
				self.movie_rotate = True
			else:
				self.movie_frame = None
		elif c == 'w':
			import random
			resolved = []
			for iw,d in self.resolution.Gw.nodes(data=True):
				if d.get('config',None) is not None:
					resolved.append(iw)
			if self.walk_workspace_path == None:
				#draw boundaries transparent now
				self.disconnectionDisplayList.markChanged()
			for iters in range(10):
				wtgt = random.choice(resolved)
				#TEMP: place widgets far away for capturing video
				far = [20]*3
				self.pointWidget.set(far)
				R,t = self.xformWidget.get()
				self.xformWidget.set(R,far[:3])

				#get current config
				if self.configs != None:
					self.robot.setConfig(self.configs[0])
				try:
					x,p = self.resolution.walkPath(wtgt)
				except Exception as e:
					import traceback
					print "Exception",e
					traceback.print_exc()
					print "WalkPath failed, trying with a new random point"
					continue
				self.walk_workspace_path = None
				if p != None:
					t = 0
					times = []
					for i,q in enumerate(p):
						times.append(t)
						if i+1 < len(p):
							t += linf_distance(q,p[i+1],self.resolution.spinJoints)
						#t += 0.1
					self.walk_workspace_path = trajectory.Trajectory(times,x)
					self.walk_path = trajectory.RobotTrajectory(self.robot,times,p)
					self.walk_progress = 0
					if self.temp_config == None:
						self.temp_config = p[0]
				break
		self.refresh()

	def idle(self):
		t0 = time.time()
		if self.movie_frame is not None:
			numframes = 180
			if self.movie_rotate:
				self.window.program.view.camera.rot[2] += math.pi*2 / numframes
			self.window.program.save_screen("frame%03d.ppm"%(self.movie_frame))
			self.movie_frame += 1
			if self.movie_rotate and self.movie_frame >= numframes:
				self.movie_frame = None
		if self.walk_path != None:
			self.walk_progress += 0.02
			if self.walk_progress >= self.walk_path.times[-1]:
				self.configs = [self.walk_path.milestones[-1]]
				self.walk_path = None
			else:
				#self.configs = [self.walk_path.eval(self.walk_progress)]
				self.configs = [self.resolution.eval(self.walk_workspace_path.eval(self.walk_progress))]
				if self.configs[0] == None:
					self.configs = []
				u = self.walk_progress / self.walk_path.times[-1]
				qstraight = self.resolution.solve(self.temp_config,vectorops.interpolate(self.walk_workspace_path.milestones[0],self.walk_workspace_path.milestones[-1],u))
				if qstraight and (self.resolution.ikTemplate.feasibilityTest==None or self.resolution.ikTemplate.feasibilityTest(qstraight)):
					self.temp_config = qstraight
			self.refresh()
			t1 = time.time()
			if t1 - t0 < 0.02:
				time.sleep(0.02 - (t1-t0))

def read_options(problem,options,settings_file=None):
	f = None
	if settings_file == None:
		import glob
		#try opening the first json file there?
		files = glob.glob(os.path.join('problems',problem,"*.json"))
		if len(files) == 0:
			raise ValueError("No valid json files in folder problems/"+problem)
		settings_file = files[0]
	else:
		if problem != None:
			settings_file = os.path.join('problems',problem,settings_file)
	f = open(settings_file,'r')
	pdef = toUtf8(json.load(f))
	f.close()
	#override with options
	pdef['problem'] = problem
	pdef.update(options)

	klampt_directory = '/home/motion/Klampt/'
	#klampt_directory = '/Program Files (x86)/Klampt/'
	filename = pdef['filename']
	if not os.path.exists(filename):
		if os.path.exists(os.path.join(klampt_directory,filename)):
			filename = os.path.join(klampt_directory,filename)
			pdef['filename'] = filename
		elif os.path.exists(os.path.join(os.path.split(settings_file)[0],filename)):
			filename = os.path.join(os.path.split(settings_file)[0],filename)
			pdef['filename'] = filename

	Nw=pdef.get('Nw',1000)
	method=pdef.get('workspace_graph','staggered_grid')
	if 'output' in pdef:
		folder = pdef['output']
	else:
		if 'output_prefix' in pdef:
			folder = "output/%s/"%(problem,)+pdef['output_prefix']
		else:
			folder = "output/%s/"%(problem,)
			if settings_file != 'settings.json':
				folder = folder + os.path.splitext(os.path.basename(settings_file))[0]
		folder = folder + "_W%d"%(Nw,)
		if method == 'grid':
			folder = folder + '_Ggrid'
		elif method == 'random':
			folder = folder + '_Grandom'
		pdef['output'] = folder
	return pdef

def parse_args(argv):
	problem = argv[1]

	#command line options
	opts = ["Nw","Nq","workspace_graph","output","output_prefix","link","links","domain","orientation"]
	stringopts = ["workspace_graph","output","output_prefix","link","orientation"]

	settings_file = None
	options = {}
	#parse command line options
	for i in argv[2:]:
		parts = i.split('=')
		if len(parts) == 1:
			if i.endswith('json'):
				settings_file = i
			elif i == 'clean':
				options['clean'] = True
			else:
				raise RuntimeError("Invalid command line argument, must be of form OPTION=VAL")
		else:
			if parts[0] not in opts:
				print "Command line option must be one of:",opts
				raise RuntimeError("Invalid command line option")
			if parts[0] not in stringopts:
				options[parts[0]] = eval(parts[1])
			else:
				options[parts[0]] = parts[1]

	if problem.endswith('json'):
		#for loading settings.json files saved in the output folders
		settings_file = problem
		problem = None
	return read_options(problem,options,settings_file)

def make_program(pdef):
	"""Loads a file from disk and returns an appropriately configured
	GLRedundancyProgram, given a JSON structure defining the problem.
	"""
	Nw=pdef.get('Nw',1000)
	method=pdef.get('workspace_graph','staggered_grid')
	filename = pdef['filename']
	folder = pdef['output']

	world = WorldModel()
	world.readFile(filename)
	robot = world.robot(0)
	program = GLRedundancyProgram(world,robot)
	program.resolution.readJsonSetup(pdef)
	program.folder = folder

	if not pdef.get('clean',False):
		try:
			print "Loading from",folder,"..."
			program.resolution.load(folder)
			#TEMP: DO PADDING FOR BAXTER
			for iq,d in program.resolution.Gw.nodes_iter(data=True):
				if 'config' in d:
					if len(d['config']) < robot.numLinks():
						d['config'] = d['config'] + [0.0]*(robot.numLinks()-len(d['config']))
		except IOError as e:
			print "Did not successfully load from folder %s, generating from scratch"%(folder,)
			program.resolution.sampleWorkspace(Nw,method=method)
	else:
		program.resolution.sampleWorkspace(Nw,method=method)

	#TEMP: make workspace grid following exact
	#if problem.startswith('planar'):
	#	for n,data in program.rr.Gw.nodes_iter(data=True):
	#		data['params'][1] = 0

	program.settings = pdef
	return program


if __name__ == "__main__":
	import sys
	if len(sys.argv) ==0 :
		print "USAGE: redundancyvisualization.py problem [settings json file] [option1=val1 option2=val2 ...]"

	opts = parse_args(sys.argv)
	program = make_program(opts)
	program.run()
	exit(0)

