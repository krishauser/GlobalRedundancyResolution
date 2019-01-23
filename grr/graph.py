from klampt import *
from klampt.model import ik,trajectory,collide
from klampt.math import vectorops,so3,se3
from ikdb.ikproblem import IKProblem,IKSolverParams
from utils.so3_grid import *
from utils.sphere_grid import *
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

DEFAULT_EDGE_CHECK_TOLERANCE = 5e-2

def from_spherical(longitude,latitude,r=1):
    dx = math.cos(longitude)
    dy = math.sin(longitude)
    z = math.sin(latitude)
    dxy = math.cos(latitude)
    return [dxy*dx*r,dxy*dy*r,z*r]

def spherical(x,radius=True):
    if radius:
        r = [vectorops.norm(x)]
    else:
        r = []
    den = math.sqrt(x[1]**2 + x[0]**2)
    if den < 1e-8:
        return [0,math.pi/2*x[2]/abs(x[2])] + r
    else:
        return [math.atan2(x[1],x[0]),math.atan2(x[2],den)] + r


class IKTask:
    """Members:"
    - objective: an IKObjective storing the current setting of the workspace parameters
    - numConstraints: the number of elements in the constraint parameter vector
    """
    def __init__(self,link,eelocal,orientation,orientationparams=None):
        """
        - link: the link to be constrained
        - eelocal: the local coordinates of the constrained point
        - orientation: a string describing the rotational IK constraint. Can be:
          - free: no constraint (3DOF space)
          - variable: 6 DOF space
          - fixed: fixed orientation (3DOF space). orientationparams must be an so3 element describing the local -> world rotation.
          - axis: local rotation about an axis is not constrained  (2DOF space). orientationparams must be a 3-vector describing the local rotation axis.
        - orientationparams: only meaningful with fixed and axis constraints
        """
        self.objective = None 
        self.orientation = orientation
        self.numConstraints = 0
        if orientation == 'free':
            self.objective = ik.objective(link,local=eelocal,world=[0,0,0])
            self.numConstraints = 3
        elif orientation == 'variable':
            T = se3.identity()
            obj = ik.objective(link,R=T[0],t=T[1])
            obj.setFixedPosConstraint(eelocal,[0,0,0])
            self.objective = obj
            self.numConstraints = 6
        elif orientation == 'fixed':
            assert isinstance(orientationparams,(list,tuple)) and len(orientationparams)==9,"must provide an so3 element as parameters for the orientation matrix"
            #template is a fixed orientation matrix 
            obj = ik.objective(link,R=orientationparams,t=[0,0,0])
            obj.setFixedPosConstraint(eelocal,[0,0,0])
            self.objective = obj
            self.numConstraints = 3
        elif orientation == 'axis':
            #template is a fixed axis matrix
            obj = ik.objective(link,local=eelocal,world=[0,0,0])
            obj.setAxialRotConstraint(orientationparams,[0,0,1])
            self.objective = obj
            self.numConstraints = 5
        else:
            raise ValueError("orientation can only be free, variable, fixed, or axis")

    def setParams(self,x):
        """Sets the current IK constraint to the workspace parameters x"""
        assert len(x) == self.numConstraints
        np = self.objective.numPosDims()
        nr = self.objective.numRotDims()
        assert np == 3 or np == 0,"Can only handle fixed or free position constraints"
        assert nr == 3 or nr == 0 or nr == 2,"Can only handle fixed, free, or axis rotation constraints"
        if np == 3:
            assert len(x) >= 3
            local,world = self.objective.getPosition()
            self.objective.setFixedPosConstraint(local,x[0:3])
        if nr == 3 and self.numConstraints == 6:
            R = so3.from_moment(x[np:np+3])
            self.objective.setFixedRotConstraint(R)
        if nr == 2:
            al,aw = self.objective.getRotationAxis()
            #latitude and longitude
            lon = x[np+0]
            lat = x[np+1]
            self.objective.setAxialRotConstraint(al,from_spherical(lon,lat))

    def getParams(self,T):
        """Gets the workspace parameters for this constraint that match T exactly"""
        np = self.objective.numPosDims()
        nr = self.objective.numRotDims()
        assert np == 3 or np == 0,"Can only handle fixed or free position constraints"
        assert nr == 3 or nr == 0 or nr == 2,"Can only handle fixed, free, or axis rotation constraints"
        if np == 3:
            #local,world = self.objective.getPosition()
            #x = se3.apply(T,local)
            x = T[1]
        else:
            x = []
        if nr == 3 and self.numConstraints == 6:
            x += so3.moment(T[0])
        if nr == 2:
            al,aw = self.objective.getRotationAxis()
            aw = so3.apply(T[0],al)
            x += spherical(aw,radius=False)
        return x

    def distance(self,a,b):
        """Returns a distance in workspace between points a and b.  If len(a)<=3, 
        this will do cartesian distance.  Otherwise, len(a)==self.numConstraints is required and
        this will do R^3 x S(2) (axis) or SE(3) (variable) distance (with a max orientation distance of 2)."""
        if len(a) <= 3:
            return vectorops.distance(a,b)
        assert len(a)==self.numConstraints,"Can only interpolate in R^n, R^3 x S(2), or SE(3)"
        dt = vectorops.distance(a[:3],b[:3])
        if self.orientation == 'variable':
            dR = so3.distance(so3.from_moment(a[3:]),so3.from_moment(b[3:]))/math.pi
        elif self.orientation == 'axis':
            xa = from_spherical(a[3],a[4])
            xb = from_spherical(b[3],b[4])
            dR = vectorops.distance(xa,xb)
        else:
            print "Warning: don't know how to handle distance() with orientation",self.orientation
        return dt + dR

    def interpolate(self,a,b,u):
        """Interpolates in workspace between points a and b.  If len(a)<=3, 
        this will do cartesian interpolation.  Otherwise, len(a)==6 is required and
        this will do SE(3) interpolation."""
        if len(a) <= 3:
            return vectorops.interpolate(a,b,u)
        assert len(a)==self.numConstraints,"Can only interpolate in R^n, R^3 x S(2), or SE(3)"
        x =vectorops.interpolate(a[:3],b[:3],u)
        if self.orientation == 'variable':
            Ra,Rb = so3.from_moment(a[3:]),so3.from_moment(b[3:])
            R = so3.interpolate(Ra,Rb,u)
            return x + so3.moment(R)
        elif self.orientation == 'axis':
            xa = from_spherical(a[3],a[4])
            xb = from_spherical(b[3],b[4])
            #TODO: handle antipodal points?
            xu = vectorops.interpolate(xa,xb,u)
            rx = spherical(vectorops.unit(xu),radius=False)
            return x + rx
        else:
            print "Warning: don't know how to handle interpolate() with orientation",self.orientation
        return x

    def sample(self,domain):
        """Samples a point from the space of workspace parameters. domain is a 3D domain for the position component"""
        x = [random.uniform(a,b) for (a,b) in zip(bmin,bmax)]
        if self.orientation == 'variable':
            #sample rotation uniformly
            if len(x) == 3:
                x += so3.moment(so3.sample())
            else:
                #assume domain is already 6D
                assert len(x) == 6
                x[3:6] = so3.moment(so3.sample())
        elif self.orientation == 'axis':
            assert len(x) == 3
            lon = random.uniform(0,math.pi*2)
            lat = random.uniform(-math.pi/2,math.pi/2)
            x += [lon,lat]
        return x

    def grid(self,domain,celldim):
        bmin,bmax = domain
        if len(bmin) > 3:
            bmin = bmin[:3]
            bmax = bmax[:3]
        divs = [max(int(math.ceil((b-a)/celldim)),1) for (a,b) in zip(bmin,bmax)]
        Pgraph = nx.Graph()
        for cell in itertools.product(*[range(n) for n in divs]):
            params = [float(i)/float(max(div-1,1)) for i,div in zip(cell,divs)]
            x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
            Pgraph.add_node(tuple(cell),params=x)
        print "Grid discretization added",Pgraph.number_of_nodes(),"workspace nodes"
        for i in Pgraph.nodes_iter():
            n = [v for v in i]
            for dim in range(3):
                n[dim] += 1
                if tuple(n) in Pgraph.node:
                    Pgraph.add_edge(i,tuple(n))
                n[dim] -= 1
        if self.orientation == 'variable':
            print "Creating rotation graph with N=",int(math.ceil(1.0/celldim)-1)
            #Rgraph = so3_staggered_grid(int(math.ceil(1.0/celldim)))
            #create rotation graph
            Rgraph = so3_grid(int(math.ceil(1.0/celldim))-1)
            for n,d in Rgraph.nodes(data=True):
                m = so3.moment(d['params'])
                d['params'] = m
            print "  ",Rgraph.number_of_nodes(),"nodes created"
            #form cross product of self.Gw and Rgraph
            CpGraph =  nx.cartesian_product(Pgraph,Rgraph)
            print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
            return CpGraph
        if self.orientation == 'axis':
            print "Creating spherical graph with N=",int(math.ceil(1.0/celldim)-1)
            Rgraph = sphere_grid(int(math.ceil(1.0/celldim))-1)
            for n,d in Rgraph.nodes(data=True):
                m = spherical(d['params'],radius=False)
                d['params'] = m
            print "  ",Rgraph.number_of_nodes(),"nodes created"
            #form cross product of self.Gw and Rgraph
            CpGraph =  nx.cartesian_product(Pgraph,Rgraph)
            print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
            return CpGraph
        return Pgraph

    def staggered_grid(self,domain,celldim):
        bmin,bmax = domain
        if len(bmin) > 3:
            bmin = bmin[:3]
            bmax = bmax[:3]
        divs = [max(int(math.ceil((b-a)/celldim)),1) for (a,b) in zip(bmin,bmax)]
        active = [(1 if b > a else 0) for (a,b) in zip(bmin,bmax)]
        assert sum(active) >= 2,"Staggered grids don't make sense in dimension <= 1"
        Pgraph = nx.Graph()
        for cell in itertools.product(*[range(n) for n in divs]):
            params = [float(i)/float(max(div-1,1)) for i,div in zip(cell,divs)]
            x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
            Pgraph.add_node(tuple(cell),params=x,center=False)
            if all(i+a < div for (i,div,a) in zip(cell,divs,active)):
                params = [(float(i)+0.5*a)/float(max(div-1,1)) for i,div,a in zip(cell,divs,active)]
                x = [a+u*(b-a) for (a,b,u) in zip(bmin,bmax,params)]
                Pgraph.add_node(tuple([cell[0]+0.5,cell[1]+0.5,cell[2]+0.5]),params=x,center=True)
        print "Grid discretization added",Pgraph.number_of_nodes(),"workspace nodes"
        print [(b-a)/div for (a,b,div) in zip(bmin,bmax,divs)]
        celldim = max((b-a)/div for (a,b,div) in zip(bmin,bmax,divs))
        print "Connecting each workspace grid node within radius",celldim
        for i,d in Pgraph.nodes_iter(data=True):
            n = [v for v in i]
            if d['center']:
                n = [int(math.floor(v)) for v in i]
                if sum(active) >= 3:
                    #connect to other centers
                    for dim,a in enumerate(active):
                        if not a: continue
                        else:
                            n[dim] += 1
                            n[dim] += 0.5
                            if tuple(n) in Pgraph.node:
                                Pgraph.add_edge(i,tuple(n))
                            n[dim] -= 0.5
                            n[dim] -= 1
                            n[dim] = int(n[dim])
                #connect to corners
                for ofs in [(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]:
                    n = [int(math.floor(v))+ofs[dim] for dim,v in enumerate(i)]
                    if tuple(n) in Pgraph.node:
                        Pgraph.add_edge(i,tuple(n))
            else:
                for dim,a in enumerate(active):
                    if not a: continue
                    else:
                        n[dim] += 1
                        if tuple(n) in Pgraph.node:
                            Pgraph.add_edge(i,tuple(n))
                        n[dim] -= 1
        
        if self.orientation == 'variable':
            #create rotation graph
            print "Creating staggered rotation graph with N=",int(math.ceil(0.7/celldim))
            Rgraph = so3_staggered_grid(int(math.ceil(0.7/celldim)))
            for n,d in Rgraph.nodes(data=True):
                m = so3.moment(d['params'])
                d['params'] = m
            print "  ",Rgraph.number_of_nodes(),"nodes created"
            #form cross product of self.Gw and Rgraph
            CpGraph = nx.cartesian_product(Pgraph,Rgraph)
            print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
            return CpGraph
        if self.orientation == 'axis':
            print "Creating staggered spherical graph with N=",int(math.ceil(0.7/celldim)-1)
            Rgraph = sphere_staggered_grid(int(math.ceil(0.7/celldim))-1)
            for n,d in Rgraph.nodes(data=True):
                m = spherical(d['params'],radius=False)
                d['params'] = m
            print "  ",Rgraph.number_of_nodes(),"nodes created"
            #form cross product of self.Gw and Rgraph
            CpGraph =  nx.cartesian_product(Pgraph,Rgraph)
            print "Cartesian product graph has",CpGraph.number_of_nodes(),"nodes and",CpGraph.number_of_edges(),"edges"
            return CpGraph
        return Pgraph
                
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
    - ikTaskSpace: stores the IKTasks for the problem under consideration.
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
        self.ikTaskSpace = []
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
        given json object (examples of the format is given in the problems/ directory)

        Note: this is unsafe to run on untrusted inputs due to the use of eval under the
        eeorientation / fixedOrientation attributes.  To make this safe, comment out the
        line that uses eval."""
        pdef = jsonObj
        orientation = pdef.get('orientation','free')
        domain = pdef['domain']
        eelocal = pdef.get('eelocal',(0,0,0))
        useCollision = pdef.get('useCollision',True)
        link = pdef.get('link',None)
        links = pdef.get('links',None)
        eeorientation = pdef.get('eeorientation',None)
        fixedOrientation = pdef.get('fixedOrientation',None)
        if fixedOrientation is not None:
            print "Warning, fixedOrientation will be deprecated. Use eeorientation instead"
            eeorientation = fixedOrientation
        if isinstance(eeorientation,str):
            #assume some klamp't code
            eeorientation = eval(eeorientation)
        eeaxis = pdef.get('eeaxis',None)

        self.domain = domain
        #ensure domain is a bunch of floats in case theres any future divisions
        self.domain = [float(v) for v in domain[0]],[float(v) for v in domain[1]]

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

        for l in (self.ikTemplate.activeDofs if self.ikTemplate.activeDofs != None else range(self.robot.numLinks())):
            if self.robot.getJointType(l) == 'spin':
                print "NOTE: using spin joint distance calculation for walk path timing"
                self.spinJoints = True
                break

        if orientation == 'fixed':
            if eeorientation is None:
                eeorientation = link.getTransform()[0]
        if orientation == 'axis':
            if eeaxis is None:
                raise ValueError("Invalid problem setup, if orientation=axis then eeaxis must be specified")

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
            #   return False
            for x in collide.self_collision_iter(movingrobotgeoms,lambda i,j:robot.selfCollisionEnabled(movinglinks[i],movinglinks[j])):
                return False
            for x in collide.group_collision_iter(movingrobotgeoms,stillrobotgeoms,lambda i,j:robot.selfCollisionEnabled(movinglinks[i],stilllinks[j])):
                return False
            for x in collide.group_collision_iter(movingrobotgeoms,worldgeoms):
                return False
            return True

        #now set up the problem
        if orientation == 'fixed':
            self.ikTaskSpace.append(IKTask(link,eelocal,orientation=orientation,orientationparams=eeorientation))
        elif orientation == 'axis':
            self.ikTaskSpace.append(IKTask(link,eelocal,orientation=orientation,orientationparams=eeaxis))
        else:
            if orientation not in ['variable','free']:
                raise ValueError("Can only handle variable, free, fixed, or axis orientations")
            self.ikTaskSpace.append(IKTask(link,eelocal,orientation=orientation))
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
        self.ikTemplate.addObjective(self.ikTaskSpace[-1].objective)
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
        numPotentialEdges = 0
        numResolvedEdges = 0
        for i,j in self.Gw.edges_iter():
            if self.isResolvedNode(i) and self.isResolvedNode(j):
                numPotentialEdges += 1
                if self.isResolvedEdge(i,j):
                    numResolvedEdges += 1
        print "Loaded problem with",self.resolvedCount,"/",self.Gw.number_of_nodes(),"resolved workspace nodes and",numResolvedEdges,"/",numPotentialEdges,"resolved edges"
        self.buildNearestNeighbors()

    def buildNearestNeighbors(self):
        #setup nearest neighbors
        assert self.domain != None,"called buildNearestNeighbors without a domain"
        bmin,bmax = self.domain
        if len(bmax) == 6:
            self.nnw = NearestNeighbors(self.workspaceDistance,'bruteforce')
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
                sumwsdistances += self.workspaceDistance(self.Gw.node[i]['params'],self.Gw.node[j]['params'])
        if numConnected != 0:
            print "  Average configuration space distance:",sumdistances / numConnected
            print "  Average configuration space distance / workspace distance:",sumdistances / sumwsdistances

    def verify(self):
        numInvalid = 0
        numValid = 0
        for iw,jw in self.Gw.edges():
            if self.isResolvedEdge(iw,jw):
                if not self.validEdge(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
                    print "Warning, edge",iw,jw,"is marked as connected, but it is actually not valid"
                    numInvalid += 1
                    del self.Gw.edge[iw][jw]['connected']
            elif self.isResolvedNode(iw) and self.isResolvedNode(jw):
                if self.validEdge(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],self.Gw.node[iw]['params'],self.Gw.node[jw]['params']):
                    print "Warning, edge",iw,jw,"is marked as disconnected, but it is actually valid!"
                    self.markConnected(iw,jw)
                    numValid += 1
        print "Verification found",numInvalid,"/",self.Gw.number_of_edges()+numInvalid,"edges infeasible and",numValid,"edges feasible"


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

    def setIKProblem(self,x):
        """Given a workspace setting x, sets the self.ikTemplate problem so that
        it corresponds with the workspace parameters x"""
        n = 0
        for obj,task in zip(self.ikTemplate.objectives,self.ikTaskSpace):
            assert(n + task.numConstraints <= len(x)),"Invalid size of workspace parameter vector"
            task.setParams(x[n:n+task.numConstraints])
            n += task.numConstraints
        assert n == len(x),"Workspace parameters must be exactly the right number of variables in template IK constraint"

    def getIKParameters(self,q):
        """Given a configuration q, gets the workspace parameters that matches q exactly."""
        x = []
        qOrig = self.robot.getConfig()
        self.robot.setConfig(q)
        for obj,task in zip(self.ikTemplate.objectives,self.ikTaskSpace):
            link = self.robot.link(obj.link())
            x += task.getParams(link.getTransform())
        self.robot.setConfig(qOrig)
        return x

    def workspaceDistance(self,xa,xb):
        assert len(xa) == len(xb)
        d = 0
        n = 0
        for task in self.ikTaskSpace:
            assert(n + task.numConstraints <= len(xa)),"Invalid size of workspace parameter vector"
            d += task.distance(xa[n:n+task.numConstraints],xb[n:n+task.numConstraints])
            n += task.numConstraints
        return d

    def workspaceInterpolate(self,xa,xb,u):
        assert len(xa) == len(xb)
        res = []
        n = 0
        for task in self.ikTaskSpace:
            assert(n + task.numConstraints <= len(xa)),"Invalid size of workspace parameter vector"
            res.append(task.interpolate(xa[n:n+task.numConstraints],xb[n:n+task.numConstraints],u))
            n += task.numConstraints
        if len(res)==1: return res[0]
        return sum(res,[])

    def solve(self,qinit,x,testfeasibility=True):
        """Solves the IK problem at x with the starting config qinit.  Returns None if
        failed, or the solution if successful."""
        self.setIKProblem(x)
        self.robot.setConfig(qinit)
        self.ikSolverParams.startRandom = False
        if not testfeasibility:
            f,self.ikTemplate.feasibilityTest = self.ikTemplate.feasibilityTest,None
            res = self.ikTemplate.solve(self.robot,self.ikSolverParams)
            self.ikTemplate.feasibilityTest = f
            return res
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
                sortlist = [(self.workspaceDistance(pt,x),pt,n) for pt,n in knn]
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
            #   print "Warning, edge",iw,jw,"is being marked connected, but the edge is not valid"
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

    def validEdgeSimple(self,qa,qb,wa,wb,epsilon=DEFAULT_EDGE_CHECK_TOLERANCE):
        ea = self.ikError(qb,wa)
        eb = self.ikError(qa,wb)
        if self.wRadius == None:
            self.autoSetWorkspaceRadius()
        emax = self.wRadius*0.5
        #midpoint fast reject test
        x = self.robot.interpolate(qa,qb,0.5)
        if self.ikError(x,self.workspaceInterpolate(wa,wb,0.5)) > emax:
            #print "Reject deviation from line too large at midpoint"
            return False
        #if self.ikError(x,wa) > ea or self.ikError(x,wb) > eb:
        #   return False
        eaold = 0
        ebold = eb
        Ndivs = int(math.ceil(self.robot.distance(qa,qb)/epsilon))
        path = []
        for i in range(Ndivs):
            u = float(i+1)/(Ndivs+1)
            x = self.robot.interpolate(qa,qb,u)
            eu = self.ikError(x,self.workspaceInterpolate(wa,wb,u))
            eau = self.ikError(x,wa)
            ebu = self.ikError(x,wb)
            #if eau < eaold:
            #   #print "Reject monotonic increase from a at",u
            #   return False
            #if ebu > ebold:
            #   #print "Reject monotonic decrease from b at",u,",",ebold,"->",ebu
            #   return False
            if eu > emax:
                #print "Reject deviation from line too large at",u
                return False
            #if eau > ea or ebu > eb:
            #   return False
            eaold = eau
            ebold = ebu
            path.append(x)
        if self.ikTemplate.feasibilityTest != None:
            for q in path:
                if not self.ikTemplate.feasibilityTest(q):
                    return False
        return True

    def validEdgeLinear(self,qa,qb,wa,wb,epsilon=DEFAULT_EDGE_CHECK_TOLERANCE):
        #ea = self.ikError(qb,wa)
        #eb = self.ikError(qa,wb)
        qprev = qa
        #TODO determine a suitable lipschitz constant for singularity avoidance
        #right now this assumes that each DOF contributes about 1 unit to the lipschitz constant
        nlinks = (self.robot.numLinks() if self.ikTemplate.activeDofs is None else len(self.ikTemplate.activeDofs))
        epsilon *= math.sqrt(nlinks)
        Ndivs = int(math.ceil(self.robot.distance(qa,qb)/epsilon))
        discontinuityThreshold = epsilon*nlinks

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
            qm = self.solve(x,self.workspaceInterpolate(wa,wb,u),testfeasibility=False)
            if qm == None:
                #print "Unable to solve",x
                return False
            d = self.robot.distance(q0,qm)
            if d > discontinuityThreshold*(ib-ia):
                #print "Discontinuity threshold exceeded mid-segment",d,discontinuityThreshold*(ib-ia)
                return False
            d = self.robot.distance(qm,q1)
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
            q = self.solve(x,self.workspaceInterpolate(wa,wb,u),testfeasibility=False)
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
        q = self.solve(x,self.workspaceInterpolate(wa,wb,u),testfeasibility=False)
        if q is None:
            return self.robot.interpolate(qa,qb,u)
        return q

    def validEdgeBisection(self,qa,qb,wa,wb,epsilon=DEFAULT_EDGE_CHECK_TOLERANCE,c=0.9):
        d0 = self.robot.distance(qa,qb)
        if d0 <= epsilon:
            return True
        wm = self.workspaceInterpolate(wa,wb,0.5)
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
        wm = self.workspaceInterpolate(wa,wb,0.5)
        qm = self.robot.interpolate(qa,qb,0.5)
        q = self.solve(qm,wm,testfeasibility=False)
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
                x = sum([task.sample(domain) for task in self.ikTaskSpace],[])
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
            if len(self.ikTaskSpace)!=1:
                raise NotImplementedError("TODO: multiple task spaces")
            task = self.ikTaskSpace[0]

            rot_dims = 0
            if task.numConstraints > 3:
                rot_dims = task.numConstraints - 3
                bmin = bmin[:3]
                bmax = bmax[:3]

            vol = (1 if rot_dims==0 else pow(math.pi,rot_dims))
            for (a,b) in zip(bmin,bmax):
                if b > a:
                    vol *= (b-a)
            cellvol = float(vol)/N
            celldim = math.pow(cellvol,1.0/dims)
            assert celldim > 0              
            TGraph = task.grid(self.domain,celldim)
            
            #renumber, convert params from pair to 6D [translation,rotation moment] vectors
            self.Gw = nx.Graph()
            nodemap = dict()
            for n,d in TGraph.nodes(data=True):
                if len(d['params']) == 2:
                    t = d['params'][0]
                    r = d['params'][1]
                    params = t + r
                else:
                    params = d['params']
                nodemap[n] = self.Gw.number_of_nodes()
                self.Gw.add_node(self.Gw.number_of_nodes(),params = params)
            for i,j in TGraph.edges():
                self.Gw.add_edge(nodemap[i],nodemap[j])
            #redo nearest neighbors structure
            self.buildNearestNeighbors()
                
        elif method == "staggered_grid":
            if len(self.ikTaskSpace)!=1:
                raise NotImplementedError("TODO: multiple task spaces")
            task = self.ikTaskSpace[0]

            rot_dims = 0
            if task.numConstraints > 3:
                rot_dims = task.numConstraints - 3
                bmin = bmin[:3]
                bmax = bmax[:3]

            vol = (1 if rot_dims==0 else pow(math.pi,rot_dims))
            for (a,b) in zip(bmin,bmax):
                if b > a:
                    vol *= (b-a)
            #x2 to account for staggering
            cellvol = 2*float(vol)/N
            if task.numConstraints > 3:
                #account for staggering in the rotational domain
                cellvol *= pow(math.pi,task.numConstraints-3)
            celldim = math.pow(cellvol,1.0/dims)
            assert celldim > 0
            TGraph = task.staggered_grid(self.domain,celldim)
            
            #renumber, convert params from pair to 6D [translation,rotation moment] vectors
            self.Gw = nx.Graph()
            nodemap = dict()
            for n,d in TGraph.nodes(data=True):
                if len(d['params']) == 2:
                    t = d['params'][0]
                    r = d['params'][1]
                    params = t + r
                else:
                    params = d['params']
                nodemap[n] = self.Gw.number_of_nodes()
                self.Gw.add_node(self.Gw.number_of_nodes(),params = params)
            for i,j in TGraph.edges():
                self.Gw.add_edge(nodemap[i],nodemap[j])
            #redo nearest neighbors structure
            self.buildNearestNeighbors()
        else:
            raise ValueError("Unknown method "+method)
            

    def autoSetWorkspaceRadius(self):
        self.wRadius = 0
        numSame = 0
        for i,d in self.Gw.nodes_iter(data=True):
            oldw = self.wRadius
            for j in self.Gw.neighbors(i):
                dist = self.workspaceDistance(d['params'],self.Gw.node[j]['params'])
                self.wRadius = max(dist,self.wRadius)
            if self.wRadius == oldw:
                numSame += 1
                if numSame > 10:
                    #assume it's a regular grid
                    break
            else:
                numSame = 1
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
            r = max(self.workspaceDistance(x,self.Gw.node[n]['params']) for (xw,n) in nw)
        #print "Point distances", [self.workspaceDistance(self.Gw.node[n[1]]['params'],x) for n in nw]
        #print [self.Gw.node[n[1]] for n in nw]

        closest = None
        maxd = float('inf')
        for (xw,iw) in nw:
            d = self.workspaceDistance(x,self.Gw.node[iw]['params'])
            if d < maxd:
                maxd = d
                closest = iw

        #if exactly on the workspace graph, just interpolate as normal
        edgeTol = 1e-3
        Wsubgraph = self.Gw.subgraph([n[1] for n in nw])
        qclosest = None
        for (iw,jw) in Wsubgraph.edges_iter():
            xi = self.Gw.node[iw]['params']
            xj = self.Gw.node[jw]['params']
            assert len(xi) == len(x)
            (d,u) = segment_point_distance(xi,xj,x)
            #print "Edge distance", segment_point_distance(xi,xj,x)
            if self.Gw.edge[iw][jw].get('connected',False):
                if d < maxd:
                    maxd = d
                    qclosest = self.interpolateEdgeLinear(self.Gw.node[iw]['config'],self.Gw.node[jw]['config'],xi,xj,u)
        if maxd < edgeTol and qclosest is not None:
            self.lastqs = [qclosest]
            return [qclosest]

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
                distances.append(self.workspaceDistance(x,self.Gw.node[iw]['params']))
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
                if self.robot.distance(a,b) > lipschitzConstant*self.workspaceDistance(xstart,x):
                    print "Big change in configuration",i,"/",len(qs),"magnitude:",self.robot.distance(a,b),"relative to workspace distance",self.workspaceDistance(xstart,x)
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


    def closestOrientation(self,orientation):
        """For variable-rotation and axis-rotation graphs, extracts the parameters of the orientation
        in the graph closest to the input orientation.

        orientation can either be an SO3 element (9-elements) or a workspace rotation parameter vector (rotation_vector or spherical coordinates)
        """
        assert len(self.ikTaskSpace) == 1,"Can only handle one task at the moment"
        task = self.ikTaskSpace[0]
        if not hasattr(self,'allMoments'):
            self.allMoments = dict()
            for i,d in self.Gw.nodes(data=True):
                m = d['params'][3:]
                m = tuple([round(v,4) for v in m])
                if m not in self.allMoments:
                    self.allMoments[m] = []
                self.allMoments[m].append(i)
            assert len(self.allMoments) * 10 < self.Gw.number_of_nodes(),"Is this a random 6D workspace graph?"
        if len(orientation)==9:
            oparams = task.getParams((orientation,[0,0,0]))[3:]
        else:
            assert len(orientation) == task.numConstraints-3
            oparams = orientation
        mbest = None
        dbest = float('inf')
        for m in self.allMoments:
            dist = task.distance([0,0,0]+list(oparams),[0,0,0]+list(m)) 
            if dist < dbest:
                dbest = dist
                mbest = m
        return mbest

    def extractOrientedSubgraph(self,orientation):
        """For variable-rotation graphs, extracts the 3D subgraph most closely matched to 
        the given orientation parameters.  The result is (Rclosest,G) where Rclosest is the closest
        orientation parameters found and G is the 3D subgraph."""
        #find the closest orientation in the graph
        mbest = self.closestOrientation(orientation)
        subgraph = self.allMoments[mbest]
        print "Extracting rotation subgraph of size",len(subgraph)
        G = nx.subgraph(self.Gw,subgraph)
        return (mbest,G)

    def computeDiscontinuities(self,useboundary=False,orientation=None,validRangeResolution=1e-1,numOptimizationPasses=5):
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
        if orientation != None:
            print "Extracting 3D graph for orientation",orientation
            Rclosest,G = self.extractOrientedSubgraph(orientation)
        else:
            assert len(bmax) <= 3,"Can't compute discontinuity boundaries for 6D workspace graph"

        def getpoint(i):
            return G.node[i]['params'][:3]

        pts = [getpoint(i) for i in G.nodes_iter()]
        ptset = set(tuple(v) for v in pts)
        assert len(ptset) == len(pts),"Duplicate workspace points! %d unique / %d"%(len(ptset),len(pts))

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
                #disregard edges straddling boundary
                if not self.isResolvedNode(i) or not self.isResolvedNode(j):
                    continue
            else:
                #keep edges straddling boundary but not entirely outside
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
                if self.isResolvedNode(i) != self.isResolvedNode(j):
                    #external discontinuity: interpolate from i toward j (or j toward i)
                    flip = self.isResolvedNode(j)
                    if flip:
                        q0 = G.node[j]['config']
                        a,b = b,a
                    else:
                        q0 = G.node[i]['config']
                    l,u = 0.0,1.0
                    while u > l+validRangeResolution:
                        m = 0.5*(l+u)
                        x = self.workspaceInterpolate(a,b,m)
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
                        while u > l+validRangeResolution:
                            m = 0.5*(l+u)
                            x = self.workspaceInterpolate(a,b,m)
                            if self.solve(q0,x) is not None:
                                l = m
                            else:
                                u = m
                        max0 = l
                    if self.solve(q1,a):
                        min1 = 0.0
                    else:
                        l,u = 0.0,1.0
                        while u > l+validRangeResolution:
                            m = 0.5*(l+u)
                            x = self.workspaceInterpolate(a,b,m)
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
        ptset = set(tuple(v) for v in pts)
        assert len(ptset) == len(pts),"Duplicate points! %d unique / %d"%(len(ptset),len(pts))

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

        #decompose the vertices + edge midpoints + centroid set into simplices...
        #we can potentially draw faces of each simplex 
        hpts = [[x[a] for a in active] for x in pts]

        try:
            delaunay = scipy.spatial.Delaunay(hpts)
        except Exception as e:
            print "Error computing Delaunay triangulation"
            print e
            return None
        faces = set()
        """
        #UNCOMMENT ME FOR DEBUGGING IN VISUALIZATION
        self.delaunay = delaunay
        self.delaunay_pts = pts
        self.delaunay_epts = epts
        self.delaunay_others = others
        """
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
                    #exactly two of the points must be a simplex point
                    if sum([1 for a in tri if a in others]) != 2:
                        continue
                    faces.add(tuple(sorted(tri)))
            else:
                for i in range(len(simplex)):
                    n = (i+1)%len(simplex)
                    a,b = simplex[i],simplex[n]
                    if not (a in epts or b in epts):
                        continue
                    #one of the points must be a simplex point
                    if not (a in others or b in others):
                        continue
                    faces.add(tuple(sorted((a,b))))
        faces = list(faces)
        mesh = MeshTopology(hpts,faces)

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
                    #   if sum(1 for u in adjacentConflicts if u in incidentVerts[w])==1:
                    #       collapseto = w
                    #       break
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

        #optimize the coordinates of epts so that they 1) lie in their specified range and 2) minimize curvature
        attached = dict()
        for v in others:
            if v not in mesh.incidentVerts: continue
            #attach centroid points to single conflict points
            adjacentConflicts = [w for w in set(mesh.incidentVerts[v]) if w in epts]
            if len(adjacentConflicts) == 1:
                if adjacentConflicts[0] not in attached:
                    attached[adjacentConflicts[0]] = []
                attached[adjacentConflicts[0]].append(v)

        for sweepCnt in xrange(numOptimizationPasses):
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
                #ensure consistent face orientation... is this necessary?
                n = vectorops.cross(vectorops.sub(verts[fmap[1]],verts[fmap[0]]),vectorops.sub(verts[fmap[2]],verts[fmap[0]]))
                if n[2] < 0:
                    fmap[1],fmap[2] = fmap[2],fmap[1]
            fmapped.append(fmap)
        return verts,fmapped

