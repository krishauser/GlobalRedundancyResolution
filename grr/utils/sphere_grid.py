from klampt.math import vectorops
import math
import networkx as nx

def sphere_grid(N):
    """Returns a nx.Graph describing a grid on S2 with resolution N. 
    The resulting graph has ~6N^2 nodes

    The node attribute 'params' gives the point.
    """
    assert N >= 1
    G = nx.Graph()
    namemap = dict()
    for ax in range(3):
        for i in xrange(N+1):
            u = 2*float(i)/N - 1.0
            u = math.tan(u*math.pi*0.25)
            for j in xrange(N+1):
                v = 2*float(j)/N - 1.0
                v = math.tan(v*math.pi*0.25)
                idx = [i,j]
                idx = idx[:ax] + [N] + idx[ax:]
                idx = tuple(idx)
                if idx in namemap:
                    continue
                else:
                    namemap[idx] = idx
                    x = [u,v]
                    x = x[:ax] + [1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(idx,params=x)
                idx = [i,j]
                idx = idx[:ax] + [0] + idx[ax:]
                idx = tuple(idx)
                if idx in namemap:
                    continue
                else:
                    namemap[idx] = idx
                    x = [u,v]
                    x = x[:ax] + [-1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(idx,params=x)
    for ax in range(3):
        for i in xrange(N+1):
            for j in xrange(N+1):
                baseidx = [i,j]
                idx = baseidx[:ax] + [N] + baseidx[ax:]
                idx = tuple(idx)
                for n in range(2):
                    nidx = baseidx[:]
                    nidx[n] += 1
                    if nidx[n] > N: continue
                    nidx = nidx[:ax] + [N] + nidx[ax:]
                    nidx = tuple(nidx)
                    G.add_edge(namemap[idx],namemap[nidx])
                idx = baseidx[:ax] + [0] + baseidx[ax:]
                idx = tuple(idx)
                for n in range(2):
                    nidx = baseidx[:]
                    nidx[n] += 1
                    if nidx[n] > N: continue
                    nidx = nidx[:ax] + [0] + nidx[ax:]
                    nidx = tuple(nidx)
                    G.add_edge(namemap[idx],namemap[nidx])
    return G

def sphere_staggered_grid(N):
    """Returns a nx.Graph describing a staggered grid on S2 with resolution N. 
    The resulting graph has ~12N^2 nodes

    The node attribute 'params' gives the so3 element.
    """
    import itertools
    assert N >= 1
    G = nx.Graph()
    namemap = dict()
    for ax in range(3):
        for i in xrange(N+1):
            u = 2*float(i)/N - 1.0
            u = math.tan(u*math.pi*0.25)
            u2 = 2*float(i + 0.5)/N - 1.0
            u2 = math.tan(u2*math.pi*0.25)
            for j in xrange(N+1):
                v = 2*float(j)/N - 1.0
                v = math.tan(v*math.pi*0.25)
                v2 = 2*float(j+ 0.5)/N - 1.0
                v2 = math.tan(v2*math.pi*0.25)
                
                idx = [i,j]
                idx = idx[:ax] + [N] + idx[ax:]
                idx = tuple(idx)
                if idx in namemap:
                    pass
                else:
                    namemap[idx] = idx
                    x = [u,v]
                    x = x[:ax] + [1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(idx,params=x)

                if i < N and j < N:
                    #add staggered point
                    sidx = [i+0.5,j+0.5]
                    sidx = sidx[:ax] + [N] + sidx[ax:]
                    sidx = tuple(sidx)
                    namemap[sidx] = sidx
                    x = [u2,v2]
                    x = x[:ax] + [1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(sidx,params=x)

                idx = [i,j]
                idx = idx[:ax] + [0] + idx[ax:]
                idx = tuple(idx)
                if idx in namemap:
                    pass
                else:
                    namemap[idx] = idx
                    x = [u,v]
                    x = x[:ax] + [-1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(idx,params=x)

                if i < N and j < N:
                    #add staggered point
                    sidx = [i+0.5,j+0.5]
                    sidx = sidx[:ax] + [0] + sidx[ax:]
                    sidx = tuple(sidx)
                    namemap[sidx] = sidx
                    x = [u2,v2]
                    x = x[:ax] + [-1.0] + x[ax:]
                    x = vectorops.unit(x)
                    G.add_node(sidx,params=x)
    for ax in range(3):
        for i in xrange(N+1):
            for j in xrange(N+1):
                baseidx = [i,j]
                idx = baseidx[:ax] + [N] + baseidx[ax:]
                idx = tuple(idx)
                for n in range(2):
                    nidx = baseidx[:]
                    nidx[n] += 1
                    if nidx[n] > N: continue
                    nidx = nidx[:ax] + [N] + nidx[ax:]
                    nidx = tuple(nidx)
                    G.add_edge(namemap[idx],namemap[nidx])

                idx = baseidx[:ax] + [0] + baseidx[ax:]
                idx = tuple(idx)
                for n in range(2):
                    nidx = baseidx[:]
                    nidx[n] += 1
                    if nidx[n] > N: continue
                    nidx = nidx[:ax] + [0] + nidx[ax:]
                    nidx = tuple(nidx)
                    G.add_edge(namemap[idx],namemap[nidx])

                #edges between staggered point and grid points
                baseidx = [i+0.5,j+0.5]
                if baseidx[0] > N or baseidx[1] > N:
                    continue

                idx = baseidx[:ax] + [N] + baseidx[ax:]
                idx = tuple(idx)
                for ofs in itertools.product(*[(-0.5,0.5)]*2):
                    nidx = vectorops.add(baseidx,ofs)
                    nidx = [int(x) for x in nidx]
                    nidx = nidx[:ax] + [N] + nidx[ax:]
                    nidx = tuple(nidx)
                    #print "Stagger-grid edge",idx,nidx
                    G.add_edge(namemap[idx],namemap[nidx])

                
                idx = baseidx[:ax] + [0] + baseidx[ax:]
                idx = tuple(idx)
                for ofs in itertools.product(*[(-0.5,0.5)]*2):
                    nidx = vectorops.add(baseidx,ofs)
                    nidx = [int(x) for x in nidx]
                    nidx = nidx[:ax] + [0] + nidx[ax:]
                    nidx = tuple(nidx)
                    #print "Stagger-grid edge",idx,nidx
                    G.add_edge(namemap[idx],namemap[nidx])
    return G




def sphere_grid_test(N=5,staggered=True):
    from klampt import vis
    from klampt.model import trajectory
    if staggered:
        G = sphere_staggered_grid(N)
    else:
        G = sphere_grid(N)
    for n in G.nodes():
        x = G.node[n]['params']
        vis.add(str(n),x)
        vis.hideLabel(str(n))
    #draw edges?
    minDist = float('inf')
    maxDist = 0.0
    for i,j in G.edges():
        xi = G.node[i]['params']
        xj = G.node[j]['params']
        vis.add(str(i)+'-'+str(j),trajectory.Trajectory([0,1],[xi,xj]))
        vis.hideLabel(str(i)+'-'+str(j))
        dist = vectorops.distance(xi,xj)
        if dist > maxDist:
            maxDist = dist
            print "Max dispersion at",i,j,":",dist
        if dist < minDist:
            minDist = dist
    print "Number of points:",G.number_of_nodes()
    print "Min/max dispersion:",minDist,maxDist
    vis.run()

if __name__ == '__main__':
    sphere_grid_test(5,True)