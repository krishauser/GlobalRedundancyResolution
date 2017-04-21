import sys

n = int(sys.argv[1])

#scale = 1.0/n
scale = 1.0;

intro = """### A %dR planar robot ###
TParent 1 0 0   0 1 0   0 0 1   0 0 0 """ %(n,)
trans = """   \\\n1 0 0   0 1 0   0 0 1   %g 0 0"""%(scale,)
axis = "0 1 0"
jointtype = "r"
qMin = "0"
qMax = "6.28319"
q ="0"
geometry = ' "objects/thincube.tri"'
geomscale = "%g"%(max(scale,0.05))
mass =  "1"
torqueMax = "1"
accMax = "1"
velMax = "1"

joint = "joint spin"
servoP = "50"
servoI = "100"
servoD = "3"

print intro,
for i in xrange(n-1):
    print trans,
print
print "axis\t"+"\t".join([axis]*n)
print "jointtype\t"+" ".join([jointtype]*n)
print "qMin\t"+" ".join([qMin]*n)
print "qMax\t"+" ".join([qMax]*n)
print "q\t"+" ".join([q]*n)
print "geometry\t"+" ".join([geometry]*n)
print "geomscale\t"+" ".join([geomscale]*n)
print "mass\t"+" ".join([mass]*n)
print "automass"
print "torqueMax\t"+" ".join([torqueMax]*n)
print "accMax\t"+" ".join([accMax]*n)
print "velMax\t"+" ".join([velMax]*n)
for i in xrange(n):
    print joint,i
print "servoP\t"+" ".join([servoP]*n)
print "servoI\t"+" ".join([servoI]*n)
print "servoD\t"+" ".join([servoD]*n)

