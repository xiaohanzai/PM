import numpy as np

N = 20000
fname = 'clump_0.txt'
G = 1
Mhalo = 20000
rhalo = 1

# parameters for the clump
N1 = 200
xc = 0.8
yc = 0
m_particle = 1
radius = 0.1
nsigma = 0.068

# generate positions
X = np.random.rand(2*N1,3)*2 - 1
ind = np.where(np.sum(X**2,axis=1)<=1)[0]
x = X[ind[0:N1],0]*radius + xc
y = X[ind[0:N1],1]*radius + yc
r = (x**2 + y**2)**0.5

# generate velocities
omega = (3*np.pi*G*N/4 + Mhalo/rhalo**3)**0.5
v_rot = r*omega
sigma = nsigma * 1*omega
v_r = np.random.randn(N1)*sigma
v_theta = np.random.randn(N1)*sigma + v_rot
v_x = v_r*x/r - v_theta*y/r
v_y = v_r*y/r + v_theta*x/r

# generate masses
m = np.ones(N1)*m_particle

# write to file
f = open(fname,'w')
for i in range(N1):
	f.write('%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % 
		(x[i], y[i], v_x[i], v_y[i], m[i], r[i], v_r[i], v_theta[i]))
f.close()
