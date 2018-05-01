import numpy as np

N = 20000
fname = 'ic_0_3d.txt'
G = 1

# generate positions
X = np.random.rand(2*N,3)*2 - 1
ind = np.where(np.sum(X**2,axis=1)<=1)[0]
x = X[ind[0:N],0]
y = X[ind[0:N],1]
r = (x**2 + y**2)**0.5

# generate velocities
omega = (3*np.pi*G*N/4)**0.5
v_rot = r*omega
sigma = 0 #0.136 * 1*omega # use 13.6% v(R_0) as the velocity dispersion
v_r = np.random.randn(N)*sigma
v_theta = np.random.randn(N)*sigma + v_rot
v_x = v_r*x/r - v_theta*y/r
v_y = v_r*y/r + v_theta*x/r

# generate masses
m = np.ones(N) # temporarily use mass = 1 for all particles

# write to file
f = open(fname,'w')
for i in range(N):
	f.write('%.6f %.6f 0 %.6f %.6f 0 %.6f\n' % 
		(x[i], y[i], v_x[i], v_y[i], m[i]))
f.close()
