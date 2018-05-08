import numpy as np

N = 20000
fname = 'ic_0.txt'
G = 1
halo = False # with or without halo

# generate positions
X = np.random.rand(2*N,3)*2 - 1
ind = np.where(np.sum(X**2,axis=1)<=1)[0]
x = X[ind[0:N],0]
y = X[ind[0:N],1]
r = (x**2 + y**2)**0.5

# generate velocities
omega = (3*np.pi*G*N/4 + halo*N)**0.5
v_rot = r*omega
# use a fraction of v(R_0) (for instace 13.6%) as the velocity dispersion
sigma = 0 #0.136 * 1*omega
# # or use equation 15 of Hockney & Hohl 1969
# sigma = 0.342*(1-r**2)**0.5*omega
v_r = np.random.randn(N)*sigma
v_theta = np.random.randn(N)*sigma + v_rot
v_x = v_r*x/r - v_theta*y/r
v_y = v_r*y/r + v_theta*x/r

# generate masses
m = np.ones(N) # temporarily use mass = 1 for all particles

# write to file
f = open(fname,'w')
for i in range(N):
	f.write('%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % 
		(x[i], y[i], v_x[i], v_y[i], m[i], r[i], v_r[i], v_theta[i]))
f.close()
