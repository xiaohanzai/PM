import numpy as np

N = 20000
fname = 'ic_0.txt'

# generate positions
icdf = lambda F: 1*np.sqrt(1 - (1-F)**(2/3))
sample = lambda n: icdf(np.random.random(n))
r = sample(N)
theta = np.random.random(N)*2*np.pi
x = r*np.cos(theta)
y = r*np.sin(theta)

# generate velocities
omega = (3*np.pi*N/4)**0.5
v_rot = r*omega
sigma = 0 #0.136 * 1*omega # use 13.6% v(R_0) as the velocity dispersion
v_r = np.random.randn(N)*sigma
v_theta = np.random.randn(N)*sigma + v_rot
v_x = v_r*np.cos(theta) - v_theta*np.sin(theta)
v_y = v_r*np.sin(theta) + v_theta*np.cos(theta)

# generate masses
m = np.ones(N) # temporarily use mass = 1 for all particles

# write to file
f = open(fname,'w')
for i in range(N):
	f.write('%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % 
		(x[i], y[i], v_x[i], v_y[i], m[i], r[i], v_r[i], v_theta[i]))
f.close()
