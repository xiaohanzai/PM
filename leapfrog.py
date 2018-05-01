import numpy as np
import matplotlib.pyplot as plt

x1 = 1
y1 = 0
vx1 = 0
vy1 = 1/2

x2 = -1
y2 = 0
vx2 = 0
vy2 = -1/2

T = 2*np.pi*2**1.5
dt = T/100

fig = plt.figure()
ax = fig.add_subplot(111)
r = ((x1 - x2)**2 + (y1 - y2)**2)**0.5
# ax.plot(0, 0.5*(vx1**2 + vx2**2 + vy1**2 + vy2**2) - 1/r,'k.')

for i in range(100):
	# kick
	r = ((x1 - x2)**2 + (y1 - y2)**2)**0.5
	ax1 = (x2 - x1)/r/r**2
	ay1 = (y2 - y1)/r/r**2
	ax2 = -ax1
	ay2 = -ay1

	vx1h = vx1 + ax1*dt/2.
	vy1h = vy1 + ay1*dt/2.
	vx2h = vx2 + ax2*dt/2.
	vy2h = vy2 + ay2*dt/2.

	# drift
	x1 += vx1h*dt
	y1 += vy1h*dt
	x2 += vx2h*dt
	y2 += vy2h*dt

	# kick
	r = ((x1 - x2)**2 + (y1 - y2)**2)**0.5
	ax1 = (x2 - x1)/r/r**2
	ay1 = (y2 - y1)/r/r**2
	ax2 = -ax1
	ay2 = -ay1

	vx1 = vx1h + ax1*dt/2.
	vy1 = vy1h + ay1*dt/2.
	vx2 = vx2h + ax2*dt/2.
	vy2 = vy2h + ay2*dt/2.

	ax.plot(x1, y1, 'k.')
	ax.plot(x2, y2, 'r.')
	# ax.plot((i+1)*dt, 0.5*(vx1**2 + vx2**2 + vy1**2 + vy2**2)-1/r,'k.')

plt.show()






