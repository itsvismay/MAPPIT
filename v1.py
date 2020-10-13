import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np 
from scipy.stats import norm, gaussian_kde
from mpl_toolkits.mplot3d import Axes3D

delta = 0.1
dt = 0.1
population = np.zeros([100,2])
population += norm.rvs(size = population.shape, scale=delta**(2*math.sqrt(dt)))
velocities = np.zeros(population.shape)
pop_mean = np.array([0,0])
mousex = 0 
mousey = 0

def on_mouse_move(event):
	mousex = event.xdata
	mousey = event.ydata


def update_population(population, velocities, pop_mean, mx, my):	
	#brownian motion
	delta = 10
	dt = 0.1
	velocities = norm.rvs(size = population.shape, scale=delta**(2*math.sqrt(dt)))

	#Liousville equation
	if(np.linalg.norm(population)==0):
		maxi = 0
		m_t = np.zeros(population.shape[0])
	else:
		kernel = gaussian_kde(population.T)
		m_t = kernel(population.T)
		maxi = np.argmax(m_t)

	pop_mean = population[maxi]

	# Update Individual
	for i in range(population.shape[0]):
		#dist from pop mean
		dir_d1 = pop_mean - population[i,:]
		d1 = np.linalg.norm(dir_d1)

		#dist from nest (0,0)
		dir_d2 = np.array([mx, my]) - population[i,:]
		d2 = np.linalg.norm(dir_d2)

		v_n = dir_d1 + 200*dir_d2
	
		velocities[i,:] += v_n[0]


	population = population + dt*velocities

	return m_t, population

z = update_population(population, velocities, pop_mean, mousex, mousey)


#setup animation
plt.style.use('seaborn-pastel')
fig = plt.figure()
ax = plt.axes(xlim=(-5, 5), ylim=(-5, 5))
dots, = ax.plot([], [], 'go')

def init():
    dots.set_data([], [])
    return dots,

def animate(i, population, velocities, pop_mean, mx, my):
	if(i==0):
		population*=0

	interval = 360/1000.0
	mx = math.cos(i*interval)*0.05
	my = math.sin(i*interval)*0.05
	print(mx, my)

	z,pop = update_population(population, velocities, pop_mean, mx, my)

	dots.set_data(pop[:,0], pop[:,1])
	# RGBz = np.multiply(RGB, z[:,np.newaxis])
	# dots.set_color(RGBz)
	return dots,

RGB = np.ones((population.shape[0], 3))
anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(population, velocities, pop_mean, mousex, mousey),
                               frames=1000, interval=50, blit=True)
plt.connect('motion_notify_event',on_mouse_move)
plt.show()

#setup 3d animation (to see the kernel function)
# def update_graph(i, population):
# 	if(i==0):
# 		population*=0
# 	print(i)
# 	z = update_population(population)
# 	graph._offsets3d = (population[:,0], population[:,1], z)
# 	text.set_text("{:d}: ".format(i)) 
	

# fig = plt.figure()
# ax = Axes3D(fig)
# title = ax.set_title('3D Test')
# text = fig.text(0, 1, "TEXT", va='top') 

# graph = ax.scatter(population[:,0], population[:,1], z)
# ani = animation.FuncAnimation(fig, update_graph, fargs=(population,), frames=20, interval=500, blit=False)
# plt.show()