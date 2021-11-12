#!/usr/bin/env python

import matplotlib.pyplot as plt
import pdb
import numpy as np
from matplotlib.lines import Line2D as line
from scipy.interpolate import interp1d

class Free_Energy_Diagram():
    def __init__(self,  ):












state_len = 0.3
state_pos = 0

#vol
hey = 0.29
diff = 0.74

vaspsol_E = [0, 0.48, 0.88, 0]

# Can make lines by inserting x and y values as lists.



def get_xspace(y_data):
	state_bounds = []
	line_bounds = []
	zero_pos = 0
	sections = np.arange(0, (len(y_data)*2)+1, 1)
	for i in sections:
		if (i % 2) == 0:
			line_bounds.append([zero_pos, zero_pos+0.3])
			zero_pos += 0.3
		else:
			state_bounds.append([zero_pos, zero_pos+0.3])
			zero_pos += 0.3
	return [state_bounds, line_bounds]



fig, ax = plt.subplots(1, 1)
vaspsol_xspace = get_xspace(vaspsol_E)

def get_barrier_coords(xspace_pos, barrier, yspace_pos):
	x_ran_mid = xspace_pos[0] + (xspace_pos[1]-xspace_pos[0])/2
	xspace_pos.insert(1, x_ran_mid)

	barrier_height = yspace_pos[0] + barrier
	yspace_pos.insert(1, barrier_height)
	return [xspace_pos, yspace_pos]



hey_yran = [0.88, 0]
diff_yran = [0.48, 0.88]

hey_bounds = get_barrier_coords(vaspsol_xspace[1][3], hey, hey_yran)
diff_bounds = get_barrier_coords(vaspsol_xspace[1][2], diff, diff_yran)

def get_barrier_spline(barrier_xran, barrier_yran):
	inter = interp1d(barrier_xran, barrier_yran, kind='quadratic')
	newx_ran = np.linspace(barrier_xran[0], barrier_xran[2], 1000)
	inter_y = inter(newx_ran)
	return [newx_ran, inter_y]

hey_inter = get_barrier_spline(hey_bounds[0], hey_bounds[1])
diff_inter = get_barrier_spline(diff_bounds[0], diff_bounds[1])

plt.plot(hey_inter[0], hey_inter[1], color='#ff8900')
plt.plot(diff_inter[0], diff_inter[1], color='#ff8900')


#print(vasp_xspace)

#ax.plot([che_xspace[0][0][0],che_xspace[0][0][1]], [che_E[0], che_E[0]], color='b', label='CHE')
# ax.plot([che_xspace[0][1][0], che_xspace[0][1][1]], [che_E[1], che_E[1]], linewidth=2.5, color='#0070c0', label='CHE')
# ax.plot([che_xspace[1][1][0], che_xspace[1][1][1]], [che_E[0], che_E[1]], color='#0070c0', ls='--')

# ax.plot([che_xspace[0][2][0], che_xspace[0][2][1]], [che_E[2], che_E[2]], linewidth=2.5, color='#0070c0')
# ax.plot([che_xspace[1][2][0], che_xspace[1][2][1]], [che_E[1], che_E[2]], color='#0070c0', ls='--')

# ax.plot([che_xspace[1][3][0], che_xspace[1][3][1]], [che_E[2], che_E[3]], color='#0070c0', ls='--')


ax.plot([vaspsol_xspace[0][0][0], vaspsol_xspace[0][0][1]], [vaspsol_E[0], vaspsol_E[0]], linewidth=2.5, color='k')

ax.plot([vaspsol_xspace[0][1][0], vaspsol_xspace[0][1][1]], [vaspsol_E[1], vaspsol_E[1]], linewidth=2.5, color='#ff8900')

#ax.plot([vaspsol_xspace[1][1][0], vaspsol_xspace[1][1][1]], [vaspsol_E[0], vaspsol_E[1]], color='#ff8900', ls='--')

ax.plot([vaspsol_xspace[0][2][0], vaspsol_xspace[0][2][1]], [vaspsol_E[2], vaspsol_E[2]], linewidth=2.5, color='#ff8900')
#ax.plot(hey_bounds[0], hey_bounds[1], 'o', diff_xnew, her_inter(diff_xnew), '-')

ax.plot([vaspsol_xspace[0][3][0], vaspsol_xspace[0][3][1]], [vaspsol_E[3], vaspsol_E[3]], linewidth=2.5, color='k')

#ax.plot([vaspsol_xspace[1][3][0], vaspsol_xspace[1][3][1]], [vaspsol_E[2], vaspsol_E[3]], color='#ff8900', ls='--')

# for i in np.arange(0, len(che_E), 1):
# 	che_x = che_xspace[0][i]
# 	che_y = che_E[i]
# 	ax.add_artist(line([che_x[0], che_x[1]],[che_y, che_y], linewidth=2, color='b', label='CHE'))

# for i in np.arange(0, len(vaspsol_E), 1):
# 	sol_x = vaspsol_xspace[0][i]
# 	sol_y = vaspsol_E[i]
# 	ax.add_artist(line([sol_x[0], sol_x[1]],[sol_y, sol_y], linewidth=2, color='k', label='VASPsol'))



# min1 = che_E[0]
# max1 = che_E[0] 
# for i in range(len(che_E)):
# 	if che_E[i] < min1:
# 		min1 = che_E[i]
# 	if che_E[i] > max1:
# 		max1 = che_E[i]

ticks = ['Clean', '$H_{s}$', '$H_{Mo}$', 'Clean']
tick_pos = []


for i in np.arange(0, 4, 1):
	var = vaspsol_xspace[0][i][0]+ 0.15
	tick_pos.append(var)


ax.set_ylim(-0.1, 1.9)
ax.set_xlim(0, vaspsol_xspace[1][-1][1])
#ax.legend()
plt.xticks(tick_pos, ticks)
plt.xlabel('Reaction Coordinate')
plt.ylabel('Gibbs Free Energy (eV)')
plt.show()





