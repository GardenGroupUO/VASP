import matplotlib.pyplot as plt
from free_energy_diagrams import FED


au_hey = FED()
au_hey.add_level(0, 'Clean', top_text='0', color='k')

au_hey.add_level(0.48, color='#0011FF', label='Ciaran')
au_hey.add_level(0.3, bottom_text='$H_{s}$', position='l', color='#FF9900', label='Charlie')

au_hey.add_level(0.88, bottom_text='$H_{Mo}$', color='#0011FF')
au_hey.add_level(1.05, position='l', color='#FF9900')

au_hey.add_level(0, 'Clean', top_text='0', color='k')

au_hey.add_barrier(1.22, 2, 4, ls='-',  color='#0011FF')
au_hey.add_barrier(1.24, 3, 5, ls='-', color='#FF9900')
au_hey.add_barrier(1.35, 4, 6, ls='-',  color='#0011FF')
au_hey.add_barrier(1.31, 5, 6, ls='-', color='#FF9900')

au_hey.add_link(0, 1, color='#0011FF', ls='--', linewidth=1.5)
au_hey.add_link(0, 2, color='#FF9900', ls='--', linewidth=1.5)

au_hey.plot()
plt.show()