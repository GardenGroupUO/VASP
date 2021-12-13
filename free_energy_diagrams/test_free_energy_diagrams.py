import matplotlib.pyplot as plt
from free_energy_diagrams import FED


au_hey = FED()
au_hey.add_level(0, 'Clean', top_text='0', color='k')

au_hey.add_level(0.48, color='k')
au_hey.add_barrier(1.22, 2, 4, ls='-', linewidth=1, color='k')
au_hey.add_level(0.3, bottom_text='$H_{s}$', position='l', color='#800080')
au_hey.add_barrier(1.24, 3, 5, ls='--', linewidth=1, color='#800080')


au_hey.add_level(0.88, bottom_text='$H_{Mo}$', color='k')
au_hey.add_barrier(1.35, 4, 6, ls='-', linewidth=1, color='k', barrier_text='1.35')
au_hey.add_level(1.05, position='l', color='#800080')
au_hey.add_barrier(1.31, 5, 6, ls='--', linewidth=1, color='#800080', barrier_text='1.31')

au_hey.add_level(0, 'Clean', top_text='0', color='k')

au_hey.plot()
plt.show()