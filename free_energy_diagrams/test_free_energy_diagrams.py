import matplotlib.pyplot as plt
from free_energy_diagrams import FED


au_hey = FED()
au_hey.add_level(0, 'Clean', top_text='0', color='k')
au_hey.add_level(0.48, '$H_{s}$', top_text='0.48', color='k')
au_hey.add_level(0.88, '$H_{Mo}$', top_text='0.88', color='k')
au_hey.add_level(0, 'Clean', top_text='0', color='k')

au_hey.add_barrier(1.22, 2, 3, ls='-', linewidth=0.5, color='k')
au_hey.add_barrier(1.35, 3, 4, ls='-', linewidth=0.5, color='k')

au_hey.plot()
plt.show()