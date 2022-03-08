#!/usr/bin/env python
""" 
This module is designed to produce free energy diagrams given a set of
thermodynamic energy levels and kinetic barriers.

It was heavily based off a script written by Giacomo Marchioro called
PyEnergyDiagrams at their GitHub: giacomomarchioro.

"""

import matplotlib.pyplot as plt
import pdb
import numpy as np
from matplotlib.lines import Line2D




class FED():
	
    golden_ratio = 1.6181
    offset_ratio = 0.02

    def __init__(self,
                reaction_coords = [],
                aspect='equal',
                level_width=3,
                barrier_width=1.5,
                dimension = 'auto',
                space = 'auto',
                offset = 'auto',
                ylabel = 'Energy (Ev)',
                xlabel = 'Reaction Progress'):

        self.reaction_coords = reaction_coords
        self.dimension = dimension
        self.space = space
        self.offset = offset
        self.aspect = aspect
        self.level_width = level_width
        self.barrier_width = barrier_width
        self.xlabel = xlabel
        self.ylabel = ylabel

        self.pos_number = 0
        self.energies = []
        self.positions = []
        self.bottom_texts = []
        self.top_texts = []
        self.left_texts = []
        self.right_texts = []
        self.colours = []
        self.links = []
        self.barriers = []
        self.legends = []
       

        self.ax = None
        self.fig = None
        
        

    def add_level(self, energy, bottom_text='', position=None, top_text='',
                left_text='', right_text='', color='k', label=''):
        '''
        This is a method of the FED class. It will save all of the data
        for a new energy level into the lists required.

        Parameters
        ----------
        energy : int
            The energy of the level in eV
        xtick_label : str
            The text on the x axis tick to label the level.
        bottom_text : str
            The text underneath the level, ie the label. (default = '')
        position : str
            The position of the level in the plot. If it is left empty then the
            level will be added after the last energy level. If position is last
            or l, then the level will be added bellow the previous energy level.
        top_text : str
            Text above the energy level, defaults to 'Energy', the energy of the
            level. (default = 'Energy')
        left_text : str
            The text to the left of the level. (default = '')
        right_text : str
            The text to the right of the level. (default = '')
        color : str
            Colour of the level. (default = 'k')

        Returns
        -------
        Append all the data about the level in question to the relevant
        class attributes. 

        '''

        if position is None:
            position = self.pos_number + 1
            self.pos_number += 1
        elif position == 'last' or position == 'l':
            position = self.pos_number
        else:
            raise ValueError(
                "Position must be None or 'last' (abrv. 'l'. It was: %s" % position
            )
        
        #if top_text == 'Energy':
        #    top_text = energy
        
        link = []
        
        self.energies.append(energy)
        self.positions.append(position)
        self.bottom_texts.append(bottom_text)
        self.top_texts.append(top_text)
        self.left_texts.append(left_text)
        self.right_texts.append(right_text)
        self.colours.append(color)
        self.links.append(link)
        self.legends.append(label)

    def add_link(self, start_level_id, end_level_id, color = 'k',
             ls='--', linewidth=1):
        '''
        This is a method of the FED class that will take in a start and end
        level id (start/end_level_id) and create a linker between those
        two levels. This will be used for levels that do not have a barrier
        between them.

        Parameters
        ----------
        start_level_id : int
            The id (index of the level) of the level that you want the
            link to start from.
        end_level_id : int
            The id (index of the level) of the level that you want the
            link to end.
        color : str
            The colour of the link.
        ls : str
            The link style.
        linewidth : int
            The width of the link.

        Returns
        -------
        Appends a tuple containin the imformation about the linker to a list in 
        the class links attribute.
        '''

        self.links[start_level_id].append((end_level_id, ls, linewidth, color))

    def auto_adjust(self):
        '''
        This method of the FED class will use the golden ratio to set the 
        best dimension and space between the levels.

        Affects
        
        self.dimension
        self.space
        self.offset
        '''
        # Max range between the energy
        Energy_variation = abs(max(self.energies) - min(self.energies))
        if self.dimension == 'auto' or self.space == 'auto':
            # Unique positions of the levels
            unique_positions = float(len(set(self.positions)))
            space_for_level = Energy_variation*FED.golden_ratio/unique_positions
            self.dimension = space_for_level*0.7
            self.space = space_for_level*0.3

        if self.offset == 'auto':
            self.offset = Energy_variation*FED.offset_ratio



    def add_barrier(self, energy, start_level_id, end_level_id,
                    ls='-', color='k'):
        '''
        This is a method of the FED class that will take in a start/end level id 
        and an energy and create a barrier between those two levels. 

        Parameters
        ----------
        start_level_id : int
            The id (index of the level) of the level that you want the
            barrier to start from.
        end_level_id : int
            The id (index of the level) of the level that you want the
            barrier to end.
        energy : int
            The energy (saddle point) of the barrier
        ls : str
            The barrier style. (Default = '-')
        linewidth : int
            The width of the link.
        color : str
            The colour of the barrier. (Default = 'k')

        Returns
        -------
        Appends a tuple containing the information about the barrier to a
        list in the class barriers attribute. 
        '''

        self.auto_adjust()

        start = self.positions[start_level_id-1]*(self.dimension+self.space)
        x1 = start + self.dimension
        y1 = self.energies[start_level_id-1]
        x2 = self.positions[end_level_id-1]*(self.dimension+self.space)
        y2 = self.energies[end_level_id-1]
        vert_x = x1 + ((x2-x1)/2)
        vert_y = energy
        a1 = (y1 - vert_y)/((x1 - vert_x)**2)
        a2 = (y2 - vert_y)/((x2 - vert_x)**2)
        left_xspace = list(np.linspace(x1, vert_x, 500))
        right_xspace = list(np.linspace(vert_x, x2, 500))
        left_yspace = []
        right_yspace = []
        for y in left_xspace:
            left_y = (a1*((y - vert_x)**2)) + vert_y
            left_yspace.append(left_y)
        for x in right_xspace:
            right_y = (a2*((x - vert_x)**2)) + vert_y
            right_yspace.append(right_y)
        overall_xspace = left_xspace + right_xspace
        overall_yspace = left_yspace + right_yspace
        self.barriers.append((overall_xspace, overall_yspace, ls, color, vert_x, vert_y))
    

    def plot(self, ax: plt.Axes = None):
        '''
        This is a method of the FED class that will plot the energy diagram.

        Parameters
        ----------
        ylabel : str
            The label for the y axis. "Energy (eV)"
        ax : plt.Axes
            The axes to plto to. If not specified, a figure and axis will
            be created.   

        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot())     
        '''
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect = self.aspect)
        else:
            self.ax = ax
            self.fig = ax.figure
            self.ax.set_aspect(self.aspect)
        
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        level_xpos_list = []

        self.auto_adjust()

        data = list(zip(self.energies,
                        self.positions,
                        self.bottom_texts,
                        self.top_texts,
                        self.colours,
                        self.right_texts,
                        self.left_texts,
                        self.legends))
        
        for level in data:
            start = level[1]*(self.dimension+self.space)
            ax.hlines(level[0], start, start+self.dimension, color=level[4], label=level[7], linewidth=self.level_width)
            level_centre = start + 0.5*self.dimension
            level_xpos_list.append(level_centre)
            
            ax.text(start+self.dimension/2,  # X
                    level[0]+self.offset,  # Y
                    level[3],  # self.top_texts
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color=level[4])        

            ax.text(start + self.dimension + self.offset,  # X
                    level[0],  # Y
                    level[5],  # self.right_text
                    horizontalalignment='left',
                    verticalalignment='center',
                    color=level[4])

            ax.text(start - self.offset,  # X
                    level[0],  # Y
                    level[6],  # self.left_text
                    horizontalalignment='right',
                    verticalalignment='center',
                    color=level[4])

            ax.text(start + self.dimension/2,  # X
                    level[0] - self.offset,  # Y
                    level[2],  # self.bottom_text
                    horizontalalignment='center',
                    verticalalignment='top',
                    color=level[4])
    
        xtick_pos = []
        for item in level_xpos_list:
            if item not in xtick_pos:
                xtick_pos.append(item)
        
        ax.set_xticks(xtick_pos) 
        ax.set_xticklabels(self.reaction_coords)
        ax.tick_params(axis='x', colors='white', labelcolor='black')

        for idx, link in enumerate(self.links):
            # here we connect the levels with the links
            # x1, x2   y1, y2
            for i in link:
                # i is a tuple: (end_level_id,ls,linewidth,color)
                start = self.positions[idx]*(self.dimension+self.space)
                x1 = start + self.dimension
                x2 = self.positions[i[0]]*(self.dimension+self.space)
                y1 = self.energies[idx]
                y2 = self.energies[i[0]]
                line = Line2D([x1, x2], [y1, y2],
                            ls=i[1],
                            linewidth=i[2],
                            color=i[3])
                ax.add_line(line)


        
        for i in self.barriers:
            ax.plot(i[0], i[1], ls=i[2],
                    linewidth=self.barrier_width, color=i[3])
        plt.legend(loc='upper left')            




                    

            
                    
            
        

