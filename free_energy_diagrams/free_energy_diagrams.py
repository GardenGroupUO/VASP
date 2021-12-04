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
from matplotlib.lines import Line2D as line




class FED():
	

    def __init__(self):
        self.pos_number = 0
        self.energies = []
        self.positions = []
        self.bottom_texts = []
        self.top_texts = []
        self.left_texts = []
        self.right_texts = []
        self.colours = []
        self.links = []
        


    
    def add_level(self, energy, bottom_text='', position=None, 
                top_text='Energy', left_text='', right_text='', color='k'):
        '''
        This is a method of the FED class. It will save all of the data
        for a new energy level into the lists required.

        Parameters
        __________
        energy : int
            The energy of the level in eV
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
        _______
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
        
        if top_text == 'Energy':
            top_text = energy
        
        link = []
        
        self.energies.append(energy)
        self.positions.append(position)
        self.bottom_texts.append(bottom_text)
        self.top_texts.append(top_text)
        self.left_texts.append(left_text)
        self.right_texts.append(right_text)
        self.colours.append(color)
        self.links.append(link)

    def add_link(self, start_level_id, end_level_id, color = 'k',
             ls='--', linewidth=1):
        '''
        This is a method of the FED class that will take in a stat and end
        level id (start/end_level_id) and create a linker between those
        two levels. This will be used for levels that do not have a barrier
        between them.

        Parameters
        __________
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
        _______
        Appends all of the imformation about the linker to a list in the
        class link attribute
        '''

        self.links[start_level_id].append((end_level_id, ls, linewidth, color))

