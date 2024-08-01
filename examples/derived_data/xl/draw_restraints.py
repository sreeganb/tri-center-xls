#!/usr/bin/env python3.11

import csv
import numpy as np

from pymol import cmd, CmdException

from pymol.cgo import *

colors = {

        'green': [0.0, 1.0, 0.5], # limegreen

        'red': [0.82, 0.0, 0.3] # dubnium

        }

radius = 0.5

def drawRestraints(fname, selection='all', prefix='xl', threshold=30.0, atom='CA', quiet=1):
    with open(fname, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        print(f'Loading restraints from {reader}')
        for row in reader:
            try:
                x1, y1, z1 = cmd.get_coords(f'{selection} and chain {row["Chain ID A"]} and resi {row["XL A"]} and name {atom}', 1)[0]
                x2, y2, z2 = cmd.get_coords(f'{selection} and chain {row["Chain ID B"]} and resi {row["XL B"]} and name {atom}', 1)[0]
            except:
                continue
            d = np.linalg.norm(np.array([x2, y2, z2]) - np.array([x1, y1, z1]))
            print('Distance:', d)
            if d <= float(threshold):
                r1, g1, b1  = colors['green']
                r2, g2, b2  = colors['green']
            else:
                r1, g1, b1  = colors['red']
                r2, g2, b2  = colors['red']
            cmd.load_cgo([CYLINDER, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2],
                f'{prefix}_{row["XL A"]}_{row["XL B"]}_{atom}')
    cmd.group(prefix, f'{prefix}_*')

cmd.extend('drawRestraints', drawRestraints)