#!/bin/env python

import shutil
import os
from pathlib import Path

outdir = '../MDP/transformation/solvated/fe_calcs'
mdpfiles = ['../MDP/templates/solvated/min_edge1.mdp', '../MDP/templates/solvated/nvt_edge1.mdp', '../MDP/templates/solvated/npt_edge1.mdp', '../MDP/templates/solvated/production_edge1.mdp']
lambdanr=16

if not os.path.isdir(outdir): os.mkdir(outdir)
for i in range(0, lambdanr):
    # Copy and edit mdp files
    for mdpfile in mdpfiles:
        mdpfilename = os.path.join(outdir, Path(mdpfile).name.replace('edge1.', 'edge1_%s.' % i))
        #breakpoint()
        shutil.copy(mdpfile, mdpfilename)
        file = open(mdpfilename, 'r')
        text = file.readlines()
        file.close()
        file = open(mdpfilename, 'w')
        for idx, line in enumerate(text):
           if 'define' in line and line[0] != ';':
                 line = line.strip() + ' -DEXCLS_ON -DEDGE1_STEP%s\n' % i
           file.write(line)
        file.close()
