#!/bin/env python

import shutil
import os
from pathlib import Path

outdir = '../MDP/transformation/solvated/fe_calcs'
mdpfiles = ['../MDP/templates/solvated/min_edge3.mdp', '../MDP/templates/solvated/nvt_edge3.mdp', '../MDP/templates/solvated/npt_edge3.mdp', '../MDP/templates/solvated/production_edge3.mdp']
lambdanr=20

if not os.path.isdir(outdir): os.mkdir(outdir)
for i in range(0, lambdanr):
    # Copy and edit mdp files
    for mdpfile in mdpfiles:
        mdpfilename = os.path.join(outdir, Path(mdpfile).name.replace('edge3.', 'edge3_%s.' % i))
        #breakpoint()
        shutil.copy(mdpfile, mdpfilename)
        file = open(mdpfilename, 'r')
        text = file.readlines()
        file.close()
        file = open(mdpfilename, 'w')
        for idx, line in enumerate(text):
           if 'init_lambda_state' in line and line[0] != ';':
                 line = 'init_lambda_state        = %s\n' % i
           file.write(line)
        file.close()
