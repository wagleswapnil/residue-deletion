#!/bin/env python

import shutil
import os
from pathlib import Path

outdir = '../MDP/wt/solvated/fe_calcs'
mdpfiles = ['../MDP/templates/solvated/min_lambda_X.mdp', '../MDP/templates/solvated/nvt_lambda_X.mdp', '../MDP/templates/solvated/npt_lambda_X.mdp', '../MDP/templates/solvated/production_lambda_X.mdp']
lambdanr=20

if not os.path.isdir(outdir): os.makedirs(outdir, exist_ok=True)
for i in range(0, lambdanr):
    # Copy and edit mdp files
    for mdpfile in mdpfiles:
        mdpfilename = os.path.join(outdir, Path(mdpfile).name.replace('_X.', '_%s.' % i))
        #mdpfilename = os.path.join(outdir, Path(mdpfile).name.replace('lambda_X.', 'edge3_%s.' % i))
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
