#!/usr/bin/python3
from collections import defaultdict
import re

with open('China_full_snp.xls', 'r') as f1, open('panel.list', 'r') as f2, open('panel.xls','w') as O:
    O.write('Drug\tFeature\tOrigin mutations\tConfidence\tImportance\n')
    dct = defaultdict(lambda: defaultdict(str))
    for l in f1.readlines()[1:]:
        u = l.strip().split('\t')
        drug, org, vartype, ci, var = u[0], u[2], u[5], u[6], u[-1]

        if re.search('stop_lost|start_lost|stop_gained|transcript_ablation|feature_ablation|frameshift_variant', vartype):
            org = var
        dct[drug.capitalize()][var] = f'{drug}\t{var}\t{org}\t{ci}'

    for l in f2.readlines()[1:]:
        Drug, Var, Imp = l.strip().split('\t')
        O.write(f'{dct[Drug][Var]}\t{Imp}\n')
