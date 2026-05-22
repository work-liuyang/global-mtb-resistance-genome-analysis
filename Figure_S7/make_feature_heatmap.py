#!/usr/bin/python3

drug_order = ["Isoniazid", "Rifampicin", "Ethambutol", "Ethionamide", "Amikacin", "Levofloxacin", "Moxifloxacin", "Kanamycin"]

from collections import defaultdict
import re


def get_features_double(panel):
    IMP = defaultdict(lambda: defaultdict(str))
    CI = defaultdict(str)

    with open(panel, 'r') as f:
        for l in f.readlines()[1:]:
            drug, _, ovar, ci, Imp = l.strip().split('\t')
            IMP[ovar][drug.capitalize()]= Imp
            if re.match('Assoc w R', ci):
                ci = "AssocwR"
            elif re.match('Uncertain', ci):
                ci = "Uncertain"
            else:
                pass           
            CI[ovar] = ci

    return (IMP, CI)


def prepare_heatmap(IMP, CI, abbr):
    with open('plot_heatmap.xls', 'w') as O:
        O.write(f'features\t{"\t".join(abbr)}\tConfidence\n')
        for k in IMP.keys():
            line = k
            for d in abbr:
                line += f'\t{IMP[k][d]}'
                
            O.write(f'{line}\t{CI[k]}\n')






IMP, CI = get_features_double('panel.xls')
import pprint
#pprint.pprint(MyDict)
#print(MyDict['katG_p.Ser315Thr']['Isoniazid']['ci'])

prepare_heatmap(IMP, CI, drug_order)
