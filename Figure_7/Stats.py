country = ["Australia", "China", "India", "Moldova", "USA"]
firstlineDrug = ["isoniazid", "rifampicin", "pyrazinamide", "ethambutol", "streptomycin"]

O = open("Sup3Stats.xls", 'w')
O.write('''#'moxifloxacin', 'ofloxacin', 'levofloxacin', 'ciprofloxacin'合并为fluoroquinolones统计
#Total: 去重后所有突变数量
#Local: 去重后date_world==date_country的突变数量
#tFirstline: Local内属于一线药物的突变数量
#Firstline(>1): 突变数目大于1的一线药物的突变数量
#Secondline(>1): 突变数目大于1的一线药物的突变数量 
Firstline_ratio(>1) : Firstline(>1)/( Firstline(>1) + Secondline(>1) ) * 100
''')
O.write(f'Country\tTotal\tLocal\tFirstline\tFirstline(>1)\tSecondline\tSecondline(>1)\tFirstline_ratio(>1)\ttSecondline_ratio(>1)\n')
import pprint
for c in country:
    local = set()
    Nonlocal = set()
    firstline = set()
    secondline = set()
    drug = {}

    with open(f"{c}_decade_snp.xls", "r") as f:
        for l in f.readlines()[1:]:
        #     #Drug    Variant Pro_change      GenomeVar       Pos     VarType Confidence      Freq    Total   Ratio   Date_World      Date_Australia  format_Var
            u = l.rstrip('\n').split('\t')            
            if u[0] in ['moxifloxacin', 'ofloxacin', 'levofloxacin', 'ciprofloxacin']:
                continue

            if u[10] == u[11]:
                local.add(u[2])
                if u[0] in firstlineDrug:
                    firstline.add(u[2])
                    drug.setdefault('first', {}).setdefault(u[0], set()).add(u[2])
                else:
                    secondline.add(u[2])
                    drug.setdefault('second', {}).setdefault(u[0], set()).add(u[2])
            else:
                Nonlocal.add(u[2])




        firstline_gt1 = 0
        secondline_gt1 = 0
        for i in drug['first'].keys():
            if len(drug['first'][i]) > 1:
                firstline_gt1 += len(drug['first'][i])
        
        for i in drug['second'].keys():
            if len(drug['second'][i]) > 1:
                secondline_gt1 += len(drug['second'][i])          
        
        #O.write(f'{c}\t{len(local)+len(Nonlocal)}\t{len(local)}\t{len(firstline)}\t{round(len(firstline)/len(local), 2)}\n')


        O.write(f'{c}\t{len(local)+len(Nonlocal)}\t{len(local)}\t{len(firstline)}\t{firstline_gt1}\t{len(secondline)}\t{secondline_gt1}\t{round(firstline_gt1/(secondline_gt1+firstline_gt1)*100, 2)}\t{round(secondline_gt1/(secondline_gt1+firstline_gt1)*100, 2)}\n')  #(secondline_gt1+firstline_gt1)

            
                