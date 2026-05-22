bed = "/project/Research_Service/software/miniconda3/envs/tbprofiler-5.0.1/share/tbprofiler/tbdb.bed"

import sys,os
snps = sys.argv[1]

Var = {}
count = {}
drug = {}
from collections import defaultdict
other = defaultdict(lambda: defaultdict(str))
drugR = defaultdict(lambda: defaultdict(int))

import re

colors = {
	"rifampicin":"#2BAE9E",
	"isoniazid":"#65A2D2",
	"ethambutol":"#D4D98A",
	"pyrazinamide":"#72B285",
	"streptomycin":"#BF95C1",
	"fluoroquinolones":"#F58E87",
	"amikacin":"",
	"capreomycin":"#EB8B39",
	"kanamycin":"",
	"cycloserine":"",
	"ethionamide":"#E36146",
	"clofazimine":"#65A2D2",
	"para-aminosalicylic_acid":"#D4D98A",
	"delamanid":"",
	"bedaquiline":"#2BAE9E",
	"linezolid":""
}

colors2 = {
    "rifampicin":"#2BAE9E",
    "isoniazid":"#0074B3",
    "ethambutol":"#892ADE",
    "pyrazinamide":"#F39F4E",
    "streptomycin":"#FF6347",
    "fluoroquinolones":"#008A99",
    "amikacin":"#FBA400",
    "capreomycin":"#A5D395",
    "kanamycin":"#F58E87",
    "cycloserine":"#61C1BF",
    "ethionamide":"#BF95C1",
    "clofazimine":"#FEDD97",
    "para-aminosalicylic_acid":"#72B285",
    "delamanid":"#E1807E",
    "bedaquiline":"#4D4D9A",
    "linezolid":"#33ABC1"
}


cols = [
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#B17BA6",
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C",
    "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02",
    "#7570B3", "#BEAED4", "#666666", "#999999", "#aa8282",
    "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", "#aeae5c",
    "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#882E72"
]

def hex_to_rgb(hex_color, alpha=None):
	"""
	将十六进制颜色代码转换为RGB元组
	
	参数:
		hex_color: 十六进制颜色字符串，如 "#FF0000" 或 "FF0000"
	
	返回:
		(r, g, b) 元组，每个值在0-255之间
	"""
	# 移除可能存在的#号
	hex_color = hex_color.lstrip('#')
	
	# 检查长度
	if len(hex_color) == 3:  # 简写格式，如 "F00"
		hex_color = ''.join([c*2 for c in hex_color])  # 扩展为完整格式
	
	# 转换为RGB
	r = int(hex_color[0:2], 16)
	g = int(hex_color[2:4], 16)
	b = int(hex_color[4:6], 16)
	
	if not alpha:
		return f"({r},{g},{b})"
	else:
		return f"({r},{g},{b},{alpha})"


# 使用示例
#    rgb = hex_to_rgb(color)


def split_interval_simple(start, end, n):
	total = end - start + 1
	base = total // n
	rem = total % n
	
	result = []
	s = start
	for i in range(n):
		length = base + (1 if i < rem else 0)
		e = s + length - 1
		result.append((s, e))
		s = e + 1
	
	return result


#from collections import OrderedDict
dedup = lambda lst: list(dict.fromkeys(lst))
gene_order = []

with open(snps, 'r') as f:
	for l in f.readlines()[1:]:
		Drug,Variant,Pvar,Gvar,Pos,VarType,Confidence,Freq,Total,Ratio,Date_World,Date_country,_ = l.rstrip('\n').split('\t')
		if Drug in ['moxifloxacin', 'ofloxacin', 'levofloxacin', 'ciprofloxacin']:
			continue
		gene, var = Pvar.split('_', maxsplit=1)
		gene_order.append(gene)
		Var.setdefault(gene, []).append(var)
#		drug.setdefault(gene, set()).add(Drug)
		drug[gene] = Drug
		other[Pvar]['ci'] = Confidence
		fRatio = float(re.sub('%', '', Ratio))/100
		other[Pvar]['Ratio'] = f'{fRatio:.2f}'
		other[Pvar]['DateWorld'] = Date_World
		other[Pvar]['DateCountry'] = Date_country

gene_order = dedup(gene_order)


with open("prepare_data/karyotype.txt", 'w') as O1, open('prepare_data/labels.txt', 'w') as O2, open('prepare_data/gene_label.txt', 'w') as O3, open('prepare_data/drug_colors.txt', 'w') as O4, open('prepare_data/drug_label.txt', 'w') as O5, open('prepare_data/confidence.txt', 'w') as O6, open('prepare_data/freq.txt', 'w') as O7, open('prepare_data/gene_color.txt', 'w') as O8, open('prepare_data/date_world.txt', 'w') as O9, open('prepare_data/date_country.txt', 'w') as O10:
#chr - tb_chr 1 0 4411532 lgrey
#	for g,vs in Var.items():
	grS = 1
	grE = 0
	drug_order = []
	it = 0
	record = []
	end = 0
	for g in gene_order:
		it += 1
		vs = Var[g]
		Hex = colors2[drug[g]]
		rgb = hex_to_rgb(Hex)
#		gc = hex_to_rgb(cols[it], 0.7)
		if drug[g] in record:
			gc = hex_to_rgb(colors2[drug[g]], 0.3)
		else:
			gc = hex_to_rgb(colors2[drug[g]], 0.6)
			record.append(drug[g])

		if len(vs) < 3:
			grE = grS + 500 - 1

			O3.write(f'H37Rv {grS} {grE} {g}\n')
			O8.write(f'H37Rv {grS} {grE} fill_color={gc}\n')
		else:
			grE = grS + 200*len(vs) - 1

			O3.write(f'H37Rv {grS} {grE} {g}\n')
			O8.write(f'H37Rv {grS} {grE} fill_color={gc}\n')

		drug_order.append(drug[g])

		if drug[g] in drugR:
			drugR[drug[g]]['E'] = grE
		else:
			drugR[drug[g]]['S'] = grS
			drugR[drug[g]]['E'] = grE

		c = 0
		for i in split_interval_simple(1,1000,len(drug[g])):
		#	Hex = colors[list(drug[g])[c]]
		#	rgb = hex_to_rgb(Hex)
		#	O4.write(f'{g} {i[0]} {i[1]} fill_color={rgb}\n')
		#	O5.write(f'{g} {i[0]} {i[1]} {list(drug[g])[c]}\n')
		#	c+=1
			pass	
		

		C = 0
		for i in split_interval_simple(grS, grE, len(vs)):
			if len(vs[C]) > 30:
				O2.write(f'H37Rv {i[0]} {i[1]} {re.sub("del.*","del",re.sub("_","-",vs[C]))} \n')

			elif len(vs[C]) > 22:
				O2.write(f'H37Rv {i[0]} {i[1]} {re.sub("_","-",vs[C])} \n')

			else:
				O2.write(f'H37Rv {i[0]} {i[1]} {re.sub("_","-",vs[C])} \n')

			fv = g + '_' + vs[C]
			fc = other[fv]["ci"]
			if fc == "":
				fc = "tbdb"
			if fc.startswith('Assoc w R'):
				O6.write(f'H37Rv {i[0]} {i[1]} 1 glyph=triangle,color=vvdred\n')
			elif fc == "Uncertain significance":
				O6.write(f'H37Rv {i[0]} {i[1]} 1 glyph=triangle,color=lgreen\n')
			else:
				O6.write(f'H37Rv {i[0]} {i[1]} 1 glyph=triangle,color=vlgreen\n')
				
			freqsize = 3+float(other[fv]["Ratio"])*11
			freqsize = 9 if freqsize > 9 else freqsize

			O7.write(f'H37Rv {i[0]} {i[1]} {other[fv]["Ratio"]} glyph_size={freqsize}\n')
			O9.write(f'H37Rv {i[0]} {i[1]} {other[fv]["DateWorld"]}\n')
			O10.write(f'H37Rv {i[0]} {i[1]} {other[fv]["DateCountry"]}\n')
			end = i[1]
			C += 1

		grS = grE + 1 + 100
	O1.write(f'chr - H37Rv H37Rv 0 {grE+100} lgrey\n')
	O10.write(f'H37Rv {end} {end+1} 2015\n')

	for d in dedup(drug_order):
		rgb = hex_to_rgb(colors2[d])
		O4.write(f'H37Rv {drugR[d]["S"]} {drugR[d]["E"]} fill_color={rgb}\n')
#		if d == "para-aminosalicylic_acid":
#			O5.write(f'H37Rv {drugR[d]["S"]} {drugR[d]["E"]} PAS\n')
#		else:
		O5.write(f'H37Rv {drugR[d]["S"]} {drugR[d]["E"]} {d}\n')
