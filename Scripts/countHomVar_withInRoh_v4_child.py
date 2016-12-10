import os,sys

inputfile = sys.argv[1] #'20160818.roh1kb.gVCF.merged.txt'

sampleID = inputfile.split('.')[2] + '.' + inputfile.split('.')[3]
familyID = sampleID.split('.')[0]
role = sampleID.split('.')[1]
outputfile = 'roh1kb.CountV4.' + sampleID + '.txt'
print outputfile

roles_idx = [-4,-3,-2,-1]
roles = ['p1','s1','mo','fa']
role_idx = roles_idx[roles.index(role)]

o = open(outputfile,'w')
header = ['#Chrom','Start','End','Size','Avg_Score',
			'Hom_PRO','Het_PRO','Ref_PRO',
			'Hom_SIB','Het_SIB','Ref_SIB',
			'Hom_MAT','Het_MAT','Ref_MAT',
			'Hom_PAT','Het_PAT','Ref_PAT',
			'HetHom_PRO','HetHom_SIB','HetHom_MAT','HetHom_PAT','QUAL']
o.write('\t'.join(header)+'\n')

funMatFile = 'funMat.' + familyID + '.txt'
tempOutFile = 'tempOut.' + familyID + '.' + sampleID + '.txt'
cmd = 'intersectBed -a ' + inputfile + ' -b ' + funMatFile + ' -wa -wb > ' + tempOutFile
print cmd
os.system(cmd)

pre_start = '0'
start = ''
end = ''
out = []

'''
1       4037665 4040926 93.9897435897   1       4039896 4039896 C       T       0       0       Intergenic      2222    0.443291
1       4037665 4040926 93.9897435897   1       4040155 4040155 C       T       0       0       Intergenic      2222    0.712
1       4049778 4053286 91.0614285714   1       4049980 4049980 C       T       0       0       Intergenic      2222    0.442292
1       4049778 4053286 91.0614285714   1       4050953 4050953 A       G       0       0       Intergenic      2222    0.704
1       4049778 4053286 91.0614285714   1       4051249 4051249 T       C       0       0       Intergenic      2222    0.445288
1       4053321 4057886 91.5019607843   1       4054067 4054067 C       T       0       0       Intergenic      2222    0.508586
'''

loop = 0

with open(tempOutFile) as fh:
	for line in fh:
		info = line.rstrip('\n').split('\t')
		chrom,start,end = info[0],info[1],info[2]
		avg_score = info[3]
		size = str(int(end) - int(start))
		GT = info[12] #2011
		print GT
		if start != pre_start: # start new ROH block
			if loop == 0: # skip the first line
				pass
			else:
				# Quality check
				tag = ''
				HetHom_PRO, HetHom_SIB, HetHom_MAT, HetHom_PAT = 0, 0 ,0 ,0
				if Het_PRO == 0 and Hom_PRO == 0:
					pass
				else:
					HetHom_PRO = Hom_PRO/float(Het_PRO+Hom_PRO)
				if Het_SIB == 0 and Hom_SIB == 0:
					pass
				else:
					HetHom_SIB = Hom_SIB/float(Het_SIB+Hom_SIB)
				if Het_MAT == 0 and Hom_MAT == 0:
					pass
				else:
					HetHom_MAT = Hom_MAT/float(Het_MAT+Hom_MAT)
				if Het_PAT == 0 and Hom_PAT == 0:
					pass
				else:
					HetHom_PAT = Hom_PAT/float(Het_PAT+Hom_PAT)
				print HetHom_PRO, HetHom_SIB, HetHom_MAT, HetHom_PAT
				if role == 'p1':
					if HetHom_PRO >= 0.95 and HetHom_MAT < 0.05 and HetHom_PAT < 0.05:
						tag = 'P'
					else:
						tag = 'F'
				else:
					if HetHom_SIB >= 0.95 and HetHom_MAT < 0.05 and HetHom_PAT < 0.05:
						tag = 'P'
					else:
						tag = 'F'
				out = out + [str(Hom_PRO),str(Het_PRO),str(Ref_PRO),str(Hom_SIB),str(Het_SIB),str(Ref_SIB),
							str(Hom_MAT),str(Het_MAT),str(Ref_MAT),str(Hom_PAT),str(Het_PAT),str(Ref_PAT),
							str(HetHom_PRO),str(HetHom_SIB),str(HetHom_MAT),str(HetHom_PAT),tag]
				line_out = '\t'.join(out) + '\n'
				o.write(line_out)
			# start new block
			out = [chrom,start,end,size,avg_score]
			Hom_PRO,Het_PRO,Ref_PRO=0,0,0
			Hom_SIB,Het_SIB,Ref_SIB=0,0,0
			Hom_MAT,Het_MAT,Ref_MAT=0,0,0
			Hom_PAT,Het_PAT,Ref_PAT=0,0,0
			# Proband
			if GT[0] == '0':
				Ref_PRO = Ref_PRO + 1
			elif GT[0] == '1':
				Het_PRO = Het_PRO + 1
			elif GT[0] == '2':
				Hom_PRO = Hom_PRO + 1
			else:
				pass
			# Sibling
			if GT[1] == '0':
				Ref_SIB = Ref_SIB + 1
			elif GT[1] == '1':
				Het_SIB = Het_SIB + 1
			elif GT[1] == '2':
				Hom_SIB = Hom_SIB + 1
			else:
				pass
			# Mother
			if GT[2] == '0':
				Ref_MAT = Ref_MAT + 1
			elif GT[2] == '1':
				Het_MAT = Het_MAT + 1
			elif GT[2] == '2':
				Hom_MAT = Hom_MAT + 1
			else:
				pass
			# Father
			if GT[3] == '0':
				Ref_PAT = Ref_PAT + 1
			elif GT[3] == '1':
				Het_PAT = Het_PAT + 1
			elif GT[3] == '2':
				Hom_PAT = Hom_PAT + 1
			else:
				pass
		else: # continute ROH block
			# Proband
			if GT[0] == '0':
				Ref_PRO = Ref_PRO + 1
			elif GT[0] == '1':
				Het_PRO = Het_PRO + 1
			elif GT[0] == '2':
				Hom_PRO = Hom_PRO + 1
			else:
				pass
			# Sibling
			if GT[1] == '0':
				Ref_SIB = Ref_SIB + 1
			elif GT[1] == '1':
				Het_SIB = Het_SIB + 1
			elif GT[1] == '2':
				Hom_SIB = Hom_SIB + 1
			else:
				pass
			# Mother
			if GT[2] == '0':
				Ref_MAT = Ref_MAT + 1
			elif GT[2] == '1':
				Het_MAT = Het_MAT + 1
			elif GT[2] == '2':
				Hom_MAT = Hom_MAT + 1
			else:
				pass
			# Father
			if GT[3] == '0':
				Ref_PAT = Ref_PAT + 1
			elif GT[3] == '1':
				Het_PAT = Het_PAT + 1
			elif GT[3] == '2':
				Hom_PAT = Hom_PAT + 1
			else:
				pass
		pre_out = out
		pre_start = start
		loop = loop + 1

o.close()

rm = 'rm ' + tempOutFile
os.system(rm)
