import os,sys,glob

files = glob.glob('roh1kb.*.txt')

match = []
for f in files:
	proband = f.split('.')[2]
	if proband not in match:
		match.append(proband)

probands = glob.glob('roh1kb*p1*txt')

match_p1 = []
for f in probands:
	proband = f.split('.')[2]
	if proband not in match_p1:
		match_p1.append(proband)

siblings = glob.glob('roh1kb*s1*txt')
match_s1 = []
for f in siblings:
	sibling = f.split('.')[2]
	if sibling not in match_s1:
		match_s1.append(sibling)

for m in match:
	if m not in match_p1:
		print 'missing proband in ' + m
	else:
		pass
	if m not in match_s1:
		print 'missing sibling in ' + m
	else:
		pass
