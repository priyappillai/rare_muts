import pandas as pd

chrms = [str(i) for i in range(1, 23)]
chrms.extend(['X', 'Y'])

out = pd.DataFrame()
for i in chrms:
	temp = pd.read_csv('fit_params_'+i+'.txt')
	temp['Chromosome'] = i
	out = pd.concat([out, temp])

out.to_csv('fit_params.txt')