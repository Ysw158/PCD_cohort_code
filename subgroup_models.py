"""Run variant-type models, miscarriage models and genotype–motility models."""
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.stats import false_discovery_control
from firth import fit
from firth_motility import fit as fit_motility

people = json.loads(Path(sys.argv[1]).read_text())
output = Path(sys.argv[2])
output.mkdir(parents=True, exist_ok=True)
config = list(csv.DictReader(Path(__file__).with_name('subgroup_pairs.csv').open()))
fields = {'S4': 'zygosity', 'S5': 'variant', 'S6': 'cbp'}
def key(value):
    return ''.join(value.lower().split())

def estimate(data, group, outcome, intervals=True, estimator=fit):
    flag = np.array([group(r) for r in data], dtype=float)
    age = np.array([r['enrollment_age'] for r in data])
    sex = np.array([r['sex']=='M' for r in data])
    y = np.array([outcome(r) for r in data])
    counts = [sum((flag==g) & (y==v)) for g,v in [(1,1),(1,0),(0,1),(0,0)]]
    result = dict(zip(['n','group_yes','group_no','other_yes','other_no'], [len(y),*counts]))
    X = np.column_stack([np.ones(len(y)),flag])
    if len(set(flag))<2 or len(set(y))<2:
        return result
    for suffix in ['raw','adj']:
        if suffix=='adj':
            if np.ptp(age): X=np.column_stack([X,(age-age.mean())/age.std()])
            if np.ptp(sex.astype(int)): X=np.column_stack([X,sex])
        z=estimator(X,y,ci=intervals)
        result.update({f'{k}_{suffix}':z[k] for k in ['or','lo','hi','p'] if k in z})
    return result

def write(name,rows):
    with (output/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader();writer.writerows(rows)

rinput=[]
for r in people:
    rinput.append(dict(age=r['enrollment_age'],sex=r['sex'], **{k:key(r[v]) for k,v in fields.items()}, **r['phenotypes']))
write('private_subgroups.csv',rinput)  # controlled input to R; keep output directory private
rows=[]
for c in config:
    family,p=c['family'],c['phenotype']
    if family=='S4' and p!='Maternal miscarriage':continue
    if family=='S6' and p=='Infertility':continue
    field=fields[family]
    data=[r for r in people if r[field] and r['phenotypes'][p] is not None and (family!='S5' or r['genes']) and r['enrollment_age'] is not None and r['sex'] in ('M','F')]
    rows.append({**c,**estimate(data,lambda r:key(r[field])==c['key'],lambda r:r['phenotypes'][p],estimator=fit_motility if family=='S6' and p!='Maternal miscarriage' else fit)})
for family,p in {(r['family'],r['phenotype']) for r in rows}:
    for suffix in ['raw','adj']:
        subset=[r for r in rows if (r['family'],r['phenotype'])==(family,p) and f'p_{suffix}' in r]
        if subset:
            for r,q in zip(subset,false_discovery_control([r[f'p_{suffix}'] for r in subset])):
                r[f'q_{suffix}']=q
write('python_subgroups.csv',rows)
data=[r for r in people if r['genes'] and r['cbp']]
assert len(data)==66
groups=['+'.join(sorted(r['genes'])) for r in data]
cells=[]
for g in sorted(set(groups)):
    if groups.count(g)<3:continue
    for pattern in sorted({r['cbp'] for r in data}):
        cells.append(dict(gene=g,pattern=pattern,**estimate(data,lambda r: r['genes']==[g],lambda r:r['cbp']==pattern,False,fit_motility)))
write('figure6f.csv',cells)
