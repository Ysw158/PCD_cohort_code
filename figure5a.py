import json,collections,csv
import numpy as np
from pathlib import Path
from firth import fit,design
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=Path, help='JSON with one record per participant')
parser.add_argument('output', type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
d = {'unique': json.loads(args.input.read_text())}
with Path(__file__).with_name('model_pairs.csv').open() as stream:
    d['old'] = list(csv.DictReader(stream))
assert len(d['unique']) == len({r['patient_id'] for r in d['unique']}) == 130
assert all(v in (0, 1, None) for r in d['unique'] for v in r['phenotypes'].values())
data=[r for r in d['unique'] if r['genes']]
out=[]
def counts(rows,g,p):return [sum((g in r['genes'])==a and r['phenotypes'][p]==b for r in rows) for a,b in [(True,1),(True,0),(False,1),(False,0)]]
for i,old in enumerate(d['old']):
    g,p=old['gene'],old['phenotype'];X,y,rows=design(data,g,p)
    rec={'gene':g,'phenotype':p,'n':len(rows),'counts':counts(rows,g,p),'covariates':'Age and sex' if X.shape[1]==4 else 'Age (sex constant)', 'model_pair':old}
    if len(set(y))<2 or len(set(X[:,1]))<2:
        rec['status']='Not estimable: no outcome or gene variation'
    else:
        try:
            rec['adj']=fit(X,y)
            Xu,yu,ru=design(data,g,p,adjust=False)
            rec['unadj']=fit(Xu,yu);rec['unadj_n']=len(ru);rec['unadj_counts']=counts(ru,g,p)
            rec['status']='Estimated'
        except Exception as e:
            rec['status']='Fit failure: '+str(e)
    rec['carrier_n']=sum(rec['counts'][:2]);rec['comparator_n']=sum(rec['counts'][2:]);rec['display_eligible']=rec['carrier_n']>=3 and rec['comparator_n']>=3
    out.append(rec)
    if (i+1)%30==0:print('Completed',i+1,flush=True)
def bh(records,key,result):
    good=[r for r in records if key in r]
    order=sorted(good,key=lambda r:r[key]['p']);m=len(order);last=1.
    for i in range(m-1,-1,-1):
        last=min(last,order[i][key]['p']*m/(i+1));order[i][key][result]=last
for p in sorted(set(r['phenotype'] for r in out)):
    sub=[r for r in out if r['phenotype']==p]
    bh(sub,'adj','q_within_phenotype');bh(sub,'unadj','q_within_phenotype')
bh(out,'adj','q_all_tests');bh(out,'unadj','q_all_tests')
lookup = {(r['gene'], r['phenotype']): r for r in out}
with Path(__file__).with_name('display_pairs.csv').open() as stream:
    selected = [lookup[r['gene'], r['phenotype']] for r in csv.DictReader(stream)]
assert len(data) == 125 and len(out) == 431
print([(r['gene'],r['phenotype'],r['status']) for r in out if r['status']!='Estimated'],flush=True)
assert all(not r['status'].startswith('Fit failure') for r in out)
assert all('adj' in r for r in selected)
result={'cohort_n':130,'genetic_n':len(data),'model_count':len([r for r in out if 'adj' in r]),'selected':selected,'results':out}
(args.output/'figure5a_results.json').write_text(json.dumps(result,indent=2))
print('MODELS',result['model_count'],'SELECTED',len(selected),'FAILURES',[r['status'] for r in out if r['status']!='Estimated'])
for r in selected:print(r['gene'],r['phenotype'],r['counts'],r['adj'])
# Independent saturated 2x2 check of the Firth point estimator.
for r in out:
    if 'unadj' not in r:continue
    a,b,c,e=r['unadj_counts'];exact=(a+.5)*(e+.5)/((b+.5)*(c+.5))
    assert np.isclose(exact,r['unadj']['or'],rtol=1e-5),(r['gene'],r['phenotype'],exact,r['unadj'])
print('All unadjusted point estimates pass independent 2x2 analytical check.')

with (args.output/'figure5a_display.csv').open('w',newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['gene','phenotype','carrier_yes','carrier_no','other_yes','other_no','n','OR','lower95','upper95','p'])
    for row in selected:
        est = row['adj']
        writer.writerow([row['gene'],row['phenotype'],*row['counts'],row['n'],est['or'],est['lo'],est['hi'],est['p']])
