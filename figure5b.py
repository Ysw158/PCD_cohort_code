import json,itertools,csv
from pathlib import Path
import numpy as np
from scipy.stats import t,f,shapiro
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
d = json.loads(args.input.read_text())
assert len({r['patient_id'] for r in d}) == len(d)
assert all(isinstance(r['multigene_patient'], bool) and r['sex'] in ('M','F') for r in d)
assert sum(r['multigene_patient'] for r in d) == 1
genes=['DNAH5','ODAD1','DNAH11','DNAAF1','RSPH4A','HYDIN','CCDC40'];a=[r for r in d if r['gene'] in genes];models={}
assert len(a)==39 and not any(r['multigene_patient'] for r in a)
assert [sum(r['gene']==g for r in a) for g in genes] == [14,6,4,4,4,4,3]
def design(data,sex=False,interaction=False):
 X=np.array([[1,*[int(r['gene']==g) for g in genes[1:]],r['age_at_spirometry_years'],*([int(r['sex']=='M')] if sex else []),*([int(r['gene']==g)*r['age_at_spirometry_years'] for g in genes[1:]] if interaction else [])] for r in data],dtype=float);return X

def fit(data,label,sex=False):
 X=design(data,sex);y=np.array([r['FEV1_percent_predicted'] for r in data]);inv=np.linalg.inv(X.T@X);b=np.linalg.lstsq(X,y,rcond=None)[0];e=y-X@b;df=len(y)-X.shape[1];rss=e@e;cov=inv*rss/df
 Z=np.array([[1,r['age_at_spirometry_years'],*([int(r['sex']=='M')] if sex else [])] for r in data]);ee=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0];F=((ee@ee-rss)/6)/(rss/df);p=f.sf(F,6,df)
 ageidx=7;se=np.sqrt(cov[ageidx,ageidx]);crit=t.ppf(.975,df)
 models[label]={'n':len(data),'global_F':float(F),'global_p':float(p),'df_num':6,'df_resid':df,'age_effect':float(b[ageidx]),'age_p':float(2*t.sf(abs(b[ageidx]/se),df)),'age_ci':[float(b[ageidx]-crit*se),float(b[ageidx]+crit*se)]}
 print(label,models[label]);return X,y,b,cov,df,e,inv
X,y,b,cov,df,e,inv=fit(a,'primary_single_gene');crit=t.ppf(.975,df);meanage=float(np.mean([r['age_at_spirometry_years'] for r in a]));rows=[]
def vec(g):return np.array([1,*[int(g==x) for x in genes[1:]],meanage])
for gene in genes:
 z=vec(gene);estimate=z@b;se=np.sqrt(z@cov@z);group=[r for r in a if r['gene']==gene];ages=[r['age_at_spirometry_years'] for r in group];ys=[r['FEV1_percent_predicted'] for r in group]
 rows.append({'gene':gene,'n':len(group),'age_mean':float(np.mean(ages)),'age_min':min(ages),'age_max':max(ages),'FEV1_mean_unadjusted':float(np.mean(ys)),'FEV1_SD_unadjusted':float(np.std(ys,ddof=1)),'adjusted_mean_at_common_age':float(estimate),'CI_low':float(estimate-crit*se),'CI_high':float(estimate+crit*se),'standardization_age':meanage})
pairrows=[]
for g1,g2 in itertools.combinations(genes,2):
 z=vec(g1)-vec(g2);estimate=z@b;se=np.sqrt(z@cov@z);p=2*t.sf(abs(estimate/se),df)
 pairrows.append({'gene_1':g1,'gene_2':g2,'adjusted_difference_percentage_points':float(estimate),'CI_low':float(estimate-crit*se),'CI_high':float(estimate+crit*se),'p_unadjusted':float(p)})
order=np.argsort([r['p_unadjusted'] for r in pairrows]);prev=0
for rank,idx in enumerate(order):
 prev=max(prev,min(1,(len(order)-rank)*pairrows[idx]['p_unadjusted']));pairrows[idx]['p_Holm']=prev
for assign in ['DNAH5','DNAH11']:
 dd=[dict(r,gene=assign) if r['multigene_patient'] else r.copy() for r in d if r['gene'] in genes or r['multigene_patient']];fit(dd,'dual_gene_assigned_'+assign)
fit(a,'age_sex_sensitivity',True)
h=np.sum((X@inv)*X,axis=1);u=e/(1-h);hc3=inv@(X.T@(u[:,None]**2*X))@inv;sub=b[1:7];wald=float(sub@np.linalg.inv(hc3[1:7,1:7])@sub/6);models['HC3_sensitivity']={'global_F':wald,'global_p':float(f.sf(wald,6,df))}
print('HC3',models['HC3_sensitivity'],'residual Shapiro',shapiro(e),'maxleverage',h.max());print('groups',rows);print('nominal',[r for r in pairrows if r['p_unadjusted']<.05]);print('holm',[r for r in pairrows if r['p_Holm']<.05])
for name,data in [('FEV1_age_adjusted_group_results.csv',rows),('FEV1_age_adjusted_pairwise_results.csv',pairrows)]:
 with open(args.output/name,'w',encoding='utf-8-sig',newline='') as out:
  dw=csv.DictWriter(out,fieldnames=list(data[0]));dw.writeheader();dw.writerows(data)
(args.output/'figure5b_results.json').write_text(json.dumps({'models':models,'groups':rows,'pairwise':pairrows,'meanage':meanage},indent=2))


# Final plotting order, independent of the contrast orientation in S4 File.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
order = ['DNAH5','ODAD1','DNAH11','DNAAF1','HYDIN','RSPH4A','CCDC40']
by_gene = {row['gene']:row for row in rows}
est = np.array([by_gene[g]['adjusted_mean_at_common_age'] for g in order])
lower = np.array([by_gene[g]['CI_low'] for g in order])
upper = np.array([by_gene[g]['CI_high'] for g in order])
plt.rcParams.update({'pdf.fonttype':42})
fig, ax = plt.subplots(figsize=(7.2,4.5),layout='constrained')
ax.errorbar(range(7),est,yerr=[est-lower,upper-est],fmt='o',color='#527EAA',capsize=0)
ax.set_xticks(range(7),[f'$\\it{{{g}}}$\n(n = {by_gene[g]["n"]})' for g in order],fontsize=9)
ax.set(ylim=(0,122.4),yticks=np.arange(0,121,20),ylabel='Age-adjusted FEV₁ (% predicted)')
ax.spines[['top','right']].set_visible(False)
fig.savefig(args.output/'Figure_5B_age_adjusted.pdf')

# S3 File also reports contrasts from a group-only model and the slope interaction.
from comparisons import holm
G=X[:,:7]
br=np.linalg.lstsq(G,y,rcond=None)[0]
resid=y-G@br
raw_df=len(y)-7
rss=resid@resid
cov_raw=np.linalg.inv(G.T@G)*rss/raw_df
raw_F=((sum((y-y.mean())**2)-rss)/6)/(rss/raw_df)
XI=design(a,interaction=True)
resid_i=y-XI@np.linalg.lstsq(XI,y,rcond=None)[0]
df_i=len(y)-XI.shape[1]
F_i=((e@e-resid_i@resid_i)/6)/(resid_i@resid_i/df_i)
models['group_only']={'global_F':raw_F,'global_p':f.sf(raw_F,6,raw_df),'df_resid':raw_df}
models['age_by_gene']={'global_F':F_i,'global_p':f.sf(F_i,6,df_i),'df_resid':df_i}
for row in pairrows:
    z=(vec(row['gene_1'])-vec(row['gene_2']))[:7]
    delta=z@br;se=np.sqrt(z@cov_raw@z);width=t.ppf(.975,raw_df)*se
    row.update(raw_difference=delta,raw_lower=delta-width,raw_upper=delta+width,raw_p=2*t.sf(abs(delta/se),raw_df))
for row,q in zip(pairrows,holm([r['raw_p'] for r in pairrows])):
    row['raw_p_Holm']=q
with (args.output/'S3_all_contrasts.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(pairrows[0]));writer.writeheader();writer.writerows(pairrows)
(args.output/'figure5b_results.json').write_text(json.dumps({'models':models,'groups':rows,'pairwise':pairrows,'meanage':meanage},indent=2))
