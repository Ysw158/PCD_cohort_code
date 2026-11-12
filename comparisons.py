"""Fig 6C–D and S5 Fig. Run with private a.json and an output directory."""
import csv,itertools,json,sys
from pathlib import Path
from collections import defaultdict
from scipy.stats import mannwhitneyu,fisher_exact
from scipy.stats.contingency import odds_ratio

def holm(p):
    order=sorted(range(len(p)),key=p.__getitem__)
    adjusted=p.copy(); max_p=0
    for rank,i in enumerate(order):
        max_p=max(max_p,min(1,(len(p)-rank)*p[i]))
        adjusted[i]=max_p
    return adjusted

def write(path,header,rows):
    with path.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(header);writer.writerows(rows)

def registry_comparison(out):
    result=[]
    for gene,n,m in [('FOXJ1',3,1),('DNAAF6',5,11),('DNAAF1',7,18),('CCNO',9,31),('ODAD1',10,37),('HYDIN',9,42)]:
        table=[[n,125-n],[m,1236-m]]
        estimate=odds_ratio(table,kind='conditional');ci=estimate.confidence_interval()
        result.append([gene,n,m,estimate.statistic,ci.low,ci.high,fisher_exact(table).pvalue])
    write(out/'s5_figure.csv',['gene','Chinese_n_of_125','registry_n_of_1236','conditional_OR','lower','upper','two_sided_p'],result)


if __name__=='__main__':
    out = Path(sys.argv[2])
    out.mkdir(exist_ok=True, parents=True)
    registry_comparison(out)
    if sys.argv[1] == "--registry-only":
        sys.exit(0)
    data = json.loads(Path(sys.argv[1]).read_text())
    for panel in ['6c','6d']:
        groups=defaultdict(list)
        for r in data:
            group='+'.join(sorted(r['genes'])) if panel=='6c' else r['cbp']
            if r['cbf'] is not None and group:groups[group].append(r['cbf'])
        rows=[]
        for g,h in itertools.combinations(sorted(g for g in groups if len(groups[g])>=3),2):
            a,b=groups[g],groups[h]
            p=mannwhitneyu(a,b,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue
            rows.append([g,h,len(a),len(b),p])
        assert len(rows)=={'6c':15,'6d':6}[panel]
        for row,q in zip(rows,holm([r[-1] for r in rows])):row.append(q)
        write(out/f'figure{panel}.csv',['group1','group2','n1','n2','p','p_holm'],rows)
