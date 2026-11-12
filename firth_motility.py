import numpy as np
from scipy.special import expit
from scipy.optimize import minimize, brentq
from scipy.stats import chi2

def fit(X,y,ci=True):
    X=np.asarray(X,float);y=np.asarray(y,float)
    def calc(b):
        eta=X@b;p=expit(eta);w=np.maximum(expit(eta)*expit(-eta),1e-100)
        info=X.T@(w[:,None]*X)
        sign,ld=np.linalg.slogdet(info)
        if sign<=0:return 1e100,np.zeros(len(b))
        if np.linalg.cond(info)>1e12:return 1e100,np.zeros(len(b))
        h=w*np.sum((X@np.linalg.inv(info))*X,axis=1)
        ll=np.sum(y*eta-np.logaddexp(0,eta))+.5*ld
        grad=X.T@(y-p+h*(.5-p))
        return -ll,-grad
    opt=minimize(calc,np.zeros(X.shape[1]),jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':1000})
    b=opt.x
    if np.max(np.abs(calc(b)[1]))>2e-5:raise RuntimeError(('score',opt.message,b,calc(b)[1]))
    def profile(v):
        free=[i for i in range(len(b)) if i!=1]
        def f(z):
            bb=b.copy();bb[1]=v;bb[free]=z
            val,g=calc(bb);return val,g[free]
        rr=minimize(f,b[free],jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':500})
        if np.max(np.abs(f(rr.x)[1]))>3e-5:
            start=minimize(lambda z:f(z)[0],b[free],method="Powell",options={"xtol":1e-10,"ftol":1e-12}).x
            rr=minimize(f,start,jac=True,method="BFGS",options={"gtol":1e-8})
        if np.max(np.abs(f(rr.x)[1]))>3e-5:raise RuntimeError(('profile',v,rr.message))
        return rr.fun
    pval=chi2.sf(max(0,2*(profile(0)-opt.fun)),1)
    result={'or':float(np.exp(b[1])),'p':float(pval),'score':float(np.max(np.abs(calc(b)[1])))}
    if ci:
        cutoff=chi2.ppf(.95,1)
        def root(v):return 2*(profile(v)-opt.fun)-cutoff
        limits=[]
        for direction in [-1,1]:
            step=1.
            while root(b[1]+direction*step)<0:
                step*=1.6
                if step>200:raise RuntimeError(('unbounded interval',direction,b.tolist()))
            a,c=sorted([b[1],b[1]+direction*step]);limits.append(float(np.exp(brentq(root,a,c,xtol=1e-8))))
        result.update(lo=limits[0],hi=limits[1])
    return result

def design(data,gene,phenotype,age='enrollment_age',adjust=True):
    data=[r for r in data if r['phenotypes'][phenotype] is not None and (not adjust or (r[age] is not None and r['sex'] in ['M','F']))]
    xx=[[1,int(gene in r['genes'])] for r in data]
    if adjust:
        aa=np.array([r[age] for r in data]);aa=(aa-aa.mean())/aa.std()
        xx=[row+[float(a),int(r['sex']=='M')] for row,a,r in zip(xx,aa,data)]
    xx=np.array(xx)
    # Constant covariates are not identifiable alongside an intercept.
    keep=[0,1]+[j for j in range(2,xx.shape[1]) if np.ptp(xx[:,j])>0]
    return xx[:,keep],np.array([r['phenotypes'][phenotype] for r in data]),data

