"""Pocket-divergence statistics (Fig 4) and regenerated Figures 3 and 5 from released-pipeline results. Run from repo root:
python code/21_pocket_stats_and_main_figures.py r20_dir  (r20_dir contains pocket_divergence_by_pair.csv)"""
import sys,json,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr,rankdata,linregress
R20=sys.argv[1]; rng=np.random.default_rng(1)
FAM={}
for f,ps in {'Kinase':'EGFR/ERBB2 CDK4/CDK6 CDK2/CDK1 JAK1/JAK2 JAK2/TYK2 PIK3CA/PIK3CD MAPK14/MAPK11 ABL1/SRC BRAF/RAF1 AURKA/AURKB FLT3/KIT KDR/FGFR1 ROCK1/ROCK2 GSK3B/GSK3A CHEK1/CHEK2',
 'Epigenetic':'HDAC1/HDAC6 HDAC1/HDAC2 HDAC1/HDAC8 BRD4/BRD2 EP300/CREBBP EZH2/EZH1',
 'Enzyme':'MAOA/MAOB CA2/CA9 CA2/CA1 CA2/CA12 PTGS1/PTGS2 PDE4D/PDE4B PDE4D/PDE5A ACHE/BCHE NOS1/NOS2 NOS1/NOS3 PARP1/PARP2 HSP90AA1/HSP90AB1 DPP4/DPP8 BACE1/BACE2 MMP2/MMP9 MMP13/MMP1 CTSK/CTSS F10/F2',
 'GPCR':'DRD2/DRD3 HTR2A/HTR2C OPRM1/OPRD1 OPRM1/OPRK1 CHRM1/CHRM3 CHRM2/CHRM3 ADRA1A/ADRA1B ADRB1/ADRB2 CNR1/CNR2 HRH3/HRH4','Ion channel':'SCN9A/SCN5A CACNA1C/CACNA1D'}.items():
    for p in ps.split(): FAM[p]=f
pd.DataFrame({'pair':list(FAM),'family':list(FAM.values())}).to_csv('results/family_map.csv',index=False)
COL={'Kinase':'#2f6db5','Enzyme':'#e8833a','GPCR':'#8a5fbf','Epigenetic':'#17a398','Ion channel':'#c9a227'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False})
M=pd.read_csv('results/models_by_pair.csv'); V=pd.read_csv('results/varcomp_by_pair.csv'); M['family']=M.pair.map(FAM); V['family']=V.pair.map(FAM)
# ---------------- Figure 4: pocket divergence
Pd=pd.read_csv(f'{R20}/pocket_divergence_by_pair.csv').merge(V[['pair','Vtot','family']],on='pair').merge(M[['pair','r2_lpo','r2_within']],on='pair',how='left')
Pd['sd']=np.sqrt(Pd.Vtot); Pd.to_csv('results/pocket_divergence_by_pair.csv',index=False)
def rho(x,y): return spearmanr(x,y)[0]
def perm(x,y,fam,n=10000):
    o=abs(rho(x,y)); c=0; fam=np.asarray(fam); x=np.asarray(x)
    for _ in range(n):
        xp=x.copy()
        for f in set(fam):
            i=np.where(fam==f)[0]; xp[i]=rng.permutation(x[i])
        c+=abs(rho(xp,y))>=o
    return (c+1)/(n+1)
def boot(x,y,n=2000):
    x=np.asarray(x);y=np.asarray(y); v=[rho(x[i],y[i]) for i in (rng.integers(0,len(x),len(x)) for _ in range(n))]; return np.nanpercentile(v,[2.5,97.5])
def partial(x,y,z):
    rx,ry,rz=rankdata(x),rankdata(y),rankdata(z)
    ex=rx-np.polyval(np.polyfit(rz,rx,1),rz); ey=ry-np.polyval(np.polyfit(rz,ry,1),rz); return np.corrcoef(ex,ey)[0,1]
S={'n_pairs':len(Pd),'pocket_div_median_by_family':Pd.groupby('family').pocket_div.median().round(3).to_dict(),'seq_div_median':float(Pd.seq_div.median())}
S['sd']={'rho':rho(Pd.pocket_div,Pd.sd),'ci':list(boot(Pd.pocket_div,Pd.sd)),'perm_p':perm(Pd.pocket_div,Pd.sd,Pd.family),'partial_seqdiv':partial(Pd.pocket_div,Pd.sd,Pd.seq_div)}
L=Pd.dropna(subset=['r2_lpo']); S['lpo']={'n':len(L),'rho':rho(L.pocket_div,L.r2_lpo),'ci':list(boot(L.pocket_div,L.r2_lpo))}
S['within']={'rho':rho(L.pocket_div,L.r2_within)}
S['pocket_vs_seqdiv']=rho(Pd.pocket_div,Pd.seq_div)
json.dump(S,open('results/pocket_divergence_summary.json','w'),indent=1); print(json.dumps(S,indent=1))
fig,ax=plt.subplots(1,3,figsize=(6.5,2.4),dpi=300)
def sc(a,x,y,d,xl,yl,t):
    for f,g in d.groupby('family'): a.scatter(g[x],g[y],s=18,color=COL[f],label=f,edgecolor='white',lw=.4)
    k=d.dropna(subset=[x,y]); s=linregress(k[x],k[y]); xs=np.array([k[x].min(),k[x].max()]); a.plot(xs,s.intercept+s.slope*xs,'--',color='#555',lw=1)
    a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(t,loc='left',fontweight='bold')
sc(ax[0],'pocket_div','sd',Pd,'Pocket divergence','SD of same-publication $\\Delta$','a'); sc(ax[1],'pocket_div','r2_lpo',L,'Pocket divergence','LPO R$^2$','b'); sc(ax[2],'pocket_div','seq_div',Pd,'Pocket divergence','Sequence divergence','c')
ax[0].legend(frameon=False,fontsize=6,loc='upper left'); fig.tight_layout(); fig.savefig('results/figure7_pocket_divergence.png')
# ---------------- Figure 3: scaffold vs LPO (+ variance share)
L2=M.dropna(subset=['r2_lpo','r2_scaffold']); fig,axs=plt.subplots(2,2,figsize=(6.5,5.0),dpi=300); ax=axs.ravel()
a=ax[0]; a.plot([-.3,1],[-.3,1],'--',color='#8a8f98',lw=1)
for f,g in L2.groupby('family'): a.scatter(g.r2_scaffold,g.r2_lpo,s=18,color=COL[f],edgecolor='white',lw=.4,label=f)
a.set_xlabel('Scaffold-split R$^2$'); a.set_ylabel('Leave-publication-out R$^2$'); a.set_title('a',loc='left',fontweight='bold'); a.legend(frameon=False,fontsize=6)
b=ax[1]; dat=[M.r2_scaffold.dropna(),M.r2_lpo.dropna(),M.r2_within.dropna()]; bp=b.boxplot(dat,widths=.55,patch_artist=True,showfliers=False)
for p,c in zip(bp['boxes'],['#e8833a','#2f6db5','#17a398']): p.set_facecolor(c); p.set_alpha(.6)
for i,d in enumerate(dat,1): b.scatter(rng.normal(i,.06,len(d)),d,s=6,color='#333',alpha=.5)
b.set_xticks([1,2,3]); b.set_xticklabels(['Scaffold','LPO','Within series']); b.set_ylabel('R$^2$'); b.set_title('b',loc='left',fontweight='bold'); b.axhline(0,color='#8a8f98',lw=.8)
c=ax[2]; k=M.dropna(subset=['r2_lpo','n_docs_lpo']); c.scatter(k.n_docs_lpo,k.r2_lpo,s=18,color='#2f6db5',edgecolor='white',lw=.4); c.set_xscale('log'); c.set_xlabel('Publications per pair (LPO)'); c.set_ylabel('LPO R$^2$'); c.set_title('c',loc='left',fontweight='bold')
S3={'spearman_lpo_vs_ndocs':[float(x) for x in spearmanr(k.n_docs_lpo,k.r2_lpo)],'within_series_spearman_median':float(M.rho_within.median())}
Vv=V.dropna(subset=['share']); d_=ax[3]; o=Vv.sort_values('share'); d_.bar(range(len(o)),o.share,color='#2f6db5',width=.85,label='Between publications'); d_.bar(range(len(o)),1-o.share,bottom=o.share,color='#d9dee6',width=.85,label='Within publications')
d_.set_xlabel('Pairs (sorted by share)'); d_.set_ylabel('Share of selectivity variance'); d_.set_xticks([]); d_.legend(frameon=False,fontsize=6,loc='lower right'); d_.set_title('d',loc='left',fontweight='bold')
json.dump(S3,open('results/figure3_stats.json','w')); print(S3)
fig.tight_layout(); fig.savefig('results/figure3_evaluation_regimes.png')
# ---------------- Figure 6: training formulations
TF=pd.concat([pd.read_csv(f'results/training_formulations_part{i}.csv') for i in (0,1)]).drop_duplicates(['pair','scheme','model'])
TF.to_csv('results/training_formulations_by_pair.csv',index=False)
schemes=['standard','series_centred','pairwise']; models=['rf','xgb','mlp']; mc={'rf':'#2f6db5','xgb':'#e8833a','mlp':'#8a5fbf'}; ml={'rf':'RF','xgb':'XGBoost','mlp':'MLP'}
fig,ax=plt.subplots(1,2,figsize=(6.5,2.7),dpi=300); b=ax[0]
for j,m in enumerate(models):
    xs=np.arange(3)+(j-1)*.26; sub=[TF[(TF.scheme==s)&(TF.model==m)] for s in schemes]; vals=np.array([x.w.median() for x in sub])
    lo=np.array([x.w.quantile(.25) for x in sub]); hi=np.array([x.w.quantile(.75) for x in sub])
    b.bar(xs,vals,width=.25,color=mc[m],alpha=.85,label=ml[m]); b.errorbar(xs,vals,yerr=[vals-lo,hi-vals],fmt='none',ecolor='#333',lw=.7)
    b.scatter(xs,[x.rho.median() for x in sub],marker='D',s=12,color='white',edgecolor='#333',zorder=3)
b.set_xticks(range(3)); b.set_xticklabels(['Standard','Series-centred','Pairwise']); b.axhline(0,color='#8a8f98',lw=.8); b.set_ylabel('Within-series R$^2$ (bars, median and IQR)\nSpearman $\\rho$ (diamonds)',fontsize=7); b.legend(frameon=False,fontsize=6); b.set_title('a',loc='left',fontweight='bold')
c=ax[1]; xx=np.linspace(0,.5,50); c.fill_between(xx,2*xx,1.0,color='#f3d9d4',alpha=.6); c.text(.02,.9,'k > 2r: R$^2_{within}$ < 0',fontsize=6,color='#a33')
mk_={'standard':'o','series_centred':'s','pairwise':'^'}
for s in schemes:
    for m in models:
        x=TF[(TF.scheme==s)&(TF.model==m)]; c.scatter(x.r.median(),x.k.median(),s=34,marker=mk_[s],color=mc[m],edgecolor='white',lw=.5)
c.plot(xx,2*xx,'--',color='#555',lw=1); c.set_xlim(0,.4); c.set_ylim(0,1.0); c.set_xlabel('Centred correlation r'); c.set_ylabel('Dispersion ratio k = $\\sigma$(pred)/$\\sigma$(obs)'); c.set_title('b',loc='left',fontweight='bold')
from matplotlib.lines import Line2D
c.legend([Line2D([],[],marker=mk_[s],color='#555',ls='') for s in schemes],['Standard','Series-centred','Pairwise'],frameon=False,fontsize=6,loc='lower right')
fig.tight_layout(); fig.savefig('results/figure6_training_formulations.png'); print('figs saved')
