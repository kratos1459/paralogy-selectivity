"""Figure 7: curated vs random pairs. Run from repo root."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
C=pd.read_csv('results/models_by_pair.csv'); R=pd.read_csv('results/random_pairs/models_by_pair.csv')
EC=pd.read_csv('results/eb_validation_by_pair.csv'); ER=pd.read_csv('results/random_pairs/eb_validation_by_pair.csv')
VC=pd.read_csv('results/varcomp_by_pair.csv'); VR=pd.read_csv('results/random_pairs/varcomp_by_pair.csv')
cc,cr='#1F5F8B','#E07A1F'
fig,ax=plt.subplots(1,3,figsize=(13,4.2),dpi=200)
a=ax[0]
for i,(lab,col,x) in enumerate([('Curated',cc,C),('Random',cr,R)]):
    for j,k in enumerate(['r2_scaffold','r2_lpo']):
        v=x[k].dropna().clip(-1.5,1); p=i*3+j
        a.boxplot(v,positions=[p],widths=.6,patch_artist=True,boxprops=dict(facecolor=col,alpha=.35 if j else .8,edgecolor=col),medianprops=dict(color='k'),showfliers=False)
        a.scatter(np.random.normal(p,.07,len(v)),v,s=7,color=col,alpha=.6)
a.set_xticks([0,1,3,4]);a.set_xticklabels(['Scaffold','LPO','Scaffold','LPO']);a.axhline(0,c='grey',lw=.8)
a.text(.5,-2.05,'Curated (n=%d)'%C.r2_lpo.notna().sum(),ha='center',color=cc);a.text(3.5,-2.05,'Random (n=%d)'%R.r2_lpo.notna().sum(),ha='center',color=cr)
a.set_ylabel('R² (clipped at −1.5)');a.set_title('a  Evaluation regime',loc='left',fontweight='bold');a.set_ylim(-1.6,1);a.tick_params(axis='x',pad=2)
b=ax[1]
b.hist(VC.share.dropna(),bins=np.linspace(0,1,15),alpha=.6,color=cc,label='Curated');b.hist(VR.share.dropna(),bins=np.linspace(0,1,15),alpha=.6,color=cr,label='Random')
b.set_xlabel('Between-publication variance share');b.set_ylabel('Pairs');b.legend(frameon=False);b.set_title('b  Variance structure',loc='left',fontweight='bold')
c=ax[2]; ms=[1,2,3,5]
for lab,E,col in [('Curated',EC,cc),('Random',ER,cr)]:
    c.plot(ms,[E[f'off{m}'].median() for m in ms],'--o',color=col,label=f'{lab} plain');c.plot(ms,[E[f'eb{m}'].median() for m in ms],'-o',color=col,label=f'{lab} EB')
c.axhline(0,c='grey',lw=.8);c.set_xlabel('Anchor compounds (m)');c.set_ylabel('Median R² (LPO)');from matplotlib.lines import Line2D;c.legend([Line2D([],[],c=cc,marker='o'),Line2D([],[],c=cr,marker='o'),Line2D([],[],c='k',ls='--'),Line2D([],[],c='k')],['Curated','Random','Plain offset','Empirical Bayes'],frameon=False,fontsize=8);c.set_title('c  Anchor calibration',loc='left',fontweight='bold')
for s in ax:s.spines[['top','right']].set_visible(False)
plt.tight_layout();plt.savefig('results/figure_random_pair_replication.png');
print(C.r2_scaffold.median(),C.r2_lpo.median(),R.r2_scaffold.median(),R.r2_lpo.median(),VC.share.median(),VR.share.median())
print(EC[[f'{p}{m}' for p in('off','eb') for m in ms]].median().round(3).to_dict())
