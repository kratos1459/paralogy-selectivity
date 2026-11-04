"""Figure: noise ceiling vs model performance, and anchor calibration (run from repo root)."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
c=pd.read_csv('results/ceiling_vs_models.csv'); v=pd.read_csv('results/eb_validation_by_pair.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
BLUE,ORANGE,TEAL,GREY='#2f6db5','#e8833a','#17a398','#8a8f98'
fig,ax=plt.subplots(1,2,figsize=(7.2,3.3),dpi=300)
a=ax[0]; lim=(-0.15,1.0); a.plot(lim,lim,color=GREY,lw=1,ls='--'); a.text(0.62,0.7,'ceiling',color=GREY,fontsize=8,rotation=38)
a.scatter(c.ceiling,c.r2_scaffold,s=26,facecolors='none',edgecolors=ORANGE,lw=1.2,label='Scaffold split')
a.scatter(c.ceiling,c.r2_lpo,s=26,color=BLUE,label='Publication hold-out')
a.set_xlim(lim); a.set_ylim(lim); a.set_xlabel('Noise ceiling on R$^2$ (1 - V$_x$/V$_{tot}$)'); a.set_ylabel('Model R$^2$'); a.legend(frameon=False,fontsize=8,loc='upper left'); a.set_title('a',loc='left',fontweight='bold')
b=ax[1]; ks=[0,1,2,3,5]
def stat(col): x=v[col].dropna(); return x.median(),x.quantile(.25),x.quantile(.75)
off=[stat('r0')]+[stat(f'off{k}') for k in ks[1:]]; eb=[stat('r0')]+[stat(f'eb{k}') for k in ks[1:]]
for s,col,lab in ((off,GREY,'Plain offset'),(eb,TEAL,'Empirical-Bayes offset')):
    m=[x[0] for x in s]; b.plot(ks,m,marker='o',ms=4,color=col,lw=1.8,label=lab); b.fill_between(ks,[x[1] for x in s],[x[2] for x in s],color=col,alpha=.15,lw=0)
b.set_xlabel('Measured anchors per new publication'); b.set_ylabel('Pooled R$^2$ (median, IQR over 47 pairs)'); b.set_xticks(ks); b.legend(frameon=False,fontsize=8,loc='lower right'); b.set_title('b',loc='left',fontweight='bold')
fig.tight_layout(); fig.savefig('results/figure_ceiling_and_calibration.png'); print('saved')
