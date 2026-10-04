#!/usr/bin/env python3
"""
Sudan measles surveillance figure (WHO monthly provisional data).
Used to justify the seasonal forcing empirically, NOT to claim a full fit.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans',
                     'axes.labelsize':12,'axes.titlesize':11,
                     'figure.dpi':300,'savefig.dpi':300,
                     'savefig.bbox':'tight','axes.grid':True,'grid.alpha':0.3})

from pathlib import Path
CSV = Path(__file__).resolve().parent.parent / "data" / "sudan_measles_monthly_2012_2022.csv"
COL = "measles_total"
s=pd.read_csv(CSV, parse_dates=["date"]).rename(columns={"date":"t"}).sort_values("t")
full=pd.date_range(s.t.min(),s.t.max(),freq="MS")
s=s.set_index("t").reindex(full)
y=s[COL].interpolate(limit_direction="both").to_numpy(float)
d=full

fig,ax=plt.subplots(1,3,figsize=(16,4.6))

# (a) monthly series
ax[0].plot(d,y,color='#2c3e50',lw=1.3)
ax[0].fill_between(d,y,alpha=.18,color='#2c3e50')
ax[0].set_xlabel('Year'); ax[0].set_ylabel('Reported measles cases / month')
ax[0].set_title('(a) Sudan monthly measles, 2012--2022 (WHO)')
ax[0].tick_params(axis='x',rotation=0)

# (b) seasonal profile + fitted sinusoid
m=pd.DataFrame({'y':y,'m':d.month}).groupby('m').y.agg(['mean','sem'])
mm=np.arange(1,13)
ax[1].bar(mm,m['mean'],yerr=m['sem'],color='#3498db',alpha=.75,
          edgecolor='white',capsize=3,label='Monthly mean $\\pm$ s.e.')
# least-squares sinusoid  a + b*sin(2pi(m-phi)/12)
A=np.column_stack([np.ones(12),np.sin(2*np.pi*mm/12),np.cos(2*np.pi*mm/12)])
c,*_=np.linalg.lstsq(A,m['mean'].to_numpy(),rcond=None)
amp=np.hypot(c[1],c[2]); b1=amp/c[0]
mf=np.linspace(1,12,200)
Af=np.column_stack([np.ones(200),np.sin(2*np.pi*mf/12),np.cos(2*np.pi*mf/12)])
ax[1].plot(mf,Af@c,'r-',lw=2.4,
           label=f'Sinusoid: $\\beta_1\\approx{b1:.2f}$')
ax[1].set_xticks(mm); ax[1].set_xlabel('Calendar month')
ax[1].set_ylabel('Mean reported cases')
ax[1].set_title('(b) Seasonal profile and fitted forcing')
ax[1].legend(fontsize=9,framealpha=.9)

# (c) variance decomposition
dd=pd.DataFrame({'y':y},index=d); dd['month']=dd.index.month; dd['year']=dd.index.year
sst=((y-y.mean())**2).sum()
seas=1-((y-dd.groupby('month').y.transform('mean'))**2).sum()/sst
yr  =1-((y-dd.groupby('year').y.transform('mean'))**2).sum()/sst
res =max(0.0,1-seas-yr)
vals=[seas*100,yr*100,res*100]
lbl=['Seasonal\n(within-year)','Interannual\n(between-year)','Residual']
col=['#3498db','#e74c3c','#95a5a6']
bars=ax[2].bar(lbl,vals,color=col,alpha=.85,edgecolor='white')
for b,v in zip(bars,vals):
    ax[2].text(b.get_x()+b.get_width()/2,v+1.2,f'{v:.0f}%',ha='center',
               fontweight='bold')
ax[2].set_ylabel('Share of total variance (%)')
ax[2].set_ylim(0,max(vals)*1.28)
ax[2].set_title('(c) Variance decomposition')

fig.tight_layout(); fig.savefig('FIGG_sudan.pdf'); plt.close(fig)
print(f"beta1 estimate from data : {b1:.3f}")
print(f"peak month {int(m['mean'].idxmax())}, trough {int(m['mean'].idxmin())}, "
      f"ratio {m['mean'].max()/m['mean'].min():.2f}")
print(f"seasonal {seas:.1%}, interannual {yr:.1%}, residual {res:.1%}")
print(f"total cases {y.sum():.0f} over {len(y)} months")
