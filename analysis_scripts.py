# Analysis scripts used for the manuscript (Python 3; SciPy, pingouin, scikit-posthocs, matplotlib)
# Run order: stats.py -> dunn.py -> extra.py -> figs.py (input: bergenin_dataset.xlsx)
-e 
# ======== stats.py ========
import pandas as pd, numpy as np, pingouin as pg, json
from scipy import stats
f='/mnt/user-data/uploads/1792186938230_bergenin_dataset.xlsx'
fbg=pd.read_excel(f,'FBG_raw'); og=pd.read_excel(f,'OGTT_raw'); bw=pd.read_excel(f,'BodyWeight_raw')
for d in (fbg,og,bw): d['Group']=d['Group'].replace({'Glibenclamide 5 mg/kg':'Glibenclamide 5 mg/kg'})
order=['Normal control','Diabetic control','Glibenclamide 5 mg/kg','Bergenin 25 mg/kg','Bergenin 50 mg/kg']
out={}
def oneway(df,col,label):
    print('\n=====',label,col)
    groups=[df[df.Group==g][col].values for g in order]
    sw={g:stats.shapiro(df[df.Group==g][col]).pvalue for g in order}
    lev=stats.levene(*groups,center='median')
    print('Shapiro min p',min(sw.values()).round(3),'Levene p',round(lev.pvalue,4))
    an=pg.anova(data=df,dv=col,between='Group',detailed=True); print(an.round(4).to_string())
    kw=stats.kruskal(*groups); print('KW',kw)
    tk=pg.pairwise_tukey(data=df,dv=col,between='Group')
    tk['pair']=tk.A+' vs '+tk.B
    print(tk[['pair','diff','p_tukey','hedges']].round(4).to_string())
    gh=pg.pairwise_gameshowell(data=df,dv=col,between='Group'); gh['pair']=gh.A+' vs '+gh.B
    print('GamesHowell'); print(gh[['pair','diff','pval']].round(4).to_string())
    return an,tk,gh
for c in ['Baseline (mg/dL)','Day 3 (mg/dL)','Day 10 (mg/dL)','Day 17 pre-OGTT (mg/dL)']: oneway(fbg,c,'FBG')
for c in ['0 min (mg/dL)','AUC (mg·min/dL)']: oneway(og,c,'OGTT')
for c in ['Baseline (g)','Day 3 (g)','Day 17 (g)']: oneway(bw,c,'BW')
# mixed anova
for name,df,cols in [('FBG',fbg,['Baseline (mg/dL)','Day 3 (mg/dL)','Day 10 (mg/dL)','Day 17 pre-OGTT (mg/dL)']),('BW',bw,['Baseline (g)','Day 3 (g)','Day 10 (g)','Day 17 (g)'])]:
    l=df.melt(id_vars=['Group','Animal_ID'],value_vars=cols,var_name='Time',value_name='v')
    m=pg.mixed_anova(data=l,dv='v',within='Time',subject='Animal_ID',between='Group',correction=True)
    print('\nMixed ANOVA',name); print(m.round(5).to_string())
# BW change
bw['d']=bw['Day 17 (g)']-bw['Day 3 (g)']; bw['pct']=100*(bw['Day 17 (g)']-bw['Baseline (g)'])/bw['Baseline (g)']
print(bw.groupby('Group',sort=False)[['d','pct']].agg(['mean','std']).round(2))
oneway(bw,'pct','BW pct'); oneway(bw,'d','BW d3-17')
# paired within-group Day3 -> Day17 FBG
for g in order:
    s=fbg[fbg.Group==g]; print(g, stats.ttest_rel(s['Day 3 (mg/dL)'],s['Day 17 pre-OGTT (mg/dL)']).pvalue.round(5), 'D3-D17 drop', (s['Day 3 (mg/dL)']-s['Day 17 pre-OGTT (mg/dL)']).mean().round(1))
# percent reductions
dc17=fbg[fbg.Group=='Diabetic control']['Day 17 pre-OGTT (mg/dL)'].mean(); dcA=og[og.Group=='Diabetic control']['AUC (mg·min/dL)'].mean()
for g in order[2:]:
    print(g,'FBG17 red vs DC %',100*(1-fbg[fbg.Group==g]['Day 17 pre-OGTT (mg/dL)'].mean()/dc17),'AUC red',100*(1-og[og.Group==g]['AUC (mg·min/dL)'].mean()/dcA))
# AUC recompute check
t=np.array([0,30,60,90,120]); cols=['0 min (mg/dL)','30 min (mg/dL)','60 min (mg/dL)','90 min (mg/dL)','120 min (mg/dL)']
print('AUC check max diff',np.abs(np.trapezoid(og[cols].values,t,axis=1)-og['AUC (mg·min/dL)']).max())
# ogtt 2-way
l=og.melt(id_vars=['Group','Animal_ID'],value_vars=cols,var_name='Time',value_name='v')
print(pg.mixed_anova(data=l,dv='v',within='Time',subject='Animal_ID',between='Group',correction=True).round(5).to_string())
# diabetes induction excl
print((fbg[fbg.Group!='Normal control']['Day 3 (mg/dL)']>200).all())
-e 
# ======== dunn.py ========
import pandas as pd, numpy as np, scikit_posthocs as sp, json
f='/mnt/user-data/uploads/1792186938230_bergenin_dataset.xlsx'
fbg=pd.read_excel(f,'FBG_raw'); og=pd.read_excel(f,'OGTT_raw')
for df,c in [(fbg,'Day 17 pre-OGTT (mg/dL)'),(og,'0 min (mg/dL)'),(og,'AUC (mg·min/dL)'),(fbg,'Day 10 (mg/dL)')]:
    p=sp.posthoc_dunn(df,val_col=c,group_col='Group',p_adjust='bonferroni')
    print('\n',c); print(p.round(4).to_string())
# Mann-Whitney / per-group shapiro on AUC
from scipy import stats
for g,s in og.groupby('Group'): print(g,'Shapiro AUC',round(stats.shapiro(s['AUC (mg·min/dL)']).pvalue,3))
-e 
# ======== extra.py ========
import pandas as pd, numpy as np, pingouin as pg, scikit_posthocs as sp
from scipy import stats
f='/mnt/user-data/uploads/1792186938230_bergenin_dataset.xlsx'
og=pd.read_excel(f,'OGTT_raw'); fbg=pd.read_excel(f,'FBG_raw')
cols=['0 min (mg/dL)','30 min (mg/dL)','60 min (mg/dL)','90 min (mg/dL)','120 min (mg/dL)']; t=np.array([0,30,60,90,120])
og['iAUC']=np.trapezoid(og[cols].values-og[['0 min (mg/dL)']].values,t,axis=1)
og['peakrise']=og[cols[1]]-og[cols[0]]
og['AUCnorm']=og['AUC (mg·min/dL)']/120
print(og.groupby('Group',sort=False)[['iAUC','peakrise','AUCnorm']].agg(['mean','std']).round(1).to_string())
for c in ['iAUC','peakrise']:
    print(c,'Levene',round(stats.levene(*[s[c] for _,s in og.groupby('Group')],center='median').pvalue,3), 'ANOVA p', pg.anova(data=og,dv=c,between='Group')['p_unc'][0].round(4), 'KW',round(stats.kruskal(*[s[c] for _,s in og.groupby('Group')]).pvalue,4))
    print(pg.pairwise_tukey(data=og,dv=c,between='Group').assign(pair=lambda d:d.A+' vs '+d.B)[['pair','diff','p_tukey']].round(4).to_string())
# paired D3->D17 with bonferroni
for g,s in fbg.groupby('Group',sort=False):
    r=stats.ttest_rel(s['Day 3 (mg/dL)'],s['Day 17 pre-OGTT (mg/dL)']); print(g,'t',round(r.statistic,2),'p',r.pvalue,'bonf',min(1,r.pvalue*5))
# group D3->D17 % reduction
for g,s in fbg.groupby('Group',sort=False): print(g, 'pct fall D3->D17', round(100*(1-s['Day 17 pre-OGTT (mg/dL)'].mean()/s['Day 3 (mg/dL)'].mean()),1))
# eta2 for KW (epsilon2)
for c,df in [('Day 17 pre-OGTT (mg/dL)',fbg),('AUC (mg·min/dL)',og),('0 min (mg/dL)',og)]:
    H=stats.kruskal(*[s[c] for _,s in df.groupby('Group')]).statistic; print(c,'H',round(H,2),'eps2',round(H/(35-1)*1,3) if False else round((H-5+1)/(35-5),3))
# glib vs bergenin 50, ber50 vs ber25 Dunn done earlier. Day17 FBG Tukey df etc
print(pg.welch_anova(data=fbg,dv='Day 17 pre-OGTT (mg/dL)',between='Group').round(4).to_string())
print(pg.welch_anova(data=og,dv='AUC (mg·min/dL)',between='Group').round(4).to_string())
-e 
# ======== figs.py ========
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, scikit_posthocs as sp, pingouin as pg
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.8,'figure.dpi':300,'savefig.dpi':300})
f='/mnt/user-data/uploads/1792186938230_bergenin_dataset.xlsx'
fbg=pd.read_excel(f,'FBG_raw'); og=pd.read_excel(f,'OGTT_raw'); bw=pd.read_excel(f,'BodyWeight_raw')
order=['Normal control','Diabetic control','Glibenclamide 5 mg/kg','Bergenin 25 mg/kg','Bergenin 50 mg/kg']
lab={'Normal control':'Normal control','Diabetic control':'Diabetic control','Glibenclamide 5 mg/kg':'Glibenclamide 5 mg/kg','Bergenin 25 mg/kg':'Bergenin 25 mg/kg','Bergenin 50 mg/kg':'Bergenin 50 mg/kg'}
col={'Normal control':'#4D4D4D','Diabetic control':'#D55E00','Glibenclamide 5 mg/kg':'#0072B2','Bergenin 25 mg/kg':'#E5B300','Bergenin 50 mg/kg':'#009E73'}
mk={'Normal control':'o','Diabetic control':'s','Glibenclamide 5 mg/kg':'^','Bergenin 25 mg/kg':'D','Bergenin 50 mg/kg':'v'}
sem=lambda x: x.std(ddof=1)/np.sqrt(len(x))
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def tukey_p(df,c):
    t=pg.pairwise_tukey(data=df,dv=c,between='Group'); d={}
    for _,r in t.iterrows(): d[(r.A,r.B)]=d[(r.B,r.A)]=r.p_tukey
    return d
def dunn_p(df,c): return sp.posthoc_dunn(df,val_col=c,group_col='Group',p_adjust='bonferroni')
def tl(ax): ax.text(-0.14,1.06,ax._lab,transform=ax.transAxes,fontsize=12,fontweight='bold')
def panel(ax,l): ax._lab=l; tl(ax)
def dotbar(ax,df,c,pm,ylabel,title=None):
    for i,g in enumerate(order):
        v=df[df.Group==g][c].values
        ax.bar(i,v.mean(),color=col[g],alpha=0.35,edgecolor=col[g],linewidth=1.2,width=0.7)
        ax.errorbar(i,v.mean(),yerr=sem(pd.Series(v)),color='k',capsize=3,lw=1)
        ax.scatter(i+np.linspace(-0.18,0.18,len(v)),v,s=14,color=col[g],edgecolor='k',linewidth=0.4,zorder=3)
    ax.set_xticks(range(5)); ax.tick_params(axis='x',labelsize=7); ax.set_xticklabels(['Normal\ncontrol','Diabetic\ncontrol','Glib\n5 mg/kg','Berg\n25 mg/kg','Berg\n50 mg/kg'])
    ax.set_ylabel(ylabel)
    ymax=max(df[c])
    # stars vs diabetic control (*) and vs normal control (#)
    for i,g in enumerate(order):
        s=''
        if g!='Diabetic control' and g!='Normal control':
            s+=star(pm(g,'Diabetic control'))
        if s: ax.text(i,df[df.Group==g][c].max()+0.02*ymax,s,ha='center',fontsize=11)
    ax.text(1,max(df[df.Group=='Diabetic control'][c].max(),0)+0.02*ymax,'###' if pm('Normal control','Diabetic control')<0.001 else '',ha='center',fontsize=9,color='#444')

# ---------- Fig 1
fig,ax=plt.subplots(1,2,figsize=(7.6,3.4),gridspec_kw={'width_ratios':[1.5,1]})
tc=['Baseline (mg/dL)','Day 3 (mg/dL)','Day 10 (mg/dL)','Day 17 pre-OGTT (mg/dL)']; x=[0,3,10,17]
for g in order:
    s=fbg[fbg.Group==g]; m=[s[c].mean() for c in tc]; e=[sem(s[c]) for c in tc]
    ax[0].errorbar(x,m,yerr=e,color=col[g],marker=mk[g],ms=5,lw=1.4,capsize=2.5,label=g)
ax[0].axvspan(3.5,17.5,color='#f1f1f1',zorder=0); ax[0].text(10.5,40,'Daily treatment (Days 4–17)',ha='center',fontsize=7.5,color='#555')
ax[0].set_xticks(x); ax[0].set_xticklabels(['Baseline','Day 3','Day 10','Day 17']); ax[0].set_xlabel('Time'); ax[0].set_ylabel('Fasting blood glucose (mg/dL)'); ax[0].set_ylim(0,460)
d17=dunn_p(fbg,tc[3])
panel(ax[0],'A')
dotbar(ax[1],fbg,tc[3],lambda g,h:d17.loc[g,h],'FBG on Day 17 (mg/dL)'); panel(ax[1],'B')
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=5,frameon=False,fontsize=7.5,bbox_to_anchor=(0.5,-0.01),handlelength=1.5,columnspacing=1.2)
fig.tight_layout(rect=[0,0.07,1,1]); fig.savefig('/home/claude/w/Fig1_FBG.png',bbox_inches='tight'); plt.close()

# ---------- Fig 2
oc=['0 min (mg/dL)','30 min (mg/dL)','60 min (mg/dL)','90 min (mg/dL)','120 min (mg/dL)']; t=[0,30,60,90,120]
fig,ax=plt.subplots(1,2,figsize=(7.6,3.4),gridspec_kw={'width_ratios':[1.5,1]})
for g in order:
    s=og[og.Group==g]; ax[0].errorbar(t,[s[c].mean() for c in oc],yerr=[sem(s[c]) for c in oc],color=col[g],marker=mk[g],ms=5,lw=1.4,capsize=2.5,label=g)
ax[0].set_xticks(t); ax[0].set_xlabel('Time after oral glucose load (min)'); ax[0].set_ylabel('Blood glucose (mg/dL)'); ax[0].set_ylim(0,560); panel(ax[0],'A')
dA=dunn_p(og,'AUC (mg·min/dL)')
dotbar(ax[1],og.assign(AUCk=og['AUC (mg·min/dL)']/1000),'AUCk',lambda g,h:dA.loc[g,h],'OGTT AUC$_{0–120}$ (×10³ mg·min/dL)'); panel(ax[1],'B')
fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',ncol=5,frameon=False,fontsize=7.5,bbox_to_anchor=(0.5,-0.01),handlelength=1.5,columnspacing=1.2)
fig.tight_layout(rect=[0,0.07,1,1]); fig.savefig('/home/claude/w/Fig2_OGTT.png',bbox_inches='tight'); plt.close()

# ---------- Fig 3
bc=['Baseline (g)','Day 3 (g)','Day 10 (g)','Day 17 (g)']
bw['pct']=100*(bw['Day 17 (g)']-bw['Baseline (g)'])/bw['Baseline (g)']
fig,ax=plt.subplots(1,2,figsize=(7.6,3.4),gridspec_kw={'width_ratios':[1.5,1]})
for g in order:
    s=bw[bw.Group==g]; ax[0].errorbar(x,[s[c].mean() for c in bc],yerr=[sem(s[c]) for c in bc],color=col[g],marker=mk[g],ms=5,lw=1.4,capsize=2.5,label=g)
ax[0].axvspan(3.5,17.5,color='#f1f1f1',zorder=0); ax[0].set_xticks(x); ax[0].set_xticklabels(['Baseline','Day 3','Day 10','Day 17']); ax[0].set_xlabel('Time'); ax[0].set_ylabel('Body weight (g)'); ax[0].set_ylim(19.5,26.5); panel(ax[0],'A')
tk=tukey_p(bw,'pct'); dotbar(ax[1],bw,'pct',lambda g,h:tk[(g,h)],'Body-weight change from baseline (%)'); ax[1].axhline(0,color='k',lw=0.6); panel(ax[1],'B')
fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',ncol=5,frameon=False,fontsize=7.5,bbox_to_anchor=(0.5,-0.01),handlelength=1.5,columnspacing=1.2)
fig.tight_layout(rect=[0,0.07,1,1]); fig.savefig('/home/claude/w/Fig3_BW.png',bbox_inches='tight'); plt.close()
