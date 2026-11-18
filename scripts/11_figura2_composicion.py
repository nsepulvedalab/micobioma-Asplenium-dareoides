"""11_figura2_composicion.py
Hacer la Figura 2 de la tesis: composición por género y por especie (porcentaje de lecturas).
Se agrupan las 10 categorías más abundantes y el resto queda como "Otros", igual que en el script 06.
Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv. Salida: resultados/Fig2_composicion.png."""
import os, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
B = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(B, '..', 'datos'))
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 9, 'mathtext.fontset': 'custom',
                     'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic'})
tax = pd.read_csv(os.path.join(D, 'ASV_taxonomy_UNITE97_60.csv'))
asv = pd.read_csv(os.path.join(D, 'ASV_table.csv'))
if 'ASV' not in tax.columns: tax = tax.rename(columns={tax.columns[0]: 'ASV'})
asv = asv.rename(columns={asv.columns[0]: 'Sample'}).set_index('Sample')
d = asv.sum(0).rename('reads').reset_index().rename(columns={'index': 'ASV'}).merge(tax[['ASV', 'taxon']], on='ASV', how='left')
fp = d.taxon.where(~(d.taxon.isna() | (d.taxon == 'Unassigned'))).str.split('|').str[0]
d['Genero'] = fp.str.split('_').str[0].fillna('Unassigned'); d['Especie'] = fp.fillna('Unassigned')
COL = ['#cccccc', '#a71b4b', '#de5925', '#f39b29', '#fbd476', '#fefdbe', '#a5e8ad', '#22c4b3', '#0090b5', '#584b9f', '#8c8c8c']
SUP = ('ales', 'aceae', 'mycota', 'mycetes')  # rangos sobre género: sin cursiva
def etiqueta(cat):
    if cat in ('Unassigned', 'Otros'): return 'No asignado' if cat == 'Unassigned' else 'Otros'
    g, _, sp = cat.partition('_')
    sup = g.endswith(SUP)
    gtxt = g if sup else r'$\it{%s}$' % g
    if not sp: return gtxt
    if sp == 'sp': return gtxt + ' sp.'
    return r'$\it{%s\ %s}$' % (g, sp)
def datos(col):
    s = d.groupby(col).reads.sum().sort_values(ascending=False)
    s = pd.concat([s.iloc[:10], pd.Series({'Otros': s.iloc[10:].sum()})]); return s, s.sum()
fig, axs = plt.subplots(2, 1, figsize=(16.6/2.54, 11.0/2.54), dpi=600)
for ax, col, let, tit in zip(axs, ['Genero', 'Especie'], 'AB', ['Género (o nivel superior)', 'Especie']):
    s, tot = datos(col)
    ax.pie(s.values, colors=COL, startangle=90, counterclock=False, center=(0, 0), radius=1,
           wedgeprops={'linewidth': 0.4, 'edgecolor': 'black'})
    ax.set_xlim(-1.25, 4.15); ax.set_ylim(-1.04, 1.04); ax.set_aspect('equal'); ax.axis('off')
    hs = [Patch(facecolor=c, edgecolor='black', linewidth=0.4) for c in COL]
    labs = ['%s  %s %%' % (etiqueta(k), ('%.1f' % (100*v/tot)).replace('.', ',')) for k, v in s.items()]
    ax.legend(hs, labs, loc='center left', bbox_to_anchor=(0.43, 0.5), frameon=False, title=tit,
              title_fontproperties={'weight': 'bold', 'size': 9}, handlelength=1.1, handleheight=1.0,
              labelspacing=0.22, alignment='left')
    ax.text(-1.25, 1.04, let, fontsize=11, fontweight='bold', va='top', ha='left')
plt.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.02, hspace=0.06)
out = os.path.join(B, '..', 'resultados', 'Fig2_composicion.png')
fig.savefig(out, dpi=600, metadata={'Software': None})
print(out)
