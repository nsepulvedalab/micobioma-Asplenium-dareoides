"""12_figuras_4a8.py
Hacer las Figuras 4 a 8 de la tesis: α-diversidad por individuo, curvas de rarefacción, ocupancia,
dendrograma de Bray-Curtis y mapa de calor de los 30 ASV más abundantes.
Se repiten los mismos cálculos de los scripts 08 y 09 (paquete vegan de R), por lo que los números son los mismos:
Shannon, Simpson (1 - D), Pielou, riqueza esperada por rarefacción, Bray-Curtis sobre abundancia relativa
y agrupamiento por enlace promedio.
Entradas: datos/ASV_table.csv y datos/ASV_taxonomy_UNITE97_60.csv. Salida: resultados/ (PNG)."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.special import gammaln
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform, pdist
B = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get('DATOS') or os.path.abspath(os.path.join(B, '..', 'datos'))
OUT = os.environ.get('SALIDA') or os.path.join(B, '..', 'resultados'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 9, 'axes.titlesize': 9.5, 'axes.titleweight': 'bold',
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
                     'axes.spines.top': False, 'axes.spines.right': False})
coma = lambda v, d=1: (f'{v:.{d}f}').replace('.', ',')
miles = FuncFormatter(lambda v, p: f'{int(v):,}'.replace(',', '.'))
asv = pd.read_csv(os.path.join(D, 'ASV_table.csv')); asv = asv.rename(columns={asv.columns[0]: 'Sample'}).set_index('Sample')
asv.index = [s.replace('sample_', '') for s in asv.index]
mat = asv.values.astype(float); muestras = list(asv.index)
tax = pd.read_csv(os.path.join(D, 'ASV_taxonomy_UNITE97_60.csv'))
if 'ASV' not in tax.columns: tax = tax.rename(columns={tax.columns[0]: 'ASV'})
fp = tax.taxon.where(~(tax.taxon.isna() | (tax.taxon == 'Unassigned'))).str.split('|').str[0]
genero = dict(zip(tax.ASV, fp.str.split('_').str[0].fillna('Unassigned')))
COL = {'ADTN1': '#a71b4b', 'ADTN2': '#de5925', 'ADTN3': '#0090b5', 'ADTN4': '#584b9f'}
def guardar(fig, nombre): fig.savefig(os.path.join(OUT, nombre), dpi=600, metadata={'Software': None}); plt.close(fig)
# ---- alfa
S = (mat > 0).sum(1); p = mat / mat.sum(1, keepdims=True)
H = -(np.where(p > 0, p*np.log(np.where(p > 0, p, 1)), 0)).sum(1); Si = 1 - (p**2).sum(1); J = H/np.log(S)
ind = {'Riqueza observada (ASV)': S, 'Índice de Shannon': H, 'Índice de Simpson': Si, 'Equidad de Pielou': J}
fig, axs = plt.subplots(1, 4, figsize=(16.6/2.54, 6.0/2.54))
dec = [0, 2, 2, 2]
for ax, (t, v), dd in zip(axs, ind.items(), dec):
    for k, m in enumerate(muestras):
        ax.scatter(0.0 + (k-1.5)*0.07, v[k], s=22, color=COL[m], zorder=3, edgecolor='black', linewidth=0.3)
    mu, sd = v.mean(), v.std(ddof=1)
    ax.errorbar(0.35, mu, yerr=sd, fmt='o', color='black', ms=4, capsize=3, lw=1)
    ax.set_title(t); ax.set_xlim(-0.35, 0.6); ax.set_xticks([])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p, dd=dd: coma(y, dd)))
    ax.text(0.35, mu + sd, '\n', ha='center')
hs = [plt.Line2D([], [], marker='o', ls='', color=COL[m], markeredgecolor='black', markeredgewidth=0.3, label=m) for m in muestras]
hs.append(plt.Line2D([], [], marker='o', ls='-', color='black', ms=4, label='Media ± DE'))
fig.legend(handles=hs, loc='lower center', ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.17, wspace=0.55)
guardar(fig, 'Fig4_alfa.png')
# ---- rarefacción (riqueza esperada; como vegan::rarefy)
def rarefy(x, n):
    x = x[x > 0]; N = x.sum()
    if n >= N: return float(len(x))
    lc = lambda a, b: gammaln(a+1) - gammaln(b+1) - gammaln(a-b+1)
    v = np.where(N - x >= n, np.exp(lc(N - x, n) - lc(N, n)), 0.0)
    return float((1 - v).sum())
mind = int(mat.sum(1).min())
fig, ax = plt.subplots(figsize=(16.6/2.54, 8.0/2.54))
for k, m in enumerate(muestras):
    N = int(mat[k].sum()); ns = list(range(1, N, 500)) + [N]
    ax.plot(ns, [rarefy(mat[k], n) for n in ns], color=COL[m], lw=1.6, label=f'{m} ({S[k]} ASV)')
ax.axvline(mind, color='grey', ls='--', lw=1)
ax.text(mind, 148, f'profundidad mínima: {mind:,} lecturas '.replace(',', '.'), color='grey', fontsize=8, va='top', ha='right')
ax.xaxis.set_major_formatter(miles); ax.set_xlabel('Número de lecturas'); ax.set_ylabel('Riqueza esperada (ASV)')
ax.set_ylim(0, 150); ax.legend(frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.03)); fig.tight_layout(); guardar(fig, 'Fig5_rarefaccion.png')
# ---- ocupancia
occ = (mat > 0).sum(0); tot = mat.sum()
fig, ax = plt.subplots(figsize=(12/2.54, 7.5/2.54))
for n in range(1, 5):
    sel = occ == n; na = int(sel.sum()); pct = 100*mat[:, sel].sum()/tot
    ax.bar(n, na, color='#a71b4b' if n == 4 else '#cccccc', edgecolor='black', lw=0.5, width=0.7)
    ax.text(n, na + 3, f'{na} ASV\n{coma(pct)} %', ha='center', va='bottom', fontsize=8.5)
ax.set_xticks(range(1, 5)); ax.set_xlabel('Número de individuos en que aparece el ASV'); ax.set_ylabel('Número de ASV')
ax.set_ylim(0, 245); fig.tight_layout(); guardar(fig, 'Fig6_ocupancia.png')
# ---- dendrograma Bray-Curtis (abundancia relativa, enlace promedio)
bc = pdist(p, metric='braycurtis'); Z = linkage(bc, 'average')
fig, ax = plt.subplots(figsize=(12/2.54, 7.0/2.54))
dn = dendrogram(Z, labels=muestras, ax=ax, color_threshold=0, above_threshold_color='black', leaf_font_size=9)
ax.set_ylabel('Disimilitud de Bray-Curtis'); ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: coma(y, 1)))
ax.spines['bottom'].set_visible(False); ax.tick_params(axis='x', length=0); fig.tight_layout(); guardar(fig, 'Fig7_dendrograma.png')
orden = dn['ivl']
# ---- mapa de calor top 30
reads_asv = mat.sum(0); top = np.argsort(-reads_asv)[:30]
rel = p[:, top]                                     # filas: muestras
o = dendrogram(linkage(pdist(rel.T, 'euclidean'), 'average'), no_plot=True)['leaves']
top = top[o]; ids = asv.columns[top]
cols = [muestras.index(m) for m in orden]
M = p[np.ix_(cols, top)].T
SUP = ('ales', 'aceae', 'mycota', 'mycetes')
def lab(a):
    g = genero.get(a, 'Unassigned')
    if g == 'Unassigned': return f'{a}  No asignado'
    return f'{a}  {g}' if g.endswith(SUP) else f'{a}  ' + r'$\it{%s}$' % g
fig, ax = plt.subplots(figsize=(14/2.54, 17/2.54))
im = ax.imshow(M, aspect='auto', cmap='viridis', vmin=0)
ax.set_yticks(range(len(ids))); ax.set_yticklabels([lab(a) for a in ids], fontsize=7.5)
ax.set_xticks(range(len(orden))); ax.set_xticklabels(orden)
ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.set_label('Abundancia relativa')
from matplotlib.ticker import MultipleLocator; cb.ax.yaxis.set_major_locator(MultipleLocator(0.1)); cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: coma(y, 1)))
fig.tight_layout(); guardar(fig, 'Fig8_mapa_calor.png')
print('BC', np.round(squareform(bc), 3).tolist(), 'S', S.tolist(), 'H', np.round(H, 3).tolist(), 'Si', np.round(Si, 4).tolist(), 'J', np.round(J, 3).tolist(), 'orden', orden)
