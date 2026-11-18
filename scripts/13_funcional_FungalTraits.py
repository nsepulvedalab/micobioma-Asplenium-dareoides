"""13_funcional_FungalTraits.py
Anotación funcional por género con la base de rasgos FungalTraits (Põlme et al., 2020) y Figura 3 de la tesis.
A cada ASV con género asignado se le asignan el estilo de vida primario, el estilo de vida secundario, la capacidad
de interacción endofítica y la capacidad patogénica para plantas de su género. El modo trófico se define agrupando
los estilos de vida primarios en saprótrofo, patótrofo y simbiótrofo. Las proporciones se ponderan por lecturas.
Entradas: datos/ASV_table.csv, datos/ASV_taxonomy_UNITE97_60.csv y la Tabla S1 «Traits of genera» de
Põlme et al. (2020) (material suplementario del artículo; no se redistribuye: descargarla y dejarla en
datos_externos/13225_2020_466_MOESM4_ESM.xlsx).
Salida: resultados/ (tablas CSV y Figura 3 en PNG)."""
import os, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get('DATOS') or os.path.abspath(os.path.join(B, '..', 'datos'))
FT = os.environ.get('FUNGALTRAITS') or os.path.abspath(os.path.join(B, '..', 'datos_externos', '13225_2020_466_MOESM4_ESM.xlsx'))
OUT = os.environ.get('SALIDA') or os.path.join(B, '..', 'resultados'); os.makedirs(OUT, exist_ok=True)

# ---- lecturas por ASV y género (mismo criterio que el script 12: texto antes del primer "_")
asv = pd.read_csv(os.path.join(D, 'ASV_table.csv')); asv = asv.rename(columns={asv.columns[0]: 'Sample'}).set_index('Sample')
lect = asv.sum(axis=0).rename('lecturas')
tax = pd.read_csv(os.path.join(D, 'ASV_taxonomy_UNITE97_60.csv'))
if 'ASV' not in tax.columns: tax = tax.rename(columns={tax.columns[0]: 'ASV'})
fp = tax.taxon.where(~(tax.taxon.isna() | (tax.taxon == 'Unassigned'))).str.split('|').str[0]
tax['Genero'] = fp.str.split('_').str[0]
a = tax[['ASV', 'Genero']].merge(lect, left_on='ASV', right_index=True)

# ---- base FungalTraits: una fila por género (la tabla trae 6 géneros duplicados; se conserva la primera fila)
ft = pd.read_excel(FT, sheet_name='data')
ft['g'] = ft.GENUS.astype(str).str.strip().str.lower()
ft = ft.drop_duplicates('g')
CAMPOS = ['primary_lifestyle', 'Secondary_lifestyle', 'Endophytic_interaction_capability_template',
          'Plant_pathogenic_capacity_template']
a['g'] = a.Genero.str.strip().str.lower()
m = a.merge(ft[['g'] + CAMPOS], on='g', how='left').drop(columns='g')
TROF = {'litter_saprotroph': 'Saprótrofo', 'wood_saprotroph': 'Saprótrofo', 'soil_saprotroph': 'Saprótrofo',
        'unspecified_saprotroph': 'Saprótrofo', 'dung_saprotroph': 'Saprótrofo', 'nectar/tap_saprotroph': 'Saprótrofo',
        'sooty_mold': 'Saprótrofo', 'plant_pathogen': 'Patótrofo', 'lichen_parasite': 'Patótrofo',
        'mycoparasite': 'Patótrofo', 'animal_parasite': 'Patótrofo', 'protistan_parasite': 'Patótrofo',
        'algal_parasite': 'Patótrofo', 'foliar_endophyte': 'Simbiótrofo', 'root_endophyte': 'Simbiótrofo',
        'ectomycorrhizal': 'Simbiótrofo', 'arbuscular_mycorrhizal': 'Simbiótrofo', 'ericoid_mycorrhizal': 'Simbiótrofo',
        'lichenized': 'Simbiótrofo', 'epiphyte': 'Simbiótrofo'}
anot = m.primary_lifestyle.notna()
falta = set(m.loc[anot, 'primary_lifestyle']) - set(TROF)
assert not falta, f'Estilos de vida sin modo trófico asignado: {falta}'
m['modo_trofico'] = m.primary_lifestyle.map(TROF)
m.to_csv(os.path.join(OUT, 'ASV_FungalTraits.csv'), index=False)

# ---- tabla por género
g = (m[m.Genero.notna()].groupby('Genero')
     .agg(lecturas=('lecturas', 'sum'), **{c: (c, 'first') for c in CAMPOS + ['modo_trofico']})
     .sort_values('lecturas', ascending=False))
T = m.lecturas.sum()
g.insert(1, 'pct_lecturas_total', (g.lecturas / T * 100).round(2))
g.to_csv(os.path.join(OUT, 'Tabla_genero_FungalTraits.csv'))

# ---- resumen (porcentajes sobre las lecturas anotadas)
an = m[anot]; TA = an.lecturas.sum()
ENDO = {'foliar_endophyte': 'Endófito foliar', 'root-associated': 'Asociado a raíz',
        'no_endophytic_capacity': 'Sin capacidad endofítica', 'root_endophyte_dark_septate': 'Endófito septado oscuro',
        'root_endophyte': 'Endófito de raíz'}
VIDA = {'litter_saprotroph': 'Saprótrofo de hojarasca', 'plant_pathogen': 'Patógeno de plantas',
        'lichen_parasite': 'Parásito de líquenes', 'foliar_endophyte': 'Endófito foliar', 'mycoparasite': 'Micoparásito'}
OTROS = 'Otros (cada uno < 0,5 %)'
A = (an.groupby('modo_trofico').lecturas.sum() / TA * 100).sort_values(ascending=False)
Bn = (an.groupby(an.Endophytic_interaction_capability_template.map(ENDO).fillna('Sin dato')).lecturas.sum() / TA * 100).sort_values(ascending=False)
C = an.groupby(an.primary_lifestyle.map(VIDA).fillna(OTROS)).lecturas.sum() / TA * 100
C = pd.concat([C.drop(OTROS).sort_values(ascending=False), C[[OTROS]]])
filas = [('Cobertura', 'Lecturas totales (sin quimeras)', T, 100.0),
         ('Cobertura', 'Asignadas a algún taxón (UNITE)', m[m.Genero.notna()].lecturas.sum(), None),
         ('Cobertura', 'Anotadas con FungalTraits', TA, None)]
res = pd.DataFrame(filas, columns=['panel', 'categoria', 'lecturas', 'pct'])
res['pct'] = (res.lecturas / T * 100).round(1)
for panel, s in (('A. Modo trófico', A), ('B. Capacidad endofítica', Bn), ('C. Estilo de vida primario', C)):
    for k, v in s.items():
        res.loc[len(res)] = [panel, k, round(v * TA / 100), round(v, 2)]
res.to_csv(os.path.join(OUT, 'Resumen_FungalTraits.csv'), index=False)
print(res.to_string(index=False))

# ---- Figura 3 (mismo formato que la figura de la tesis: 4251 × 6141 px, 600 dpi)
plt.rcParams['font.family'] = ['Carlito', 'Liberation Sans', 'DejaVu Sans']
pal = ['#A51C48', '#EA820C', '#FFE08A', '#B7E4A8', '#7FC8C8', '#9AA9D6']
fig, axs = plt.subplots(3, 1, figsize=(4251/600, 6141/600), dpi=600, gridspec_kw={'height_ratios': [3, 6, 6]})
for ax, s, t in zip(axs, (A, Bn, C), ('Modo trófico', 'Capacidad endofítica', 'Estilo de vida primario')):
    y = range(len(s))[::-1]
    ax.barh(list(y), s.values, color=pal[:len(s)], edgecolor='black', linewidth=0.6, height=0.75)
    ax.set_yticks(list(y)); ax.set_yticklabels(s.index, fontsize=11, color='#444444')
    for yy, v in zip(y, s.values):
        ax.text(v + 1, yy, f'{v:.1f} %'.replace('.', ','), va='center', fontsize=10.5, fontweight='bold')
    ax.set_xlim(0, 100); ax.set_xticks([0, 30, 60, 90]); ax.tick_params(axis='x', labelsize=10, colors='#444444', length=0)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#E5E5E5', linewidth=1); ax.set_axisbelow(True)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_xlabel('% de lecturas (asignadas)', fontsize=10.5)
    ax.set_title(t, fontsize=14, fontweight='bold', loc='left')
plt.subplots_adjust(left=0.32, right=0.93, top=0.945, bottom=0.05, hspace=0.55)
for ax, letra in zip(axs, 'ABC'):
    bb = ax.get_position(); fig.text(0.03, bb.y1 + 0.03, letra, fontsize=16)
fig.savefig(os.path.join(OUT, 'Figura3_FungalTraits.png'), dpi=600, metadata={'Software': None})
