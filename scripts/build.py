#!/usr/bin/env python3
"""TTW CFR public access defibrillator gap map - monthly build.

Reads the newest defibrillator register extract in data/defibrillators/ (.csv or .xlsx, The Circuit
column names), models coverage against Census 2021 population and writes the web map to docs/.

    python scripts/build.py                       # build docs/ (GitHub Pages site)
    python scripts/build.py --standalone out.html # also write a single self-contained HTML file
"""
import argparse, base64, csv, datetime as dt, glob, io, json, math, os, re, subprocess, sys
import numpy as np, pandas as pd
from PIL import Image
from pyproj import Transformer
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from shapely import contains_xy
from shapely.geometry import shape, Point
from shapely.prepared import prep

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATIC = os.path.join(ROOT, 'data', 'static')
INBOX = os.path.join(ROOT, 'data', 'defibrillators')
DOCS = os.path.join(ROOT, 'docs')
OUTD = os.path.join(DOCS, 'data')

# ---------------- model settings ----------------
CELL = 100.0              # grid size, metres
DECAY_M = 100.0           # retrieval probability halves every 100 m walked (Deakin et al., Heart 2018)
DETOUR = 1.25             # straight line -> walking route
AVAIL = {'public247': 1.0, 'publicvar': 0.35, 'restricted': 0.0}
AGE_WEIGHT = 3.5          # residents 65+ vs under 65 (cardiac arrest incidence)
DEP_MAX = 0.5             # up to +50% for the most deprived areas (IMD 2025 rank)
BUFFER_M = 3000           # keep AEDs this far outside the area (they cover border streets)
N_SITES = 20
N_UPGRADES = 12
KERNEL_R = 6              # cells (600 m) considered for a new AED's reach

T = Transformer.from_crs(4326, 27700, always_xy=True)
Ti = Transformer.from_crs(27700, 4326, always_xy=True)
def decay(d): return np.power(0.5, (np.asarray(d) * DETOUR) / DECAY_M)


def newest_input():
    files = [f for f in glob.glob(os.path.join(INBOX, '*')) if f.lower().endswith(('.csv', '.xlsx', '.xls', '.zip'))]
    if not files:
        sys.exit('No .csv, .xlsx or .zip file found in data/defibrillators/')
    def when(f):
        try:
            t = subprocess.run(['git', 'log', '-1', '--format=%ct', '--', f], cwd=ROOT, capture_output=True, text=True).stdout.strip()
            if t: return int(t)
        except Exception:
            pass
        return int(os.path.getmtime(f))
    return max(files, key=lambda f: (when(f), f))


def read_register(path):
    if path.lower().endswith('.zip'):
        import zipfile, tempfile
        z = zipfile.ZipFile(path)
        inner = [n for n in z.namelist() if n.lower().endswith(('.csv', '.xlsx', '.xls')) and not n.startswith('__MACOSX')]
        if not inner: sys.exit('The zip file has no .csv or .xlsx inside')
        tmp = os.path.join(tempfile.mkdtemp(), os.path.basename(inner[0])); open(tmp, 'wb').write(z.read(inner[0]))
        return read_register(tmp)
    if path.lower().endswith('.csv'):
        df = pd.read_csv(path, encoding_errors='replace', low_memory=False)
    else:
        sheets = pd.read_excel(path, sheet_name=None)
        df = next((s for s in sheets.values() if any(str(c).strip().lower() in ('lat', 'latitude') for c in s.columns)), None)
        if df is None: sys.exit('No sheet with a lat/latitude column in ' + path)
    df.columns = [str(c).strip().lower().replace(' ', '_') for c in df.columns]
    ren = {'latitude': 'lat', 'longitude': 'long', 'lng': 'long', 'lon': 'long'}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    need = ['lat', 'long', 'defibrillators_access_type', 'defibrillators_availability']
    miss = [c for c in need if c not in df.columns]
    if miss: sys.exit('Missing columns in register file: ' + ', '.join(miss))
    for c in ['unique_identifier', 'location_name', 'address_line1', 'address_line2', 'address_line3', 'address_city', 'address_post_code']:
        if c not in df.columns: df[c] = ''
    df['lat'] = pd.to_numeric(df.lat, errors='coerce'); df['long'] = pd.to_numeric(df.long, errors='coerce')
    return df.dropna(subset=['lat', 'long'])


def classify(r):
    acc = str(r.defibrillators_access_type).strip().lower()
    av = str(r.defibrillators_availability).strip().lower()
    if acc.startswith('public') and ('24' in av): return 'public247'
    if acc.startswith('public'): return 'publicvar'
    return 'restricted'


def clean(s):
    s = '' if s is None or (isinstance(s, float) and math.isnan(s)) else str(s).strip()
    return '' if s.lower() == 'nan' else s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--standalone'); ap.add_argument('--input')
    ap.add_argument('--prune', action='store_true', help='delete older register files once the build succeeds')
    args = ap.parse_args()
    src = args.input or newest_input()
    print('Register file:', os.path.relpath(src, ROOT))
    os.makedirs(OUTD, exist_ok=True)

    area_ll = shape(json.load(open(os.path.join(STATIC, 'area.geojson')))['geometry'])
    dist_fc = json.load(open(os.path.join(STATIC, 'postcode_districts.geojson')))
    import shapely
    to_bng = lambda g: shapely.transform(g, lambda c: np.c_[T.transform(c[:, 0], c[:, 1])])
    area = to_bng(area_ll)
    districts = {f['properties']['district']: to_bng(shape(f['geometry'])) for f in dist_fc['features']}
    P = pd.read_csv(os.path.join(STATIC, 'population_by_postcode.csv'))
    places = json.load(open(os.path.join(STATIC, 'places.json')))

    # ---------------- defibrillators ----------------
    D = read_register(src)
    w, s, e, n = area_ll.bounds
    D = D[D.lat.between(s - 0.06, n + 0.06) & D.long.between(w - 0.09, e + 0.09)].copy()
    D['x'], D['y'] = T.transform(D.long.values, D.lat.values)
    near = prep(area.buffer(BUFFER_M)); inside = prep(area)
    D = D[[near.contains(Point(a, b)) for a, b in zip(D.x, D.y)]].copy()
    D['in_area'] = [inside.contains(Point(a, b)) for a, b in zip(D.x, D.y)]
    D['cls'] = D.apply(classify, axis=1)
    D['w'] = D.cls.map(AVAIL)
    D['district'] = ''
    for d, g in districts.items():
        pg = prep(g)
        D.loc[[pg.contains(Point(a, b)) for a, b in zip(D.x, D.y)], 'district'] = d
    print('AEDs in area:', D.in_area.sum(), D[D.in_area].cls.value_counts().to_dict())

    # ---------------- population grid ----------------
    P['dep'] = 1.0 + DEP_MAX * (1 - P.imd_rank / 33755.0)
    P['burden'] = ((P.residents - P.aged65) + AGE_WEIGHT * P.aged65) * P.dep
    minx, miny, maxx, maxy = area.bounds
    x0 = math.floor(minx / CELL) * CELL - 1000; y0 = math.floor(miny / CELL) * CELL - 1000
    nx = int((maxx - x0 + 1000) / CELL) + 1; ny = int((maxy - y0 + 1000) / CELL) + 1
    cx = x0 + (np.arange(nx) + 0.5) * CELL; cy = y0 + (np.arange(ny) + 0.5) * CELL
    XX, YY = np.meshgrid(cx, cy)
    mask = contains_xy(area, XX, YY)
    iy = ((P.northing - y0) // CELL).astype(int).values; ix = ((P.easting - x0) // CELL).astype(int).values
    popg = np.zeros((ny, nx)); burg = np.zeros((ny, nx)); p65g = np.zeros((ny, nx))
    np.add.at(popg, (iy, ix), P.residents); np.add.at(burg, (iy, ix), P.burden); np.add.at(p65g, (iy, ix), P.aged65)
    burg *= popg.sum() / burg.sum()

    R = KERNEL_R
    def window(i, j): return (slice(max(i - R, 0), i + R + 1), slice(max(j - R, 0), j + R + 1))
    def cell_of(x, y): return int((y - y0) // CELL), int((x - x0) // CELL)
    cov = np.zeros((ny, nx))
    for _, r in D[D.w > 0].iterrows():
        sl = window(*cell_of(r.x, r.y))
        cov[sl] = np.maximum(cov[sl], r.w * decay(np.hypot(XX[sl] - r.x, YY[sl] - r.y)))

    # ---------------- headline coverage (straight-line distance) ----------------
    tot = float(P.residents.sum())
    def dist_to(pts):
        if len(pts) == 0: return np.full(len(P), 1e9)
        return cKDTree(pts).query(P[['easting', 'northing']].values)[0]
    d247 = dist_to(D[D.cls == 'public247'][['x', 'y']].values)
    dall = dist_to(D[['x', 'y']].values)
    stats = {str(t): {'p247': float(P.residents[d247 <= t].sum() / tot), 'pall': float(P.residents[dall <= t].sum() / tot)} for t in (100, 200, 500)}
    bydist = []
    for d in sorted(districts, key=lambda k: int(k[2:])):
        m = (P.district == d).values
        pp = float(P.residents[m].sum())
        dd = D[D.in_area & (D.district == d)]
        bydist.append({'district': d, 'residents': round(pp), 'aeds': int(len(dd)), 'public247': int((dd.cls == 'public247').sum()),
                       **{f'c{t}': float(P.residents[m & (d247 <= t)].sum() / pp) for t in (100, 200, 500)}})

    # ---------------- marginal benefit + recommendations ----------------
    offs = [(di, dj, float(decay(math.hypot(di, dj) * CELL))) for di in range(-R, R + 1) for dj in range(-R, R + 1) if math.hypot(di, dj) <= R]
    def shift(a, di, dj):
        out = np.zeros_like(a)
        out[max(-di, 0):ny + min(-di, 0), max(-dj, 0):nx + min(-dj, 0)] = a[max(di, 0):ny + min(di, 0), max(dj, 0):nx + min(dj, 0)]
        return out
    def benefit(c):
        b = np.zeros((ny, nx))
        for di, dj, wgt in offs: b += shift(burg * np.maximum(0, wgt - c), di, dj)
        b[~mask] = 0
        return b
    def name_of(r):
        nm = clean(r.location_name)
        return nm or clean(r.address_line1) or 'Defibrillator'

    covk = cov.copy(); ups = []
    cand = D[D.in_area & (D.cls != 'public247')].copy()
    cand['key'] = (cand.x // 30).astype(int).astype(str) + '_' + (cand.y // 30).astype(int).astype(str)
    cand = cand.sort_values('cls').drop_duplicates('key')
    for _ in range(N_UPGRADES):
        best = None
        for idx, r in cand.iterrows():
            if idx in [u['_idx'] for u in ups]: continue
            sl = window(*cell_of(r.x, r.y)); d2 = np.hypot(XX[sl] - r.x, YY[sl] - r.y)
            g = float((burg[sl] * np.maximum(0, decay(d2) - covk[sl])).sum())
            if best is None or g > best[0]: best = (g, idx, r, sl, d2)
        if best is None: break
        g, idx, r, sl, d2 = best
        ppl = float(popg[sl][(d2 <= 200) & (covk[sl] < decay(200))].sum())
        if ppl < 150: break
        covk[sl] = np.maximum(covk[sl], decay(d2))
        ups.append({'_idx': idx, 'name': name_of(r), 'address': ', '.join(v for v in [clean(r.address_line1), clean(r.address_city), clean(r.address_post_code)] if v),
                    'access': r.cls, 'residents_gained': int(round(ppl, -1)), 'x': r.x, 'y': r.y, 'district': r.district})
    B = benefit(covk); need_surface = B.copy()
    pc_tree = cKDTree(P[['easting', 'northing']].values)
    t247 = cKDTree(D[D.cls == 'public247'][['x', 'y']].values) if (D.cls == 'public247').any() else None
    sites = []
    for k in range(N_SITES):
        i, j = np.unravel_index(np.argmax(B), B.shape)
        sx, sy = cx[j], cy[i]
        sl = window(i, j); d2 = np.hypot(XX[sl] - sx, YY[sl] - sy)
        newly = (d2 <= 200) & (covk[sl] < decay(200))
        _, pi = pc_tree.query([sx, sy])
        sites.append({'rank': k + 1, 'postcode': P.postcode.iloc[pi], 'district': P.district.iloc[pi],
                      'residents_gained': int(round(float(popg[sl][newly].sum()), -1)), 'aged65_gained': int(round(float(p65g[sl][newly].sum()), -1)),
                      'nearest_247_m': int(round(t247.query([sx, sy])[0], -1)) if t247 else None, 'x': float(sx), 'y': float(sy)})
        covk[sl] = np.maximum(covk[sl], decay(d2)); B = benefit(covk)

    def share(extra):
        pts = np.r_[D[D.cls == 'public247'][['x', 'y']].values, np.array(extra).reshape(-1, 2)]
        return float(P.residents[dist_to(pts) <= 200].sum() / tot)
    upxy = [[u['x'], u['y']] for u in ups]; sxy = [[s_['x'], s_['y']] for s_ in sites]
    impact = {'now': stats['200']['p247'], 'ups': share(upxy), 'sites10': share(upxy + sxy[:10]), 'sites20': share(upxy + sxy[:20])}

    # ---------------- localities ----------------
    pl = [p for p in places if p[3] in ('town', 'village', 'suburb', 'hamlet', 'neighbourhood')]
    plx, ply = T.transform([p[0] for p in pl], [p[1] for p in pl]); ptree = cKDTree(np.c_[plx, ply])
    towns = [p for p in pl if p[3] == 'town']; tx_, ty_ = T.transform([p[0] for p in towns], [p[1] for p in towns]); ttree = cKDTree(np.c_[tx_, ty_])
    def locality(x, y):
        d, ii = ptree.query([x, y], k=6)
        name, typ = pl[ii[0]][2], pl[ii[0]][3]
        for dd, k_ in zip(d, ii):
            if pl[k_][3] in ('suburb', 'neighbourhood', 'village', 'town') and dd < 1500: name, typ = pl[k_][2], pl[k_][3]; break
        if typ in ('village', 'town'): return name
        dt_, it = ttree.query([x, y]); town = towns[it][2]
        return name if (town == name or name.endswith(town) or dt_ > 4000) else f'{name}, {town}'

    # ---------------- write outputs ----------------
    def feat(x, y, props):
        lon, lat = Ti.transform(x, y)
        return {'type': 'Feature', 'properties': props, 'geometry': {'type': 'Point', 'coordinates': [round(lon, 6), round(lat, 6)]}}
    lab = {'public247': 'Public, 24/7', 'publicvar': 'Public, limited hours', 'restricted': 'Restricted access'}
    aed_fc = {'type': 'FeatureCollection', 'features': [feat(r.x, r.y, {
        'id': clean(r.unique_identifier), 'name': name_of(r),
        'address': ', '.join(v for v in [clean(r.address_line1), clean(r.address_line2), clean(r.address_city), clean(r.address_post_code)] if v),
        'access': r.cls, 'access_label': lab[r.cls], 'availability': clean(r.defibrillators_availability), 'in_area': bool(r.in_area), 'district': r.district}) for _, r in D.iterrows()]}
    site_fc = {'type': 'FeatureCollection', 'features': [feat(s_['x'], s_['y'], {**{k: v for k, v in s_.items() if k not in ('x', 'y')}, 'locality': locality(s_['x'], s_['y'])}) for s_ in sites]}
    up_fc = {'type': 'FeatureCollection', 'features': [feat(u['x'], u['y'], {**{k: v for k, v in u.items() if k not in ('x', 'y', '_idx')}, 'rank': i + 1, 'access_label': lab[u['access']], 'locality': locality(u['x'], u['y'])}) for i, u in enumerate(ups)]}
    for name, fc in [('defibrillators', aed_fc), ('priority_sites', site_fc), ('quick_wins', up_fc)]:
        json.dump(fc, open(os.path.join(OUTD, name + '.geojson'), 'w'), separators=(',', ':'))
    # plain CSV copies (for spreadsheets, Google My Maps or ArcGIS Online)
    def to_csv(fc, name):
        rows = [{**{k: v for k, v in f_['properties'].items() if not isinstance(v, (dict, list))}, 'longitude': f_['geometry']['coordinates'][0], 'latitude': f_['geometry']['coordinates'][1]} for f_ in fc['features']]
        pd.DataFrame(rows).to_csv(os.path.join(OUTD, name + '.csv'), index=False)
    to_csv(aed_fc, 'defibrillators'); to_csv(site_fc, 'priority_sites'); to_csv(up_fc, 'quick_wins')
    for f in ('postcode_districts.geojson', 'output_areas.geojson', 'area.geojson'):
        open(os.path.join(OUTD, f), 'w').write(open(os.path.join(STATIC, f)).read())

    # rasters (Web Mercator aligned PNG overlays)
    lons, lats = Ti.transform([x0, x0 + nx * CELL, x0, x0 + nx * CELL], [y0, y0, y0 + ny * CELL, y0 + ny * CELL])
    W, S, E, N = min(lons), min(lats), max(lons), max(lats)
    my = lambda lat: math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    PX = 1600; PY = int(PX * (my(N) - my(S)) / math.radians(E - W))
    LON = W + (np.arange(PX) + 0.5) / PX * (E - W)
    LAT = np.degrees(2 * np.arctan(np.exp(my(N) - (np.arange(PY) + 0.5) / PY * (my(N) - my(S)))) - math.pi / 2)
    LL, AA = np.meshgrid(LON, LAT); BX, BY = T.transform(LL, AA)
    J = ((BX - x0) // CELL).astype(int); I = ((BY - y0) // CELL).astype(int)
    ok = (I >= 0) & (I < ny) & (J >= 0) & (J < nx); I = np.clip(I, 0, ny - 1); J = np.clip(J, 0, nx - 1)
    hx = lambda h, a: [int(h[1:3], 16), int(h[3:5], 16), int(h[5:7], 16), a]
    def ramp(stops):
        pos = np.array([q[0] for q in stops]); cols = np.array([q[1] for q in stops], float)
        return lambda t: np.stack([np.interp(np.clip(t, 0, 1), pos, cols[:, k]) for k in range(4)], -1)
    def png(grid, f, vmax, smooth=0.0, gamma=1.0, name=None):
        a = gaussian_filter(grid.astype(float), smooth) if smooth else grid.astype(float)
        a = np.where(mask, a, 0)
        img = f((a / vmax) ** gamma)[I, J]; img[~ok] = 0; img[~mask[I, J]] = 0
        b = io.BytesIO(); Image.fromarray(img.astype(np.uint8), 'RGBA').save(b, 'PNG', optimize=True)
        open(os.path.join(OUTD, name + '.png'), 'wb').write(b.getvalue())
        return 'data:image/png;base64,' + base64.b64encode(b.getvalue()).decode()
    heat = ramp([(0, hx('#C5D1E0', 0)), (0.10, hx('#C5D1E0', 0)), (0.22, hx('#8BA3C1', 170)), (0.45, hx('#51769F', 215)), (0.75, hx('#174882', 235)), (1, hx('#312783', 245))])
    q = lambda a: float(np.percentile(a[mask & (a > 0)], 99.6)) if (mask & (a > 0)).any() else 1.0
    rasters = {
        'need': png(need_surface, heat, q(gaussian_filter(need_surface, 1.0)), 1.0, 0.6, 'need'),
        'unmet': png(burg * (1 - cov), heat, q(gaussian_filter(burg * (1 - cov), 1.5)), 1.5, 0.8, 'unmet'),
        'coverage': png(cov, ramp([(0, hx('#83A090', 0)), (0.06, hx('#83A090', 50)), (0.25, hx('#83A090', 130)), (1, hx('#084121', 200))]), 1.0, name='coverage'),
    }
    dens = gaussian_filter(popg, 1.0) * 100
    vmax_d = q(dens)
    rasters['density'] = png(dens, ramp([(0, hx('#9893C1', 0)), (0.03, hx('#9893C1', 60)), (0.3, hx('#9893C1', 160)), (0.7, hx('#645DA2', 210)), (1, hx('#312783', 230))]), vmax_d, gamma=0.7, name='density')

    m = re.search(r'(20\d\d)[-_]?(\d\d)[-_]?(\d\d)', os.path.basename(src))
    extract_date = f'{m.group(1)}-{m.group(2)}-{m.group(3)}' if m else None
    counts = D[D.in_area].cls.value_counts().to_dict()
    summary = {'built': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'), 'source_file': os.path.basename(src), 'extract_date': extract_date,
               'residents': round(tot), 'aged65': round(float(P.aged65.sum())), 'counts': {k: int(counts.get(k, 0)) for k in AVAIL},
               'stats': stats, 'districts': bydist, 'impact': impact, 'bounds': [[S, W], [N, E]], 'density_max': vmax_d,
               'settings': {'decay_m': DECAY_M, 'detour': DETOUR, 'availability_weights': AVAIL, 'age_weight': AGE_WEIGHT, 'deprivation_max': DEP_MAX}}

    # history (one row per register extract)
    hpath = os.path.join(OUTD, 'history.csv')
    hist = pd.read_csv(hpath) if os.path.exists(hpath) else pd.DataFrame()
    row = {'extract': extract_date or summary['built'][:10], 'source_file': summary['source_file'], 'aeds': int(D.in_area.sum()),
           'public_247': summary['counts']['public247'], 'within_100m': round(stats['100']['p247'], 4), 'within_200m': round(stats['200']['p247'], 4), 'within_500m': round(stats['500']['p247'], 4)}
    if len(hist): hist = hist[hist.source_file != row['source_file']]
    hist = pd.concat([hist, pd.DataFrame([row])], ignore_index=True).sort_values('extract')
    hist.to_csv(hpath, index=False)
    summary['history'] = hist.to_dict('records')
    json.dump(summary, open(os.path.join(OUTD, 'summary.json'), 'w'), indent=1)

    # site page
    tpl = open(os.path.join(ROOT, 'scripts', 'map_template.html')).read()
    os.makedirs(os.path.join(DOCS, 'brand'), exist_ok=True)
    for f in ('ttw_symbol.png', 'favicon.png'):
        open(os.path.join(DOCS, 'brand', f), 'wb').write(open(os.path.join(ROOT, 'data', 'brand', f), 'rb').read())
    open(os.path.join(DOCS, 'index.html'), 'w').write(tpl.replace('/*EMBED*/', '').replace('__LOGO__', 'brand/ttw_symbol.png').replace('__FAVICON__', 'brand/favicon.png'))
    print(f"Residents {tot:,.0f}; within 200 m of a 24/7 public AED: {stats['200']['p247']:.1%}")

    if args.prune:
        for f in glob.glob(os.path.join(INBOX, '*')):
            if os.path.abspath(f) != os.path.abspath(src) and f.lower().endswith(('.csv', '.xlsx', '.xls', '.zip')):
                os.remove(f); print('Removed older register file:', os.path.basename(f))

    if args.standalone:
        emb = {'summary': summary, 'aeds': aed_fc, 'sites': site_fc, 'ups': up_fc,
               'districts': dist_fc, 'oas': json.load(open(os.path.join(STATIC, 'output_areas.geojson'))), 'rasters': rasters,
               'base': json.load(open(os.path.join(STATIC, 'basemap_vector.json')))}
        lcss = open(os.path.join(ROOT, 'scripts', 'leaflet.css')).read() if os.path.exists(os.path.join(ROOT, 'scripts', 'leaflet.css')) else ''
        b64 = lambda f: 'data:image/png;base64,' + base64.b64encode(open(os.path.join(ROOT, 'data', 'brand', f), 'rb').read()).decode()
        tpl = tpl.replace('__LOGO__', b64('ttw_symbol.png')).replace('__FAVICON__', b64('favicon.png'))
        html = tpl.replace('/*EMBED*/', 'window.EMBED=' + json.dumps(emb, separators=(',', ':')) + ';').replace('/*LEAFLET_CSS*/', lcss)
        # page body only (the Claude artifact host supplies the document shell)
        html = re.sub(r'(?is)<link rel="icon"[^>]*>|<link rel="stylesheet" href="https://unpkg[^>]*>|<!doctype html>|</?html[^>]*>|</?head>|</?body>|<meta charset="utf-8">|<meta name="viewport"[^>]*>', '', html)
        open(args.standalone, 'w').write(html)
        print('Standalone page:', args.standalone, f'{len(html) / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
