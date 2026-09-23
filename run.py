"""Usage: python run.py   (paths are relative to this script)."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from security import SETS, bound

ROOT = Path(__file__).resolve().parent

def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')

def inputs(root, s):
    objs, hashes = {}, {}
    for kind in ('indcpa','dp','dq'):
        path = root/f"n{s['n']}_q{s['q']}_{kind}.json"
        objs[kind] = json.loads(path.read_text(encoding='utf-8'))
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        p = objs[kind]['parameters']
        for k in ('n','rows'):
            if p[k] != s[k]:
                raise ValueError(f'{path.name}: wrong {k}')
        if kind!='indcpa' and p['beta'] != 2**s['beta']:
            raise ValueError(f'{path.name}: wrong beta')
    c = objs['indcpa']
    if c['parameters']['q_log2']!=s['q'] or c['message_bits']!=s['mu']:
        raise ValueError('Wrong q or message length')
    b = c['adjusted_indcpa']['clsf']
    expected = c['security']['ciphertext_side']['models']['clsf']['raw_minimum']-math.log2(24+s['rows'])
    if not math.isclose(b, expected, abs_tol=1e-9):
        raise ValueError('Unexpected CPA hybrid loss')
    # One outward binary64 step when importing the stored interval endpoint.
    dp = math.nextafter(objs['dp']['log2_dP_upper']['upper_float'], math.inf)
    dq = math.nextafter(objs['dq']['log2_dQ_upper']['upper_float'], math.inf)
    return b, dp, dq, hashes

def plot(s, rows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9,5.6), layout='constrained')
    ns = [r['N'] for r in rows]
    marks = [1,2,4,8,16,32,64,128] + ([3] if s['key']=='L1' else [])
    indices = [n-1 for n in sorted(marks)]
    ax.plot(ns,[r['leech_bits'] for r in rows], marker='o',markevery=indices,
            label=f"LeechKEM-{s['key'][1]} (n={s['n']}, q=2^{s['q']})")
    ax.plot(ns,[r['frodo_bits'] for r in rows], marker='s',markevery=indices,
            label=f"FrodoKEM-{s['frodo']} (C-model inputs)")
    ax.axhline(s['target'],color='tab:red',ls='--',label=f"Leech design target: {s['target']} bit")
    ax.axvline(s['design'],color='gray',ls=':',label=f"Leech design N={s['design']}")
    y = rows[s['design']-1]['leech_bits']
    ax.annotate(f"{y:.3f} bit at N={s['design']}",xy=(s['design'],y),
                xytext=(.48,.87),textcoords='axes fraction',
                bbox=dict(facecolor='white',edgecolor='.8',alpha=.95),
                arrowprops=dict(arrowstyle='->'))
    ax.set_xscale('log',base=2)
    ax.set_xticks([1,2,4,8,16,32,64,128],labels=['1','2','4','8','16','32','64','128'])
    ax.set_xlabel('Challenges N under one public key')
    ax.set_ylabel('Conditional classical-ROM bound (bits)')
    ax.set_title('Core-SVP/C-LSF: normalized IND-CCA reduction estimates')
    ax.grid(alpha=.25)
    ax.legend(loc='lower left',fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(out/f"{s['key']}_comparison.{ext}",dpi=220)
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,default=ROOT/'data/raw')
    ap.add_argument('--output',type=Path,default=ROOT/'results')
    ap.add_argument('--no-plots',action='store_true')
    a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    summary=[]
    for s in SETS:
        b,dp,dq,hashes=inputs(a.inputs,s)
        leech=[bound(s,n,b,dp,dq) for n in range(1,129)]
        frodo=[bound(s,n,frodo=True) for n in range(1,129)]
        rows=[dict(N=l['N'],leech_bits=l['security_bits'],frodo_bits=f['security_bits']) for l,f in zip(leech,frodo)]
        with (a.output/f"{s['key']}_comparison.csv").open('w',newline='',encoding='utf-8') as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        write_json(a.output/f"{s['key']}_audit.json",dict(parameters=s,input_sha256=hashes,
                   bCPA=b,log2_dP=dp,log2_dQ=dq,leech=leech,frodo=frodo))
        if not a.no_plots: plot(s,rows,a.output)
        v=leech[s['design']-1]['security_bits']
        summary.append(dict(set=s['key'],n=s['n'],N=s['design'],security_bits=v,margin_bits=v-s['target']))
        print(f"{s['key']}: n={s['n']}, N={s['design']}, {v:.9f} bits, margin={v-s['target']:.9f}")
    write_json(a.output/'summary.json',summary)
    print(f'Results: {a.output.resolve()}')

if __name__=='__main__': main()
