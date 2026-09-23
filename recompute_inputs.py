"""Optional fixed-parameter input generation. Never searches over n."""
import argparse
import json
import sys
from pathlib import Path
from security import SETS
from run import write_json

ROOT=Path(__file__).resolve().parent
BUNDLED_ESTIMATES=ROOT/'third_party/frodo_estimates'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['dp','dq','lattice','all'],required=True)
    p.add_argument('--set',choices=['L1','L2','L3','all'],default='all')
    p.add_argument(
        '--estimates-dir',
        type=Path,
        default=BUNDLED_ESTIMATES,
        help=(
            'Frodo estimates directory containing cost_models.py and '
            'lattice_estimator (default: bundled third_party copy)'
        ),
    )
    p.add_argument('--output',type=Path,default=ROOT/'data/regenerated')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if a.stage in ('lattice','all'):
        estimates_dir=a.estimates_dir.resolve()
        required=(estimates_dir/'cost_models.py',estimates_dir/'lattice_estimator')
        missing=[str(path) for path in required if not path.exists()]
        if missing:
            p.error('missing estimator dependency: '+', '.join(missing))
        try:
            import sage.all  # noqa: F401
        except ModuleNotFoundError:
            p.error(
                'the lattice stage requires a SageMath-enabled Python; '
                'activate the Sage environment and run this command again'
            )
        sys.path.insert(0,str(estimates_dir))
    for s in SETS:
        if a.set not in ('all',s['key']): continue
        prefix=f"n{s['n']}_q{s['q']}"
        if a.stage in ('dp','all'):
            from generators import leech_dp_chernoff as dp
            spec=dp.Preset(s['key'],s['key'],s['n'],2**s['beta'],f"frodo{s['frodo']}",s['rows'])
            write_json(a.output/f'{prefix}_dp.json',dp.compute_preset(spec,90))
        if a.stage in ('dq','all'):
            from generators import leech_dq_chernoff as dq
            spec=dq.Preset(s['key'],s['n'],2**s['beta'],s['sigma'],s['rows'])
            write_json(a.output/f'{prefix}_dq.json',dq.compute_preset(spec,90,120,160))
        if a.stage in ('lattice','all'):
            import math
            from generators.lattice_backend import estimate_one_side
            side=estimate_one_side(s['key'],s['n'],2**s['q'],s['n']+24,s['sigma'])
            result=dict(parameters=dict(n=s['n'],q_log2=s['q'],beta_log2=s['beta'],rows=s['rows']),
                        message_bits=s['mu'],hybrid_loss_bits=math.log2(24+s['rows']),
                        adjusted_indcpa={k:v['raw_minimum']-math.log2(24+s['rows']) for k,v in side['models'].items()},
                        security=dict(ciphertext_side=side))
            write_json(a.output/f'{prefix}_indcpa.json',result)
    print('Generated inputs:',a.output.resolve())

if __name__=='__main__': main()
