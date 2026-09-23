"""Conditional normalized classical-ROM bounds; all sums in log2 space."""
import math

SETS = [
    dict(key='L1', n=880, q=18, beta=13, rows=1, mu=132, salt=256,
         target=128, design=128, alpha=200, D=3.24e-5, sigma=2.8, frodo=640,
         fmu=128, fb=134.5, fd=-138.7),
    dict(key='L2', n=1136, q=17, beta=13, rows=2, mu=216, salt=384,
         target=192, design=128, alpha=500, D=1.4e-5, sigma=2.3, frodo=976,
         fmu=192, fb=195.8, fd=-199.6),
    dict(key='L3', n=1640, q=17, beta=12, rows=2, mu=264, salt=512,
         target=256, design=64, alpha=1000, D=2.64e-5, sigma=1.4, frodo=1344,
         fmu=256, fb=250.8, fd=-252.5),
]
NEG = float('-inf')
LOG_T, LOG_Q = 32, 14  # t=2^32; qG=qH=t/2^18

def add(values):
    values = list(values)
    m = max(values)
    return m if m == NEG else m + math.log2(sum(2**(v-m) for v in values))

def bound(s, N, b=None, dp=None, dq=None, frodo=False):
    if N < 1:
        raise ValueError('N must be positive')
    ln = math.log2(N)
    mu = s['fmu'] if frodo else s['mu']
    n, nb, mb = (s['frodo'], 8, 8) if frodo else (s['n'], 24, s['rows'])
    samples = 2*n*nb + N*(2*mb*n+mb*nb)
    loss = samples*s['D']/math.log(2)
    p = 1-1/s['alpha']
    collision = math.log2(N*(N-1))-mu-s['salt'] if N>1 else NEG
    if frodo:
        # Frodo Theorem 1, Eq.(3), evaluated directly at t=2^32.
        inner = dict(six_q_over_M=math.log2(6)+LOG_Q-mu,
                     statistical=-512, N_over_M=ln-mu,
                     q_delta=LOG_Q+s['fd'], two_cpa=1+ln+LOG_T-s['fb'])
        outer = dict(main=p*(add(inner.values())+loss), dP=s['fd'],
                     qH_over_M=LOG_Q-mu, message_salt_collision=collision)
    else:
        # GHS OW conversion: 4*Adv_unnorm_CPA = 2*Adv_norm_CPA.
        inner = dict(qG_dQ=LOG_Q+dq, N_over_M=ln-mu,
                     two_T_qG_over_M=1+LOG_Q-mu, two_cpa=1+ln+LOG_T-b)
        outer = dict(main=p*(add(inner.values())+loss), dP=dp,
                     qH_over_M=LOG_Q-mu, message_salt_collision=collision,
                     salt_bad=math.log2(N*(N-1)/2)-s['salt'] if N>1 else NEG)
    log_bound = min(0.0, 1+add(outer.values()))  # normalized advantage <= 1
    clean = lambda d: {k: (None if v==NEG else v) for k,v in d.items()}
    return dict(N=N, security_bits=LOG_T-log_bound, log2_advantage_bound=log_bound,
                samples=samples, renyi_loss_bits=loss,
                inner_log2_terms=clean(inner), outer_log2_terms=clean(outer))
