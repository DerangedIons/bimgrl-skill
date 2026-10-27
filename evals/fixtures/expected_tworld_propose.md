# Proposed blocks

Validation: valid (gotranx validate_blocks)

## Implicit blocks

- C0, C1, C2, I_, O_: Markov chain (exact); rule 1
- ryr_CaRI, ryr_I, ryr_O, ryr_R: Markov chain (closure); rule 1
- ryr_CaRI_p, ryr_I_p, ryr_O_p, ryr_R_p: Markov chain (closure); rule 1
- buffers_SLHsl, buffers_SLLsl, casl: pool with its buffers, stiff mode Re λ -163 /ms owned by casl; rule 2
- buffers_Csqn, casr, buffers_SLHj, buffers_SLLj, caj: pool with its buffers, stiff mode Re λ -259.6 /ms owned by buffers_Csqn, casr; rule 2 + pool with its buffers, stiff mode Re λ -141.8 /ms owned by buffers_SLLj, caj; rule 2; joined so that the release conserves Ca (rule 3)

## Sub-step blocks

- v, h, hL, hLp, h_P, hp, hp_P, j, j_P, jp, jp_P, m, mL, m_P (proportional, h_min 0.001, triggers [{'state': 'v', 'atol': 0.5}]): V with the INa/INaL gates, refined in proportion to |dV/dt| (rule 4)
- buffers_Csqn, casr, buffers_SLHj, buffers_SLLj, caj, ryr_CaRI, ryr_I, ryr_O, ryr_R, ryr_CaRI_p, ryr_I_p, ryr_O_p, ryr_R_p (binary, h_min 0.005, triggers [{'state': 'caj', 'rtol': 0.02}, {'state': 'ryr_O', 'atol': 0.001}, {'state': 'ryr_O_p', 'atol': 0.001}]): junctional Ca with its buffers, the release channel and the SR, refined while release is active (rule 5)

## Left explicit (by rule)

- chain contraction_TmBlocked, contraction_XS, contraction_XW does not feed V (contraction?): left explicit
- binding cluster buffers_CaM, buffers_Myosin_ca, buffers_SRB, buffers_TnCHc, cai owns no stiff mode (threshold -20 /ms): left explicit
- binding cluster buffers_NaBj, naj owns no stiff mode (threshold -20 /ms): left explicit
- binding cluster buffers_NaBsl, nasl owns no stiff mode (threshold -20 /ms): left explicit

## Judgment calls for the user

- the voltage-dependent chain C0, C1, C2, I_, O_ is an implicit block outside the V block; nest it (--nest-v-markov) if its states go negative at large steps (Winslow 1999 L-type)

## gotranx suggest-blocks at dt 0.1, for comparison

- implicit: ['C0', 'C1', 'C2', 'I_', 'O_']
- implicit: ['buffers_Csqn', 'casr']
- implicit: ['buffers_Myosin_ca', 'buffers_Myosin_mg']
- implicit: ['buffers_SLLj', 'caj', 'ryr_CaRI', 'ryr_I', 'ryr_O', 'ryr_R']
- implicit: ['buffers_SLLsl', 'casl']
- implicit: ['buffers_TnCHc', 'buffers_TnCHm']
- implicit: ['cai', 'contraction_Ca_TRPN']
- implicit: ['contraction_TmBlocked', 'contraction_XS', 'contraction_XW']
- implicit: ['ryr_CaRI_p', 'ryr_I_p', 'ryr_O_p', 'ryr_R_p']
- sub-step: ['v', 'h', 'hp', 'j', 'jp', 'm'] (proportional, h_min 0.001, [{'state': 'v', 'atol': 0.5, 'rtol': 0.0}])

Reasons: buffers_Myosin_ca,buffers_Myosin_mg: linear subsystem (the states are affine in each other); buffers_TnCHc,buffers_TnCHm: linear subsystem (the states are affine in each other); contraction_TmBlocked,contraction_XS,contraction_XW: linear subsystem (the states are affine in each other); ryr_CaRI_p,ryr_I_p,ryr_O_p,ryr_R_p: linear subsystem (the states are affine in each other); C0,C1,C2,I_,O_: linear subsystem (the states are affine in each other); stiff coupled pairs (kappa > 0.5 and dt*(|J_ii| + |J_jj|) > 1 at some sample); buffers_Csqn,casr: stiff coupled pairs (kappa > 0.5 and dt*(|J_ii| + |J_jj|) > 1 at some sample); buffers_SLLj,caj,ryr_CaRI,ryr_I,ryr_O,ryr_R: linear subsystem (the states are affine in each other); stiff coupled pairs (kappa > 0.5 and dt*(|J_ii| + |J_jj|) > 1 at some sample); buffers_SLLsl,casl: stiff coupled pairs (kappa > 0.5 and dt*(|J_ii| + |J_jj|) > 1 at some sample); cai,contraction_Ca_TRPN: stiff coupled pairs (kappa > 0.5 and dt*(|J_ii| + |J_jj|) > 1 at some sample); v,h,hp,j,jp,m: V/INa block (membrane potential v; gates with |J_gg|*dt > 0.1 and |J_Vg|*range_g*dt > 1 mV within 2 ms of the upstroke at t = 0.88 ms)

Differences: only here: ['buffers_SLHsl', 'buffers_SLLsl', 'casl']; only here: ['ryr_CaRI', 'ryr_I', 'ryr_O', 'ryr_R']; only here: ['buffers_Csqn', 'buffers_SLHj', 'buffers_SLLj', 'caj', 'casr']; only gotranx: ['contraction_TmBlocked', 'contraction_XS', 'contraction_XW']; only gotranx: ['buffers_Csqn', 'casr']; only gotranx: ['buffers_Myosin_ca', 'buffers_Myosin_mg']; only gotranx: ['buffers_SLLsl', 'casl']; only gotranx: ['cai', 'contraction_Ca_TRPN']; only gotranx: ['buffers_SLLj', 'caj', 'ryr_CaRI', 'ryr_I', 'ryr_O', 'ryr_R']; only gotranx: ['buffers_TnCHc', 'buffers_TnCHm']
gotranx does not detect the dyad; its V block comes from the numeric fast-gate rule, which is the cross-check of the INa gate list.
