import pandas as pd

import re
import math
from collections import Counter

def extract_features_batch_robust(df):
    def safe_gc(seq):
        return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0

    def shannon_entropy(seq):
        freq = Counter(seq)
        total = len(seq)
        return -sum((c / total) * math.log2(c / total) for c in freq.values()) if total > 0 else 0

    def calculate_mfe(seq, struct):
        stack = []
        mfe = 0
        for i, s in enumerate(struct):
            if s == '(':
                stack.append(i)
            elif s == ')' and stack:
                j = stack.pop()
                pair = seq[j] + seq[i]
                if pair in ['GC', 'CG']: mfe -= 3
                elif pair in ['AU', 'UA']: mfe -= 2
                elif pair in ['GU', 'UG']: mfe -= 1
        return mfe

    def pairing_fraction(seq, struct):
        return (struct.count('(') + struct.count(')')) / len(seq) if seq and struct and len(seq) == len(struct) else 0

    def dinucleotide_frequencies(seq):
        counts = Counter([seq[i:i+2] for i in range(len(seq) - 1)])
        total = sum(counts.values())
        return {f'di_{k}': v / total for k, v in counts.items()} if total > 0 else {}

    rows = []
    for _, row in df.iterrows():
        try:
            seq, struct = row['hairpin_seq'], row['hairpin_structure']
            caps = list(re.finditer(r'[A-Z]+', seq))
            if len(caps) != 2:
                continue
            span_5p = caps[0].span()
            span_3p = caps[1].span()
            seq_5p = seq[span_5p[0]:span_5p[1]]
            seq_3p = seq[span_3p[0]:span_3p[1]]
            struct_5p = struct[span_5p[0]:span_5p[1]]
            struct_3p = struct[span_3p[0]:span_3p[1]]
            combined = seq_5p + seq_3p

            feats = {
                'length': len(combined),
                'gc_content': safe_gc(combined),
                'au_content': (combined.count('A') + combined.count('U')) / len(combined),
                'five_p_mfe': calculate_mfe(seq_5p, struct_5p),
                'three_p_mfe': calculate_mfe(seq_3p, struct_3p),
                'gc_5p': safe_gc(seq_5p),
                'gc_3p': safe_gc(seq_3p),
                'shannon_entropy_5p': shannon_entropy(seq_5p),
                'shannon_entropy_3p': shannon_entropy(seq_3p),
                'pairing_frac_5p': pairing_fraction(seq_5p, struct_5p),
                'pairing_frac_3p': pairing_fraction(seq_3p, struct_3p),
                'length_diff': len(seq_5p) - len(seq_3p),
                'mfe_per_nt_5p': calculate_mfe(seq_5p, struct_5p) / len(seq_5p),
                'mfe_per_nt_3p': calculate_mfe(seq_3p, struct_3p) / len(seq_3p),
                '5p_A': seq_5p.count('A') / len(seq_5p),
                '5p_U': seq_5p.count('U') / len(seq_5p),
                '5p_G': seq_5p.count('G') / len(seq_5p),
                '5p_C': seq_5p.count('C') / len(seq_5p),
                '3p_A': seq_3p.count('A') / len(seq_3p),
                '3p_U': seq_3p.count('U') / len(seq_3p),
                '3p_G': seq_3p.count('G') / len(seq_3p),
                '3p_C': seq_3p.count('C') / len(seq_3p),
            }
            feats['GC_diff'] = feats['gc_5p'] - feats['gc_3p']
            feats['GC_ratio_5p3p'] = feats['gc_5p'] / feats['gc_3p'] if feats['gc_3p'] else 0
            feats['mfe_diff'] = feats['five_p_mfe'] - feats['three_p_mfe']
            feats['mfe_ratio'] = abs(feats['five_p_mfe'] / feats['three_p_mfe']) if feats['three_p_mfe'] else 0
            feats['mfe_per_nt_diff'] = feats['mfe_per_nt_5p'] - feats['mfe_per_nt_3p']

            for m in ['UGU', 'GAG', 'AAU', 'UAA', 'GUG', 'CGC', 'CCU']:
                feats[f'5p_has_{m}'] = 1 if m in seq_5p else 0
                feats[f'3p_has_{m}'] = 1 if m in seq_3p else 0

            feats.update(dinucleotide_frequencies(combined))

            for region, seq_str in zip(['5p', '3p'], [seq_5p, seq_3p]):
                for pos, idx in [('first', 0), ('last', -1)]:
                    base = seq_str[idx] if len(seq_str) > 0 else ''
                    for nt in 'AUGC':
                        feats[f"{region}_{pos}_nt_{nt}"] = 1 if base == nt else 0
                    feats[f"{region}_{pos}_nt_"] = 1 if base not in 'AUGC' else 0

            rows.append(feats)
        except:
            continue

    return pd.DataFrame(rows)
