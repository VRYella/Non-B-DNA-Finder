import re
from .utils import wrap, gc_content

def g4hunter_score(seq):
    seq = seq.upper()
    vals = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == 'G':
            run_len = 1
            while i + run_len < n and seq[i + run_len] == 'G':
                run_len += 1
            val = min(run_len, 4)
            vals.extend([val] * run_len)
            i += run_len
        elif seq[i] == 'C':
            run_len = 1
            while i + run_len < n and seq[i + run_len] == 'C':
                run_len += 1
            val = -min(run_len, 4)
            vals.extend([val] * run_len)
            i += run_len
        else:
            vals.append(0)
            i += 1
    return round(np.mean(vals), 2) if vals else 0

def find_gquadruplex(seq):
    pattern = r'(G{3,}[ATGC]{1,7}){3}G{3,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Quadruplex",
            Subtype="G-Quadruplex",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score=g4hunter_score(region),
            ScoreMethod="G4Hunter"
        ))
    return results

def find_imotif(seq):
    pattern = r'(C{3,}[ATGC]{1,7}){3}C{3,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Quadruplex",
            Subtype="i-Motif",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score=-g4hunter_score(region.replace("G", "C").replace("C", "G")),
            ScoreMethod="G4Hunter"
        ))
    return results

def find_bipartite_g4(seq):
    pattern = r'(G{3,}[ATGC]{1,7}){3}G{3,}[ATGC]{0,100}(G{3,}[ATGC]{1,7}){3}G{3,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Quadruplex",
            Subtype="Bipartite_G-Quadruplex",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score=g4hunter_score(region),
            ScoreMethod="G4Hunter"
        ))
    return results

def find_gtriplex(seq, g4_spans):
    pattern = r'(G{3,}[ATGC]{1,7}){2}G{3,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        s, e = m.start(), m.end()
        overlaps_g4 = any((s < g4e and e > g4s) for g4s, g4e in g4_spans)
        if not overlaps_g4:
            results.append(dict(
                Class="Quadruplex",
                Subtype="G-Triplex",
                Start=s+1,
                End=e,
                Length=len(region),
                Sequence=wrap(region),
                GC=f"{gc_content(region):.1f}",
                Score=g4hunter_score(region) * 0.75,
                ScoreMethod="G4Hunter scaled"
            ))
    return results
