import re
from .utils import wrap, gc_content

def find_direct_repeats(seq):
    results = []
    for unit_len in range(10, 31):
        pattern = rf'([ATGC]{{{unit_len}}})([ATGC]{{0,100}})\1'
        for m in re.finditer(pattern, seq):
            region = seq[m.start():m.end()]
            results.append(dict(
                Class="Direct Repeat",
                Subtype="Slipped_DNA",
                Start=m.start()+1,
                End=m.end(),
                Length=len(region),
                Sequence=wrap(region),
                GC=f"{gc_content(region):.1f}",
                Score=unit_len,
                ScoreMethod="Unit length"
            ))
    return results
import re
from .utils import wrap, gc_content

def find_local_bends(seq):
    pattern = r'A{6,7}|T{6,7}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Local Bend",
            Subtype="A/T-tract",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score="NA",
            ScoreMethod="A/T-tract"
        ))
    return results

def find_local_flexible(seq):
    pattern = r'(?:CA){4,}|(?:TG){4,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Local Flexibility",
            Subtype="CA/TG_dinucleotide",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score="NA",
            ScoreMethod="Dinucleotide"
        ))
    return results
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
from .utils import wrap, gc_content

def find_str(seq):
    results = []
    n = len(seq)
    for unit_len in range(1, 10):
        i = 0
        while i < n - unit_len:
            unit = seq[i:i+unit_len]
            reps = 1
            j = i + unit_len
            while j + unit_len <= n and seq[j:j+unit_len] == unit:
                reps += 1
                j += unit_len
            total_len = reps * unit_len
            if reps >= 2 and total_len >= 10:
                region = seq[i:j]
                results.append(dict(
                    Class="STR",
                    Subtype=f"STR_{unit_len}bp",
                    Start=i+1,
                    End=j,
                    Length=len(region),
                    Sequence=wrap(region),
                    GC=f"{gc_content(region):.1f}",
                    Score=reps,
                    ScoreMethod="Repeat count"
                ))
                i = j
            else:
                i += 1
    return results
import re
from .utils import wrap, gc_content

def find_hdna(seq):
    pattern = r'([AG]{10,}|[CT]{10,})([ATGC]{0,8})([AG]{10,}|[CT]{10,})'
    results = []
    for m in re.finditer(pattern, seq):
        left, spacer, right = m.group(1), m.group(2), m.group(3)
        if left == right[::-1] and len(left) >= 10 and len(right) >= 10:
            region = seq[m.start():m.end()]
            results.append(dict(
                Class="Triplex",
                Subtype="H-DNA",
                Start=m.start()+1,
                End=m.end(),
                Length=len(region),
                Sequence=wrap(region),
                GC=f"{gc_content(region):.1f}",
                Score="NA",
                ScoreMethod="Mirror repeat"
            ))
    return results

def find_sticky_dna(seq):
    pattern = r'(?:GAA){5,}|(?:TTC){5,}'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        results.append(dict(
            Class="Triplex",
            Subtype="Sticky_DNA",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score=region.count("GAA")+region.count("TTC"),
            ScoreMethod="Repeat count"
        ))
    return results
import re
from .utils import wrap, gc_content

def z_seeker_score(seq):
    return len(re.findall(r'GC|CG|GT|TG|AC|CA', seq))

def find_zdna(seq):
    pattern = r'((?:GC|CG|GT|TG|AC|CA){6,})'
    results = []
    for m in re.finditer(pattern, seq):
        region = seq[m.start():m.end()]
        score = z_seeker_score(region)
        results.append(dict(
            Class="Z-DNA",
            Subtype="Z-DNA",
            Start=m.start()+1,
            End=m.end(),
            Length=len(region),
            Sequence=wrap(region),
            GC=f"{gc_content(region):.1f}",
            Score=score,
            ScoreMethod="Z-Seeker"
        ))
    return results
