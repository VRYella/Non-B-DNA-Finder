import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

def non_overlapping_finditer(pattern, seq):
    regex = re.compile(pattern)
    i = 0
    while i < len(seq):
        m = regex.search(seq, i)
        if not m:
            break
        yield m
        i = m.end()

def find_gquadruplex(seq):
    pattern = r"(G{3,}(?:[ATGC]{0,7}G{3,}){3})"
    return [
        dict(Class="Quadruplex", Subtype="Canonical_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_relaxed_gquadruplex(seq):
    pattern = r"(G{3,}(?:[ATGC]{0,12}G{3,}){3})"
    return [
        dict(Class="Quadruplex", Subtype="Relaxed_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_bulged_gquadruplex(seq):
    pattern = r"(G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,})"
    return [
        dict(Class="Quadruplex", Subtype="Bulged_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter (bulge)", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_imotif(seq):
    pattern = r"(C{3,}(?:[ATGC]{0,7}C{3,}){3})"
    return [
        dict(Class="Quadruplex", Subtype="i-Motif", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{-g4hunter_score(m.group().replace('C','G')):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_gtriplex(seq):
    pattern = r"(G{3,}(?:[ATGC]{0,7}G{3,}){2})"
    return [
        dict(Class="Triplex", Subtype="G-Triplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_bipartite_gquadruplex(seq):
    pattern = r"(G{3,}(?:[ATGC]{0,30}G{3,}){3})"
    return [
        dict(Class="Quadruplex", Subtype="Bipartite_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_multimeric_gquadruplex(seq):
    pattern = r"(G{3,}(?:[ATGC]{0,12}G{3,}){4,})"
    return [
        dict(Class="Quadruplex", Subtype="Multimeric_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def all_motifs(seq):
    motif_funcs = [
        find_gquadruplex, find_relaxed_gquadruplex, find_bulged_gquadruplex,
        find_imotif, find_gtriplex, find_bipartite_gquadruplex, find_multimeric_gquadruplex
    ]
    seen = set()
    all_hits = []
    for func in motif_funcs:
        for hit in func(seq):
            key = (hit['Start'], hit['End'], hit['Subtype'])
            if key not in seen:
                seen.add(key)
                all_hits.append(hit)
    return all_hits
