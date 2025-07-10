import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern)
    for m in regex.finditer(seq):
        yield m

def find_gquadruplex(seq):
    # Allow loops of length 0–7 (so uninterrupted G runs are detected)
    pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))"
    return [
        dict(Class="Quadruplex", Subtype="Canonical_G-Quadruplex", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group(0)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_relaxed_gquadruplex(seq):
    pattern = r"(?=(G{3,}(?:[ATGC]{0,12}G{3,}){3}))"
    return [
        dict(Class="Quadruplex", Subtype="Relaxed_G-Quadruplex", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group(0)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bulged_gquadruplex(seq):
    # Allow up to 3 non-Gs in G runs
    pattern = r"(?=(G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}))"
    return [
        dict(Class="Quadruplex", Subtype="Bulged_G-Quadruplex", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="G4Hunter (bulge)", Score=f"{g4hunter_score(m.group(1)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_imotif(seq):
    pattern = r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))"
    return [
        dict(Class="Quadruplex", Subtype="i-Motif", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="G4Hunter", Score=f"{-g4hunter_score(m.group(0).replace('C','G')):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_gtriplex(seq):
    pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))"
    return [
        dict(Class="Triplex", Subtype="G-Triplex", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group(0)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bipartite_gquadruplex(seq):
    pattern = r"(?=(G{3,}(?:[ATGC]{0,30}G{3,}){3}))"
    return [
        dict(Class="Quadruplex", Subtype="Bipartite_G-Quadruplex", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group(0)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_multimeric_gquadruplex(seq):
    # At least 5 G runs with up to 12 bases in between
    pattern = r"(?=((G{3,}(?:[ATGC]{0,12}G{3,}){4,})))"
    return [
        dict(Class="Quadruplex", Subtype="Multimeric_G-Quadruplex", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group(1)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_zdna(seq):
    pattern = r"(?=((?:CG){6,}))"
    return [
        dict(Class="Z-DNA", Subtype="CG_Repeat", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="ZSeeker", Score=f"{zseeker_score(m.group(1)):.2f}")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_hdna(seq):
    pattern = r"(?=(T{3,}[ATGC]{1,7}A{3,}))"
    return [
        dict(Class="H-DNA", Subtype="T-A", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_sticky_dna(seq):
    pattern = r"(?=(CTGCTGCTGCTG))"
    return [
        dict(Class="Sticky_DNA", Subtype="CTG", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_slipped_dna(seq):
    pattern = r"(?=((?:AT){6,}))"
    return [
        dict(Class="Slipped_DNA", Subtype="AT_Slippage", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_cruciform(seq):
    pattern = r"(?=(A{4,}TTTT))"
    return [
        dict(Class="Cruciform", Subtype="A-T", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bent_dna(seq):
    # Local bent regions: polyA or polyT (6 or 7)
    pattern = r"(?=(A{6,7}|T{6,7}))"
    return [
        dict(Class="Bent_DNA", Subtype="Poly-A/T", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_apr(seq):
    pattern = r"(?=(?:AAATT){2,})"
    return [
        dict(Class="A-Phased_Repeat", Subtype="APR", Start=m.start()+1, End=m.start()+len(m.group(0)), Length=len(m.group(0)),
             Sequence=wrap(m.group(0)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_mirror_repeat(seq):
    pattern = r"(?=(ATCGCGAT))"
    return [
        dict(Class="Mirror_Repeat", Subtype="ATCGCGAT", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_quadruplex_triplex_hybrid(seq):
    g4_pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))"
    triplex_pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))"
    hits = []
    g4_hits = [(m.start(), m.end()) for m in overlapping_finditer(g4_pattern, seq)]
    for m in overlapping_finditer(triplex_pattern, seq):
        for g4_start, g4_end in g4_hits:
            if (m.start() < g4_end and m.end() > g4_start):
                hits.append(dict(Class="Hybrid", Subtype="G4-Triplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
                                 Sequence=wrap(seq[m.start():m.end()]), ScoreMethod="None", Score="0"))
                break
    return hits

def find_cruciform_triplex_junction(seq):
    cruciform_pattern = r"(?=(A{4,}TTTT))"
    triplex_pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))"
    hits = []
    cruciform_hits = [(m.start(), m.end()) for m in overlapping_finditer(cruciform_pattern, seq)]
    for m in overlapping_finditer(triplex_pattern, seq):
        for c_start, c_end in cruciform_hits:
            if (m.start() < c_end and m.end() > c_start):
                hits.append(dict(Class="Junction", Subtype="Cruciform-Triplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
                                 Sequence=wrap(seq[m.start():m.end()]), ScoreMethod="None", Score="0"))
                break
    return hits

def find_g4_imotif_hybrid(seq):
    g4_pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))"
    imotif_pattern = r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))"
    hits = []
    g4_hits = [(m.start(), m.end()) for m in overlapping_finditer(g4_pattern, seq)]
    for m in overlapping_finditer(imotif_pattern, seq):
        for g4_start, g4_end in g4_hits:
            if (m.start() < g4_end and m.end() > g4_start):
                hits.append(dict(Class="Hybrid", Subtype="G4-i-Motif", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
                                 Sequence=wrap(seq[m.start():m.end()]), ScoreMethod="None", Score="0"))
                break
    return hits

def find_polyG(seq):
    pattern = r"(?=(G{6,}))"
    return [
        dict(Class="Direct_Repeat", Subtype="Poly-G", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def find_local_bent(seq):
    # Poly-A/T as local bend, but let's keep this for clarity
    pattern = r"(?=(A{6,7}|T{6,7}))"
    return [
        dict(Class="Bent_DNA", Subtype="Poly-A/T", Start=m.start()+1, End=m.start()+len(m.group(1)), Length=len(m.group(1)),
             Sequence=wrap(m.group(1)), ScoreMethod="None", Score="0")
        for m in overlapping_finditer(pattern, seq)
    ]

def all_motifs(seq):
    motif_funcs = [
        find_gquadruplex, find_relaxed_gquadruplex, find_bulged_gquadruplex,
        find_imotif, find_gtriplex, find_bipartite_gquadruplex, find_multimeric_gquadruplex,
        find_zdna, find_hdna, find_sticky_dna, find_slipped_dna, find_cruciform,
        find_bent_dna, find_apr, find_mirror_repeat,
        find_quadruplex_triplex_hybrid, find_cruciform_triplex_junction, find_g4_imotif_hybrid,
        find_polyG,  # Poly-G direct repeat
        find_local_bent
    ]
    all_hits = []
    for func in motif_funcs:
        all_hits.extend(func(seq))
    return all_hits

def find_hotspots(seq, motif_hits, window=100, min_count=3):
    n = len(seq)
    hotspots = []
    motif_positions = [(hit["Start"], hit["End"]) for hit in motif_hits]
    for i in range(1, n-window+2):
        region_start = i
        region_end = i + window - 1
        count = sum(mstart <= region_end and mend >= region_start for mstart, mend in motif_positions)
        if count >= min_count:
            hotspots.append(dict(RegionStart=region_start, RegionEnd=region_end, MotifCount=count))
    # Remove duplicates by region
    unique_hotspots = {(h['RegionStart'], h['RegionEnd']): h for h in hotspots}
    return list(unique_hotspots.values())
