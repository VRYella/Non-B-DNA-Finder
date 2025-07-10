import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

def overlapping_finditer(pattern, seq):
    return re.compile(pattern).finditer(seq)

def create_motif_dict(cls, subtype, match, seq, score_method="None", score="0", group=0):
    """Helper to create standardized motif dictionary"""
    sequence = match.group(group)
    return {
        "Class": cls, "Subtype": subtype, "Start": match.start() + 1,
        "End": match.start() + len(sequence), "Length": len(sequence),
        "Sequence": wrap(sequence), "ScoreMethod": score_method, "Score": score
    }

def find_motif(seq, pattern, cls, subtype, score_method="None", score_func=None, group=0):
    """Generic motif finder"""
    results = []
    for m in overlapping_finditer(pattern, seq):
        if score_func:
            score = f"{score_func(m.group(group)):.2f}"
        else:
            score = "0"
        results.append(create_motif_dict(cls, subtype, m, seq, score_method, score, group))
    return results

# G-Quadruplex variants
def find_gquadruplex(seq):
    return find_motif(seq, r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))", 
                     "Quadruplex", "Canonical_G-Quadruplex", "G4Hunter", g4hunter_score)

def find_relaxed_gquadruplex(seq):
    return find_motif(seq, r"(?=(G{3,}(?:[ATGC]{0,12}G{3,}){3}))", 
                     "Quadruplex", "Relaxed_G-Quadruplex", "G4Hunter", g4hunter_score)

def find_bulged_gquadruplex(seq):
    return find_motif(seq, r"(?=(G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}))", 
                     "Quadruplex", "Bulged_G-Quadruplex", "G4Hunter (bulge)", g4hunter_score, 1)

def find_bipartite_gquadruplex(seq):
    return find_motif(seq, r"(?=(G{3,}(?:[ATGC]{0,30}G{3,}){3}))", 
                     "Quadruplex", "Bipartite_G-Quadruplex", "G4Hunter", g4hunter_score)

def find_multimeric_gquadruplex(seq):
    return find_motif(seq, r"(?=((G{3,}(?:[ATGC]{0,12}G{3,}){4,})))", 
                     "Quadruplex", "Multimeric_G-Quadruplex", "G4Hunter", g4hunter_score, 1)

def find_imotif(seq):
    return find_motif(seq, r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))", 
                     "Quadruplex", "i-Motif", "G4Hunter", 
                     lambda x: -g4hunter_score(x.replace('C','G')))

def find_gtriplex(seq):
    return find_motif(seq, r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))", 
                     "Triplex", "G-Triplex", "G4Hunter", g4hunter_score)

def find_zdna(seq):
    return find_motif(seq, r"(?=((?:CG){6,}))", 
                     "Z-DNA", "CG_Repeat", "ZSeeker", zseeker_score, 1)

# Simple motifs (no scoring)
SIMPLE_MOTIFS = [
    (r"(?=(T{3,}[ATGC]{1,7}A{3,}))", "H-DNA", "T-A"),
    (r"(?=(CTGCTGCTGCTG))", "Sticky_DNA", "CTG"),
    (r"(?=((?:AT){6,}))", "Slipped_DNA", "AT_Slippage"),
    (r"(?=(A{4,}TTTT))", "Cruciform", "A-T"),
    (r"(?=(A{6,7}|T{6,7}))", "Bent_DNA", "Poly-A/T"),
    (r"(?=(?:AAATT){2,})", "A-Phased_Repeat", "APR"),
    (r"(?=(ATCGCGAT))", "Mirror_Repeat", "ATCGCGAT"),
    (r"(?=(G{6,}))", "Direct_Repeat", "Poly-G"),
]

def find_simple_motifs(seq):
    """Find all simple motifs that don't require scoring"""
    results = []
    for pattern, cls, subtype in SIMPLE_MOTIFS:
        group = 1 if pattern.count('(') > 1 else 0
        results.extend(find_motif(seq, pattern, cls, subtype, group=group))
    return results

def find_local_bent(seq):
    return find_motif(seq, r"(?=(A{6,7}|T{6,7}))", "Bent_DNA", "Poly-A/T", group=1)

def find_overlap_hybrid(seq, pattern1, pattern2, cls, subtype):
    """Generic function to find overlapping motifs"""
    hits = []
    hits1 = [(m.start(), m.end()) for m in overlapping_finditer(pattern1, seq)]
    for m in overlapping_finditer(pattern2, seq):
        for start1, end1 in hits1:
            if m.start() < end1 and m.end() > start1:
                hits.append(create_motif_dict(cls, subtype, m, seq))
                break
    return hits

def find_quadruplex_triplex_hybrid(seq):
    return find_overlap_hybrid(seq, r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))", 
                              r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))", 
                              "Hybrid", "G4-Triplex")

def find_cruciform_triplex_junction(seq):
    return find_overlap_hybrid(seq, r"(?=(A{4,}TTTT))", 
                              r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))", 
                              "Junction", "Cruciform-Triplex")

def find_g4_imotif_hybrid(seq):
    return find_overlap_hybrid(seq, r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))", 
                              r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))", 
                              "Hybrid", "G4-i-Motif")

def find_polyG(seq):
    return find_motif(seq, r"(?=(G{6,}))", "Direct_Repeat", "Poly-G", group=1)

def all_motifs(seq):
    """Find all motifs in sequence"""
    motif_funcs = [
        find_gquadruplex, find_relaxed_gquadruplex, find_bulged_gquadruplex,
        find_imotif, find_gtriplex, find_bipartite_gquadruplex, find_multimeric_gquadruplex,
        find_zdna, find_simple_motifs, find_quadruplex_triplex_hybrid, 
        find_cruciform_triplex_junction, find_g4_imotif_hybrid, find_polyG, find_local_bent
    ]
    return [hit for func in motif_funcs for hit in func(seq)]

def find_hotspots(seq, motif_hits, window=100, min_count=3):
    """Find regions with high motif density"""
    motif_positions = [(hit["Start"], hit["End"]) for hit in motif_hits]
    hotspots = {}
    
    for i in range(1, len(seq) - window + 2):
        region_end = i + window - 1
        count = sum(start <= region_end and end >= i for start, end in motif_positions)
        if count >= min_count:
            hotspots[(i, region_end)] = {"RegionStart": i, "RegionEnd": region_end, "MotifCount": count}
    
    return list(hotspots.values())
