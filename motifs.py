import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

def non_overlapping_finditer(pattern, seq):
    i = 0
    n = len(seq)
    regex = re.compile(pattern)
    while i < n:
        m = regex.search(seq, i)
        if not m: break
        yield m
        i = m.end()

# --- NEW, more specific G4 finders ---
def find_gquadruplex(seq):
    # Canonical G4: (G3+N1-7)3G3+
    pattern = r"(G{3,}[ATGC]{1,7}){3}G{3,}"
    return [
        dict(Class="Quadruplex", Subtype="Canonical_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_relaxed_gquadruplex(seq):
    # Relaxed G4: (G3+N1-12)3G3+
    pattern = r"(G{3,}[ATGC]{1,12}){3}G{3,}"
    return [
        dict(Class="Quadruplex", Subtype="Relaxed_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_bulged_gquadruplex(seq):
    # Bulged G4: allow bulges (up to 3nt) in G-runs
    pattern = r"G{3,}(?:[ATGC]{0,3}G{3,}){3,}"
    return [
        dict(Class="Quadruplex", Subtype="Bulged_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter (bulge)", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

# --- Existing motif finders from your file, unchanged for backward compatibility ---
def find_imotif(seq):
    pattern = r"(C{3,}[ATGC]{1,7}){3}C{3,}"
    return [
        dict(Class="Quadruplex", Subtype="i-Motif", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="G4Hunter", Score=f"{-g4hunter_score(m.group().replace('C','G')):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_gtriplex(seq):
    g4_pat = re.compile(r"(G{3,}[ATGC]{1,7}){3}G{3,}")
    pattern = re.compile(r"(G{3,}[ATGC]{1,7}){2}G{3,}")
    found = []
    i = 0
    n = len(seq)
    while i < n:
        m = pattern.search(seq, i)
        if not m: break
        region = m.group()
        if not g4_pat.match(region):
            found.append(dict(
                Class="Quadruplex", Subtype="G-Triplex", Start=m.start()+1, End=m.end(),
                Length=m.end()-m.start(), Sequence=wrap(region), ScoreMethod="G4Hunter*0.7",
                Score=f"{g4hunter_score(region)*0.7:.2f}"
            ))
        i = m.end()
    return found

def find_bipartite_gquadruplex(seq):
    pattern = r"((G{3,}[ATGC]{1,7}){3}G{3,})[ATGC]{0,100}((G{3,}[ATGC]{1,7}){3}G{3,})"
    return [
        dict(Class="Quadruplex", Subtype="Bipartite_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="Mean G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_multimeric_gquadruplex(seq):
    pattern = r"((G{3,}[ATGC]{1,7}){3}G{3,}(?:[ATGC]{1,50}(G{3,}[ATGC]{1,7}){3}G{3,})+)"
    return [
        dict(Class="Quadruplex", Subtype="Multimeric_G-Quadruplex", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="Mean G4Hunter", Score=f"{g4hunter_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_zdna(seq):
    pattern = r"((?:GC|CG|GT|TG|AC|CA){5,})"
    return [
        dict(Class="Z-DNA", Subtype="Z-DNA", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="Z-Seeker", Score=f"{zseeker_score(m.group()):.2f}")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_hdna(seq):
    pattern = r"([AG]{10,}|[CT]{10,})([ATGC]{0,8})([AG]{10,}|[CT]{10,})"
    return [
        dict(Class="Triplex", Subtype="H-DNA", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="NA", Score="NA")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_sticky_dna(seq):
    pattern = r"(?:GAA){5,}|(?:TTC){5,}"
    return [
        dict(Class="Triplex", Subtype="Sticky_DNA", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="RepeatCount", Score=str(len(m.group())//3))
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_slipped_dna(seq):
    pattern = r"([ATGC]{10,25})([ATGC]{0,10})\1"
    return [
        dict(Class="Direct Repeat", Subtype="Slipped_DNA", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="RepeatUnit", Score=str(len(m.group())))
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_cruciform(seq):
    matches = []
    n = len(seq)
    for arm in range(6, 21):
        for loop in range(0, 101):
            pattern = rf"([ATGC]{{{arm}}})([ATGC]{{0,{loop}}})([ATGC]{{{arm}}})"
            regex = re.compile(pattern)
            i = 0
            while i < n:
                m = regex.search(seq, i)
                if not m: break
                left = m.group(1)
                right = m.group(3)
                if reverse_complement(left) == right:
                    region = seq[m.start():m.end()]
                    matches.append(dict(
                        Class="Inverted Repeat", Subtype="Cruciform_DNA", Start=m.start()+1, End=m.end(), Length=len(region),
                        Sequence=wrap(region), ScoreMethod="Arm length", Score=arm
                    ))
                    i = m.end()
                else:
                    i += 1
    return matches

def find_bent_dna(seq):
    pattern = r"(A{3,11})([ATGC]{3,11})(A{3,11})([ATGC]{3,11})(A{3,11})"
    return [
        dict(Class="Bent DNA", Subtype="Bent_DNA", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="Tract length", Score=str(len(m.group(1))))
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_apr(seq):
    pattern = r"(A{3,11})([ATGC]{3,11})(A{3,11})([ATGC]{3,11})(A{3,11})"
    return [
        dict(Class="Bent DNA", Subtype="APR", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="NA", Score="NA")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_mirror_repeat(seq):
    pattern = r"([ATGC]{10,})([ATGC]{0,100})([ATGC]{10,})"
    results = []
    for m in non_overlapping_finditer(pattern, seq):
        left = m.group(1)
        right = m.group(3)
        if left == right[::-1]:
            region = seq[m.start():m.end()]
            results.append(dict(
                Class="Mirror Repeat", Subtype="Mirror_Repeat", Start=m.start()+1, End=m.end(), Length=len(region),
                Sequence=wrap(region), ScoreMethod="NA", Score="NA"
            ))
    return results

def find_quadruplex_triplex_hybrid(seq):
    pattern = r"(G{3,}[ATGC]{1,7}){3}G{3,}[ATGC]{0,100}([AG]{10,}|[CT]{10,})"
    return [
        dict(Class="Quadruplex-Triplex Hybrid", Subtype="Quadruplex-Triplex_Hybrid", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="NA", Score="NA")
        for m in non_overlapping_finditer(pattern, seq)
    ]

def find_cruciform_triplex_junction(seq):
    pattern = r"([ATGC]{10,})([ATGC]{0,100})([ATGC]{10,})([ATGC]{0,100})([AG]{10,}|[CT]{10,})"
    return [
        dict(Class="Cruciform-Triplex Junction", Subtype="Cruciform-Triplex_Junctions", Start=m.start()+1, End=m.end(), Length=m.end()-m.start(),
             Sequence=wrap(m.group()), ScoreMethod="NA", Score="NA")
        for m in non_overlapping_finditer(pattern, seq)
    ]

# --- Updated hybrid finder (stricter window) ---
def find_g4_imotif_hybrid(seq):
    # G4 and i-motif within 20 nt (updated from 100)
    g4s = [(m.start(), m.end()) for m in re.finditer(r"(G{3,}[ATGC]{1,7}){3}G{3,}", seq)]
    imots = [(m.start(), m.end()) for m in re.finditer(r"(C{3,}[ATGC]{1,7}){3}C{3,}", seq)]
    motifs = []
    for g4s_start, g4s_end in g4s:
        for ims_start, ims_end in imots:
            if abs(g4s_end - ims_start) <= 20 or abs(ims_end - g4s_start) <= 20:
                region = seq[min(g4s_start, ims_start):max(g4s_end, ims_end)]
                motifs.append(dict(
                    Class="Hybrid", Subtype="G4-iMotif_Hybrid",
                    Start=min(g4s_start, ims_start)+1, End=max(g4s_end, ims_end),
                    Length=len(region), Sequence=wrap(region),
                    ScoreMethod="NA", Score="NA"
                ))
    return motifs
