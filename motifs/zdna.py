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
