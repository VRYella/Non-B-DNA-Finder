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
