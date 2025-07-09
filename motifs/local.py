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
