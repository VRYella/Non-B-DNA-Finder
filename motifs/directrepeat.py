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
