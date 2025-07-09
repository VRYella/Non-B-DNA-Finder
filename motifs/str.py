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
