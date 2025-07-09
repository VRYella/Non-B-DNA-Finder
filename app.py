import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from motifs import (
    find_gquadruplex, find_imotif, find_gtriplex,
    find_bipartite_gquadruplex, find_multimeric_gquadruplex,
    find_zdna, find_hdna, find_sticky_dna,
    find_slipped_dna, find_cruciform, find_bent_dna,
    find_apr, find_mirror_repeat,
    find_quadruplex_triplex_hybrid,
    find_cruciform_triplex_junction,
    find_g4_imotif_hybrid
)
from utils import parse_fasta, wrap

EXAMPLE_FASTA = """>Example
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
"""

st.set_page_config(page_title="Non-B DNA Motif Finder", layout="wide")

# --- Home Page ---
st.title("Non-B DNA Motif Finder")
try:
    st.image("nbd.PNG", use_container_width=True)
except Exception:
    st.warning("Logo image (nbd.PNG) not found. Place it in the app folder.")

st.markdown("""
**Ultra-fast detection and visualization of non-B DNA motifs in DNA sequences.**
- Upload FASTA or paste sequence
- Finds: G-quadruplexes, i-motif, triplex, Z-DNA, cruciform, bent DNA, and more.
- Export CSV/Excel, see motif maps.
""")

st.subheader("Paste sequence, upload FASTA, or use example:")
col1, col2 = st.columns([2,1])
with col1:
    fasta_file = st.file_uploader("Upload FASTA", type=["fa", "fasta", "txt"])
    seq_text = st.text_area("Paste sequence (FASTA or raw)", height=120)
with col2:
    if st.button("Use Example Sequence"):
        seq_text = EXAMPLE_FASTA

seq = ""
if fasta_file is not None:
    try:
        seq = parse_fasta(fasta_file.read().decode("utf-8"))
    except Exception:
        st.error("Could not read FASTA file.")
elif seq_text:
    seq = parse_fasta(seq_text)

if seq:
    st.info(f"Sequence loaded ({len(seq):,} bp)")
    if st.button("Find Motifs"):
        with st.spinner("Analyzing..."):
            results = []
            results += find_gquadruplex(seq)
            results += find_imotif(seq)
            results += find_gtriplex(seq)
            results += find_bipartite_gquadruplex(seq)
            results += find_multimeric_gquadruplex(seq)
            results += find_zdna(seq)
            results += find_hdna(seq)
            results += find_sticky_dna(seq)
            results += find_slipped_dna(seq)
            results += find_cruciform(seq)
            results += find_bent_dna(seq)
            results += find_apr(seq)
            results += find_mirror_repeat(seq)
            results += find_quadruplex_triplex_hybrid(seq)
            results += find_cruciform_triplex_junction(seq)
            results += find_g4_imotif_hybrid(seq)
            # Remove overlapping motifs (keep highest score)
            results = sorted(results, key=lambda x: (x['Start'], -float(str(x.get('Score', 0)))))
            nonoverlap = []
            covered = set()
            for r in results:
                covered_range = set(range(r['Start'], r['End'] + 1))
                if not covered_range & covered:
                    nonoverlap.append(r)
                    covered |= covered_range
            df = pd.DataFrame(nonoverlap)
            if df.empty:
                st.warning("No non-B DNA motifs detected.")
            else:
                st.success(f"Found {len(df)} motif regions.")
                st.dataframe(df, use_container_width=True)
                # Export
                csv = df.to_csv(index=False).encode('utf-8')
                st.download_button("Download as CSV", data=csv, file_name="motifs.csv", mime="text/csv")
                try:
                    import io
                    import xlsxwriter
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                        df.to_excel(writer, index=False)
                    st.download_button("Download as Excel", data=output.getvalue(),
                        file_name="motifs.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                except ImportError:
                    st.info("Install 'xlsxwriter' for Excel export.")
                # Visualization
                st.subheader("Motif Map")
                motif_types = df['Subtype'].unique()
                color_palette = sns.color_palette('husl', n_colors=len(motif_types))
                color_map = {typ: color_palette[i] for i, typ in enumerate(motif_types)}
                y_map = {typ: i+1 for i, typ in enumerate(motif_types)}
                fig, ax = plt.subplots(figsize=(10, len(motif_types)*0.7+2))
                for _, motif in df.iterrows():
                    motif_type = motif['Subtype']
                    y = y_map[motif_type]
                    color = color_map[motif_type]
                    ax.hlines(y, motif['Start'], motif['End'], color=color, linewidth=8)
                ax.set_yticks(list(y_map.values()))
                ax.set_yticklabels(list(y_map.keys()))
                ax.set_xlim(0, len(seq)+1)
                ax.set_xlabel('Position (bp)')
                ax.set_title('Motif Map')
                sns.despine(left=False, bottom=False)
                plt.tight_layout()
                st.pyplot(fig)
else:
    st.info("Paste a sequence or upload a FASTA to begin.")

st.markdown("""
---
**Developed by [Your Name] & Team** | [GitHub](https://github.com/yourrepo)
""")
