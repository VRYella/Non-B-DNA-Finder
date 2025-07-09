import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
from datetime import datetime

from motifs import (
    parse_fasta, find_gquadruplex, find_imotif, find_bipartite_g4, find_gtriplex,
    find_zdna, find_cruciform, find_hdna, find_sticky_dna, find_direct_repeats,
    find_local_bends, find_local_flexible, find_str
)

EXAMPLE_FASTA = """>Example
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"""

st.set_page_config(page_title="Non-B DNA Motif Finder", layout="wide")

PAGES = ["Home", "Upload & Analyze", "Results", "Visualization", "Download"]
page = st.sidebar.radio("Navigation", PAGES)

def collect_all_motifs(seq):
    g4 = find_gquadruplex(seq)
    g4_spans = [(m['Start']-1, m['End']) for m in g4]
    motifs = (
        g4
        + find_imotif(seq)
        + find_bipartite_g4(seq)
        + find_gtriplex(seq, g4_spans)
        + find_zdna(seq)
        + find_cruciform(seq)
        + find_hdna(seq)
        + find_sticky_dna(seq)
        + find_direct_repeats(seq)
        + find_local_bends(seq)
        + find_local_flexible(seq)
        + find_str(seq)
    )
    motifs.sort(key=lambda x: (x['Start'], -x['Length']))
    mask = set()
    nonoverlap = []
    for m in motifs:
        s, e = m['Start'], m['End']
        if not any(i in mask for i in range(s, e+1)):
            nonoverlap.append(m)
            mask.update(range(s, e+1))
    return nonoverlap

if page == "Home":
    st.title("Non-B DNA Motif Finder")
    st.image("nbd.PNG", use_container_width=True)
    st.markdown("""
    Comprehensive, modular, and fast non-B DNA motif finder.
    **Motifs:** G-quadruplex, i-Motif, Bipartite G4, G-Triplex, Z-DNA (Z-Seeker), Cruciform, H-DNA, Sticky DNA, Direct/Mirror Repeats, STRs, local bends, and more.
    """)

elif page == "Upload & Analyze":
    st.header("Input Sequence")
    col1, col2 = st.columns([1,1])
    with col1:
        fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
        if fasta_file:
            try:
                seq = parse_fasta(fasta_file.read().decode("utf-8"))
                st.session_state['seq'] = seq
                st.success("FASTA file loaded!")
            except Exception:
                st.error("Could not parse file as UTF-8 or FASTA.")
    with col2:
        if st.button("Use Example Sequence"):
            st.session_state['seq'] = parse_fasta(EXAMPLE_FASTA)
        seq_input = st.text_area("Paste sequence (FASTA or raw)", value=st.session_state.get('seq', ""), height=120)
        if seq_input:
            try:
                seq = parse_fasta(seq_input)
                st.session_state['seq'] = seq
            except Exception:
                st.error("Paste a valid FASTA or sequence.")

    if st.button("Run Analysis"):
        seq = st.session_state.get('seq', "")
        if not seq or not re.match("^[ATGC]+$", seq):
            st.error("Please upload or paste a valid DNA sequence (A/T/G/C only).")
        else:
            with st.spinner("Analyzing sequence ..."):
                results = collect_all_motifs(seq)
                st.session_state['motif_results'] = results
                st.session_state['df'] = pd.DataFrame(results)
            if not results:
                st.warning("No non-B DNA motifs detected in this sequence.")
            else:
                st.success(f"Detected {len(results)} motif region(s) in {len(seq):,} bp.")

elif page == "Results":
    st.header("Motif Detection Results")
    df = st.session_state.get('df', pd.DataFrame())
    if df.empty:
        st.info("No results yet. Go to 'Upload & Analyze' and run analysis.")
    else:
        st.markdown(f"**Sequence length:** {len(st.session_state['seq']):,} bp")
        st.dataframe(df[['Class', 'Subtype', 'Start', 'End', 'Length', 'GC', 'ScoreMethod', 'Score', 'Sequence']],
            use_container_width=True, hide_index=True)
        with st.expander("Motif Class Summary"):
            motif_counts = df["Subtype"].value
