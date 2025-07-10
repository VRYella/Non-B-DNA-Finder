import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re, io
from datetime import datetime
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

# ======= STYLES AND FONTS =======
st.markdown("""
    <link href="https://fonts.googleapis.com/css?family=Montserrat:400,700&display=swap" rel="stylesheet">
    <style>
        html, body, [class*="css"]  {
            font-family: 'Montserrat', sans-serif !important;
        }
        section[data-testid="stSidebar"] {
            background: linear-gradient(135deg, #e0eafc 0%, #cfdef3 100%);
            color: #222;
        }
        div[data-testid="stSidebar"] label {
            font-size: 22px !important;
            color: #1A5276 !important;
            font-family: 'Montserrat', sans-serif !important;
            font-weight: bold !important;
            padding: 6px 12px !important;
            border-radius: 8px !important;
        }
        div[data-testid="stSidebar"] .stRadio [role="radio"][aria-checked="true"] label {
            background: linear-gradient(90deg,#a1c4fd,#c2e9fb);
            color: #fff !important;
            box-shadow: 0 4px 16px rgba(161,196,253,0.2);
        }
        div[data-testid="stSidebar"] .stRadio [role="radio"]:hover label {
            background: #cfdef3;
            color: #2874A6 !important;
        }
        .sidebar-title {
            font-size: 28px;
            font-family: 'Montserrat', sans-serif;
            color: #3B5998;
            font-weight: bold;
            margin-bottom: 20px;
            margin-top: 10px;
            text-align: center;
            letter-spacing: 1px;
        }
    </style>
""", unsafe_allow_html=True)
# ======= END STYLES =======

st.set_page_config(page_title="Non-B DNA Motif Finder", layout="wide")

if 'seq' not in st.session_state:
    st.session_state['seq'] = ""
if 'df' not in st.session_state:
    st.session_state['df'] = pd.DataFrame()
if 'motif_results' not in st.session_state:
    st.session_state['motif_results'] = []
if 'analysis_status' not in st.session_state:
    st.session_state['analysis_status'] = ""
if 'stop_analysis' not in st.session_state:
    st.session_state['stop_analysis'] = False

PAGES = [
    "Home",
    "Upload & Analyze",
    "Results",
    "Visualization",
    "Download",
    "Additional Information"
]

st.sidebar.markdown('<div class="sidebar-title">🧬 Navigation</div>', unsafe_allow_html=True)
page = st.sidebar.radio("", PAGES, key="nav_radio")

def collect_all_motifs(seq, status_callback=None, stop_flag=None):
    results = []
    motif_steps = [
        ("Searching for G-Quadruplex motifs...", find_gquadruplex),
        ("Searching for i-Motif motifs...", find_imotif),
        ("Searching for G-Triplex motifs...", find_gtriplex),
        ("Searching for Bipartite G-Quadruplex motifs...", find_bipartite_gquadruplex),
        ("Searching for Multimeric G-Quadruplex motifs...", find_multimeric_gquadruplex),
        ("Searching for Z-DNA motifs...", find_zdna),
        ("Searching for H-DNA motifs...", find_hdna),
        ("Searching for Sticky DNA motifs...", find_sticky_dna),
        ("Searching for Slipped DNA motifs...", find_slipped_dna),
        ("Searching for Cruciform motifs...", find_cruciform),
        ("Searching for Bent DNA motifs...", find_bent_dna),
        ("Searching for APR motifs...", find_apr),
        ("Searching for Mirror Repeat motifs...", find_mirror_repeat),
        ("Searching for Quadruplex-Triplex Hybrid motifs...", find_quadruplex_triplex_hybrid),
        ("Searching for Cruciform-Triplex Junction motifs...", find_cruciform_triplex_junction),
        ("Searching for G4-iMotif Hybrid motifs...", find_g4_imotif_hybrid),
    ]
    for status, func in motif_steps:
        if stop_flag is not None and stop_flag():
            if status_callback is not None:
                status_callback("Motif search stopped by user.")
            return []
        if status_callback is not None:
            status_callback(status)
        results += func(seq)
    # Remove overlapping motifs (keep highest score if numeric)
    results = sorted(results, key=lambda x: (x['Start'], -float(x['Score']) if str(x['Score']).replace('.','',1).isdigit() else 0))
    nonoverlap = []
    covered = set()
    for r in results:
        covered_range = set(range(r['Start'], r['End'] + 1))
        if not covered_range & covered:
            nonoverlap.append(r)
            covered |= covered_range
    return nonoverlap

if page == "Home":
    st.title("Non-B DNA Motif Finder")
    try:
        st.image("nbd.PNG", use_container_width=True)
    except Exception:
        st.warning("Logo image (nbd.PNG) not found. Place it in the app folder.")

    st.markdown("""
    **Comprehensive, fast, and reference-grade non-B DNA motif finder.**
    - G-quadruplex, i-Motif, Bipartite G4, G-Triplex, Z-DNA (Z-Seeker), Cruciform, H-DNA, Sticky DNA, Direct/Mirror Repeats, STRs, local bends, flexible regions, and more.
    - Export to CSV/Excel, motif visualization included.
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

    # Run Analysis and Stop/Reset Buttons
    colA, colB, colC = st.columns([1, 0.6, 0.8])
    with colA:
        run_analysis = st.button("Run Analysis")
    with colB:
        stop_analysis = st.button("Stop/Reset")
    analysis_placeholder = st.empty()

    if run_analysis:
        st.session_state['stop_analysis'] = False
        seq = st.session_state.get('seq', "")
        if not seq or not re.match("^[ATGC]+$", seq):
            st.error("Please upload or paste a valid DNA sequence (A/T/G/C only).")
        else:
            progress_placeholder = st.empty()
            def update_status(msg):
                progress_placeholder.info(msg)
            def stop_flag():
                return st.session_state.get('stop_analysis', False)
            with st.spinner("Analyzing sequence ..."):
                results = collect_all_motifs(seq, status_callback=update_status, stop_flag=stop_flag)
                st.session_state['motif_results'] = results
                st.session_state['df'] = pd.DataFrame(results)
            progress_placeholder.empty()
            if st.session_state.get('stop_analysis', False):
                st.warning("Motif search was stopped.")
            elif not results:
                st.warning("No non-B DNA motifs detected in this sequence.")
            else:
                st.success(f"Detected {len(results)} motif region(s) in {len(seq):,} bp.")

    if stop_analysis:
        st.session_state['stop_analysis'] = True
        st.session_state['motif_results'] = []
        st.session_state['df'] = pd.DataFrame()
        analysis_placeholder.warning("Motif search has been reset. You may run a new analysis.")

elif page == "Results":
    st.header("Motif Detection Results")
    df = st.session_state.get('df', pd.DataFrame())
    if df.empty:
        st.info("No results yet. Go to 'Upload & Analyze' and run analysis.")
    else:
        st.markdown(f"**Sequence length:** {len(st.session_state['seq']):,} bp")
        st.dataframe(df[['Class', 'Subtype', 'Start', 'End', 'Length', 'Sequence', 'ScoreMethod', 'Score']],
            use_container_width=True, hide_index=True)
        with st.expander("Motif Class Summary"):
            motif_counts = df["Subtype"].value_counts().reset_index()
            motif_counts.columns = ["Motif Type", "Count"]
            st.dataframe(motif_counts, use_container_width=True, hide_index=True)

elif page == "Visualization":
    st.header("Motif Visualization")
    df = st.session_state.get('df', pd.DataFrame())
    seq = st.session_state.get('seq', "")
    if df.empty:
        st.info("No results to visualize. Run analysis first.")
    else:
        st.subheader("Motif Map (Full Sequence)")
        motif_types = sorted(df['Subtype'].unique())
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
        ax.set_xlabel('Position on Sequence (bp)')
        ax.set_title('Motif Map (Full Sequence)')
        sns.despine(left=False, bottom=False)
        plt.tight_layout()
        st.pyplot(fig)

        st.subheader("Motif Type Distribution (Pie Chart)")
        counts = df['Subtype'].value_counts()
        fig2, ax2 = plt.subplots()
        ax2.pie(counts, labels=counts.index, autopct='%1.1f%%', startangle=140)
        ax2.axis('equal')
        st.pyplot(fig2)

        st.subheader("Motif Counts (Bar Chart)")
        fig3, ax3 = plt.subplots()
        counts.plot.bar(ax=ax3)
        ax3.set_ylabel("Count")
        ax3.set_xlabel("Motif Type")
        plt.tight_layout()
        st.pyplot(fig3)

elif page == "Download":
    st.header("Download Motif Report")
    df = st.session_state.get('df', pd.DataFrame())
    if df.empty:
        st.info("No results to download. Run analysis first.")
    else:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download Results as CSV",
            data=csv,
            file_name=f"motif_results_{datetime.now().strftime('%Y%m%d-%H%M%S')}.csv",
            mime="text/csv"
        )
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False)
        st.download_button(
            label="Download Results as Excel",
            data=output.getvalue(),
            file_name=f"motif_results_{datetime.now().strftime('%Y%m%d-%H%M%S')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

elif page == "Additional Information":
    st.header("Additional Information")
    st.markdown("""
    - This app detects a broad array of non-B DNA motifs using reference algorithms.
    - For details, visit [GitHub](https://github.com/VRYella/Non-B-DNA-Finder).
    - Developed by Dr. Venkata Rajesh Yella & Chandrika Gummadi.
    """)
    st.markdown("---")
    st.subheader("How Are Motifs Predicted? (Technical, Stepwise Details)")
    st.markdown("""
<b>G-Quadruplex:</b><br>
Pattern: Four runs of G (≥3) separated by 1–7 bases.<br>
Steps: Regex scan. Each match scored by G4Hunter.<br>
<br>
<b>i-Motif:</b><br>
Pattern: Four runs of C (≥3) separated by 1–7 bases.<br>
Steps: Regex scan. Each match scored by G4Hunter (as negative value).<br>
<br>
<b>G-Triplex:</b><br>
Pattern: Three runs of G (≥3) separated by 1–7 bases, not matching G4 pattern.<br>
Steps: Regex scan, exclude G4s. Score = 0.7 × G4Hunter.<br>
<br>
<b>Bipartite G4:</b><br>
Pattern: Two G-quadruplexes within 100 nt.<br>
Steps: Regex for two G4s separated by ≤100 bases. Score: mean G4Hunter.<br>
<br>
<b>Multimeric G4:</b><br>
Pattern: Multiple tandem G4s (≥2) separated by ≤50 nt.<br>
Steps: Regex for multiple G4s. Score: mean G4Hunter.<br>
<br>
<b>Z-DNA:</b><br>
Pattern: Alternating purine/pyrimidine dinucleotides (>10 bp).<br>
Steps: Regex for (GC|CG|GT|TG|AC|CA) repeats. Score: Z-Seeker.<br>
<br>
<b>H-DNA:</b><br>
Pattern: Homopurine/homopyrimidine mirror repeats (≥10 nt) with ≤8 nt spacer.<br>
Steps: Regex for repeats, no score.<br>
<br>
<b>Sticky DNA:</b><br>
Pattern: ≥5 GAA or TTC repeats.<br>
Steps: Regex for runs, score = number of repeats.<br>
<br>
<b>Slipped DNA:</b><br>
Pattern: Direct repeats of 10–25 nt, separated by ≤10 nt.<br>
Steps: Regex for repeats, score = repeat region length.<br>
<br>
<b>Cruciform DNA:</b><br>
Pattern: Inverted repeats (arms ≥6 nt, loop ≤100 nt).<br>
Steps: For arm/loop sizes, scan for inverted repeats; verify reverse complement.<br>
Score: arm length.<br>
<br>
<b>Bent DNA:</b><br>
Pattern: ≥3 A-tracts (3–11 nt) with 3–11 nt spacers.<br>
Steps: Regex for periodic A-tracts, score = tract length.<br>
<br>
<b>APR:</b><br>
Pattern: Same as Bent DNA.<br>
Steps: As above.<br>
<br>
<b>Mirror Repeats:</b><br>
Pattern: Two arms (≥10 nt) flanking ≤100 nt loop, arms are mirrors.<br>
Steps: Regex for arms, verify left is reverse of right.<br>
<br>
<b>Quadruplex-Triplex Hybrid:</b><br>
Pattern: G4 and triplex within 100 nt.<br>
Steps: Regex for G4 plus triplex, no score.<br>
<br>
<b>Cruciform-Triplex Junction:</b><br>
Pattern: Cruciform arms with triplex nearby.<br>
Steps: Regex for pattern, no score.<br>
<br>
<b>G4-iMotif Hybrid:</b><br>
Pattern: G4 and i-motif within 100 nt.<br>
Steps: Find all G4s/i-motifs, report if within 100 nt.<br>
""", unsafe_allow_html=True)

st.markdown("""
---
**Developed by [Dr. Venkata Rajesh Yella] & Chandrika Gummadi** | [GitHub](https://github.com/VRYella/Non-B-DNA-Finder)
""")
