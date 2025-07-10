import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re, io
from datetime import datetime
from motifs import (
    all_motifs, find_hotspots
)
from utils import parse_fasta, wrap

EXAMPLE_FASTA = ">Example\nATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC\nATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT\nGGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA\nGAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"

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
    "Additional Information",
    "Motif Definitions Glossary"
]

st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", PAGES)

st.title("Non-B DNA Motif Finder")

def collect_all_motifs(seq, status_callback=None, stop_flag=None):
    if status_callback:
        status_callback("Scanning for non-B DNA motifs using non-overlapping regex patterns...")
    if stop_flag and stop_flag():
        return []
    return all_motifs(seq)

if page == "Home":
    st.markdown("""
    <style>
    .home-header {
        font-size: 2.5em;
        font-weight: bold;
        color: #1A5276;
        text-align: center;
        margin-bottom: 20px;
    }
    </style>
    <div class='home-header'>Welcome to Non-B DNA Motif Finder</div>
    """, unsafe_allow_html=True)
    try:
        st.image("nbd.PNG", use_container_width=True)
    except Exception:
        st.warning("Logo image (nbd.PNG) not found. Place it in the app folder.")
    st.markdown("""
    This application allows you to identify and visualize various non-B DNA motifs including:
    - Canonical and non-canonical G-quadruplexes
    - i-Motifs, Z-DNA, Cruciforms, and more
    
    Upload a sequence or use the example to begin.
    """)

elif page == "Upload & Analyze":
    st.markdown("<h2 style='color:#1A5276;'>Upload & Analyze</h2>", unsafe_allow_html=True)
    fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
    if fasta_file:
        seq = parse_fasta(fasta_file.read().decode("utf-8"))
        st.session_state['seq'] = seq
        st.success("FASTA loaded.")

    if st.button("Use Example Sequence"):
        st.session_state['seq'] = parse_fasta(EXAMPLE_FASTA)

    seq_input = st.text_area("Paste Sequence (FASTA or Raw)", value=st.session_state.get('seq', ''), height=150)
    if seq_input:
        try:
            st.session_state['seq'] = parse_fasta(seq_input)
        except Exception:
            st.error("Invalid FASTA or sequence.")

    if st.button("Run Motif Analysis"):
        seq = st.session_state.get('seq', '')
        if not seq or not re.match("^[ATGC]+$", seq):
            st.error("Please input a valid A/T/G/C sequence.")
        else:
            with st.spinner("Running motif analysis..."):
                results = collect_all_motifs(seq)
                st.session_state['motif_results'] = results
                st.session_state['df'] = pd.DataFrame(results)
            if results:
                st.success(f"Found {len(results)} motifs in sequence.")
            else:
                st.warning("No motifs found.")

elif page == "Results":
    st.markdown("<h2 style='color:#1A5276;'>Detected Motifs</h2>", unsafe_allow_html=True)
    df = st.session_state.get('df', pd.DataFrame())
    if df.empty:
        st.info("No results. Please run analysis first.")
    else:
        st.dataframe(df[['Class', 'Subtype', 'Start', 'End', 'Length', 'Sequence', 'ScoreMethod', 'Score']], use_container_width=True)

        st.subheader("Motif Type Distribution")
        counts = df['Subtype'].value_counts()
        fig, ax = plt.subplots()
        counts.plot(kind='bar', ax=ax)
        ax.set_xlabel("Motif Type")
        ax.set_ylabel("Count")
        st.pyplot(fig)

        st.subheader("Hotspot Regions (≥3 motifs in 100 nt)")
        hotspots = find_hotspots(st.session_state['seq'], st.session_state['motif_results'], window=100, min_count=3)
        if hotspots:
            st.dataframe(pd.DataFrame(hotspots))
        else:
            st.info("No hotspot regions found.")

elif page == "Visualization":
    st.markdown("<h2 style='color:#1A5276;'>Motif Visualization</h2>", unsafe_allow_html=True)
    df = st.session_state.get('df', pd.DataFrame())
    seq = st.session_state.get('seq', '')
    if df.empty:
        st.info("Run analysis to generate motifs first.")
    else:
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

elif page == "Download":
    st.markdown("<h2 style='color:#1A5276;'>Download Results</h2>", unsafe_allow_html=True)
    df = st.session_state.get('df', pd.DataFrame())
    if df.empty:
        st.info("No results to download. Run analysis first.")
    else:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("Download CSV", data=csv, file_name="motif_results.csv", mime="text/csv")

        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False)
        st.download_button("Download Excel", data=excel_buffer.getvalue(), file_name="motif_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

elif page == "Additional Information":
    st.markdown("<h2 style='color:#1A5276;'>Additional Information</h2>", unsafe_allow_html=True)
    st.markdown("""
This app detects a broad array of non-B DNA motifs using reference algorithms.\n
Visit [GitHub](https://github.com/VRYella/Non-B-DNA-Finder) for source code and documentation.\n
Developed by Dr. Venkata Rajesh Yella & Chandrika Gummadi.
    """)

elif page == "Motif Definitions Glossary":
    st.markdown("<h2 style='color:#1A5276;'>Motif Definitions Glossary</h2>", unsafe_allow_html=True)
    st.markdown("""
- **G-Quadruplex (G4):** Four runs of guanines forming a square planar structure stabilized by Hoogsteen hydrogen bonding.
- **i-Motif:** Cytosine-rich regions that form intercalated four-stranded structures in acidic pH.
- **Z-DNA:** Left-handed DNA helix often found in alternating purine-pyrimidine sequences.
- **Cruciform:** Hairpin-like structures in palindromic sequences.
- **Triplex DNA:** Three-stranded structures often formed by homopurine or homopyrimidine sequences.
- **Sticky DNA:** Repetitive purine/pyrimidine tracts that hybridize with themselves.
- **Bent DNA:** PolyA or PolyT tracts that locally bend the helix.
- **APR:** A-Phased Repeats associated with nucleosome positioning.
- **Mirror Repeats:** Symmetric sequences within the same DNA strand.
- **Hybrid Structures:** Junctions of G4, i-Motif, triplex, and cruciform motifs co-localized in sequence.
    """)
