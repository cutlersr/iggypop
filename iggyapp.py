import sys
import os
import subprocess
import shutil
import streamlit as st

class rsc:
    @staticmethod
    def get_image(filename):
        return f"./png/{filename}"

def setup_page():
    st.set_page_config(
        page_title="iggypop",
        page_icon=rsc.get_image("best_icon.png"),
        layout="centered",
        initial_sidebar_state="auto",
    )
    sys.tracebacklimit = 0
    dark_theme_style = """
    <style>
    </style>
    """
    st.markdown(dark_theme_style, unsafe_allow_html=True)

    with st.sidebar:
        st.markdown('### What is iggypop?')
        st.markdown(
            'iggypop is a pipeline for designing and synthesizing genes from oligonucleotide pools. '
            'Input sequences are fragmented into segments that can be amplified using gene-specific primers '
            'and reassembled by Golden Gate cloning. This website is a simplified, user-friendly version '
            'of the pipeline. To gain access to the full capabilities of the package visit our GitHub repository at: '
            'https://github.com/cutlersr/iggypop/tree/main. Please, do not forget to cite the iggypop paper if you found '
            'it useful for your research.'
        )
        st.markdown(
            "<div style='font-size:10px; margin-top:50px;'>"
            "Dvir, G., Xing, Z., Beldman, I., Rivera, A., Wheeldon, I., & Cutler, S. R. (2025). "
            "Synthesis of large single-transcript pathways from oligonucleotide pools: Design of STARBURST, "
            "an autobioluminescent reporter. Manuscript submitted for publication to PNAS (April 2025)."
            "</div>",
            unsafe_allow_html=True,
        )

def go_to(page_name):
    st.session_state.page = page_name

# The user can introduce custom arguments for iggypop on both runtypes
def collect_common_options():
    cmd_args = []

    st.subheader("Mode")
    mode = st.radio("", ["chisel", "no_mods", "no_hinge"], index=0)
    cmd_args += ["--mode", mode]

    st.subheader("Codon optimization method")
    codon_method = st.radio(
        "",
        ["use_best_codon", "match_codon_usage", "harmonize_rca", "hybrid", "none"],
        index=4,
    )
    cmd_args += ["--codon_opt", codon_method]

    if codon_method != "none":
        st.subheader("Species")
        species = st.text_input("Species", value="arabidopsis")
        cmd_args += ["--species", species]
        if codon_method == "harmonize_rca":
            st.subheader("Original species")
            orig_species = st.text_input("Original species", value="arabidopsis")
            cmd_args += ["--original_species", orig_species]

    if codon_method == "hybrid":
        st.subheader("Target sequence divergence for hybrid optimization")
        pct = st.number_input("Percent divergence", value=20.0, min_value=0.0, max_value=100.0)
        cmd_args += ["--pct", str(pct)]

    st.subheader("Repeats")
    repeats = st.number_input("", min_value=1, value=1)
    cmd_args += ["--repeats", str(repeats)]

    st.subheader("Oligo length")
    oligo_length = st.number_input("", min_value=1, value=250)
    cmd_args += ["--oligo_length", str(oligo_length)]

    st.subheader("External overhangs")
    over5 = st.text_input("External overhang 5′", value="AATG")
    over3 = st.text_input("External overhang 3′", value="GCTT")
    cmd_args += ["--ext_overhangs", over5, over3]

    st.subheader("Base 5′ end")
    base5 = st.text_input("Base 5′ end", value="AATGCGGTCTCTA")
    st.subheader("Base 3′ end")
    base3 = st.text_input("Base 3′ end", value="GCTTAGAGACCGCTT")
    cmd_args += ["--base_5p_end", base5, "--base_3p_end", base3]

    st.subheader("Two-step assemblies")
    two_step = st.radio("", ["On", "Off"], index=1)
    if two_step == "On":
        cmd_args += ["--two_step", "on"]

    st.subheader("PCR 5′ CUT")
    pcr5 = st.text_input("PCR 5′ CUT", value="CGTCTCA")
    st.subheader("PCR 3′ CUT")
    pcr3 = st.text_input("PCR 3′ CUT", value="AGAGACG")
    cmd_args += ["--pcr_5p_cut", pcr5, "--pcr_3p_cut", pcr3]

    st.subheader("Primer index")
    primer_idx = st.number_input(
        "This is the row of indexsets file to start from", min_value=1, value=1
    )
    cmd_args += ["--primer_index", str(primer_idx)]

    st.subheader("Number of tries")
    n_tries = st.number_input("", min_value=1, value=50)
    cmd_args += ["--n_tries", str(n_tries)]

    st.subheader("Radius")
    radius = st.number_input("", min_value=1, value=8)
    cmd_args += ["--radius", str(radius)]

    st.subheader("Maximum fragments per PCR")
    max_frags = st.number_input("", min_value=1, value=18)
    cmd_args += ["--max_fragments", str(max_frags)]

    st.subheader("Reports")
    reports = st.radio("",  ["On", "Off"], index = 0, key = "report")
    if reports == "Off":
        cmd_args += ["--no-reports"]

    return cmd_args

def run_streamlit():
    try:
        setup_page()

        if "page" not in st.session_state:
            st.session_state.page = "main"

        # MAIN MENU
        if st.session_state.page == "main":
            st.title("Welcome to Iggypop")
            st.write("Choose a runtype to continue:")
            c1, c2 = st.columns(2)
            c1.button("CDS", on_click=go_to, args=("cds",))
            c2.button("GB", on_click=go_to, args=("gb",))

        # CDS PAGE
        elif st.session_state.page == "cds":
            st.title("CDS")
            uploaded_fasta = st.file_uploader(
                "Upload your FASTA file", type=["fasta", "fa", "txt"]
            )
            if uploaded_fasta:
                st.write(f"Uploaded file: **{uploaded_fasta.name}**")

            with st.expander("Advanced mode"):
                common_args = collect_common_options()

            col1, col2, col3 = st.columns(3)
            pop_clicked = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=("main",))

            def _go_to_yaml_cds():
                if not uploaded_fasta:
                    st.error("Please upload a FASTA file first.")
                    return
                # save upload for YAML run
                upload_dir = "./uploads"
                os.makedirs(upload_dir, exist_ok=True)
                input_path = os.path.join(upload_dir, uploaded_fasta.name)
                with open(input_path, "wb") as f:
                    f.write(uploaded_fasta.getbuffer())
                stem = os.path.splitext(uploaded_fasta.name)[0]
                st.session_state.yaml_type = "cds"
                st.session_state.input_path = input_path
                st.session_state.last_stem = stem
                go_to("yaml")

            col3.button("YAML", on_click=_go_to_yaml_cds)

            if pop_clicked:
                if not uploaded_fasta:
                    st.error("Please upload a FASTA file first.")
                else:
                    upload_dir = "./uploads"
                    os.makedirs(upload_dir, exist_ok=True)
                    input_path = os.path.join(upload_dir, uploaded_fasta.name)
                    with open(input_path, "wb") as f:
                        f.write(uploaded_fasta.getbuffer())
                    stem = os.path.splitext(uploaded_fasta.name)[0]
                    out_dir = os.path.join("out", stem)
                    if os.path.isdir(out_dir):
                        shutil.rmtree(out_dir)
                    cmd = ["./iggypop.py", "cds", "--i", input_path, "--o", stem] + common_args
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    if result.returncode == 0:
                        zip_path = shutil.make_archive(out_dir, "zip", out_dir)
                        st.session_state.result_zip = zip_path
                        st.session_state.last_stem = stem
                        st.session_state.last_page = "cds"
                        st.session_state.page = "results"
                        st.rerun()
                    else:
                        st.error(f"iggypop failed (exit {result.returncode}):\n```\n{result.stderr}\n```")

        # GB PAGE
        elif st.session_state.page == "gb":
            st.title("GB")
            uploaded_gb = st.file_uploader(
                "Upload your GenBank file", type=["gb", "gbk"]
            )
            if uploaded_gb:
                st.write(f"Uploaded file: **{uploaded_gb.name}**")

            with st.expander("Advanced mode"):
                common_args = collect_common_options()

            col1, col2, col3 = st.columns(3)
            pop_clicked = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=("main",))

            def _go_to_yaml_gb():
                if not uploaded_gb:
                    st.error("Please upload a GenBank file first.")
                    return
                upload_dir = "./uploads"
                os.makedirs(upload_dir, exist_ok=True)
                input_path = os.path.join(upload_dir, uploaded_gb.name)
                with open(input_path, "wb") as f:
                    f.write(uploaded_gb.getbuffer())
                stem = os.path.splitext(uploaded_gb.name)[0]
                formatted_path = os.path.join(upload_dir, f"{stem}_formatted.gb")
                fmt_cmd = ["./iggypop.py", "format", "--i", input_path, "--o", formatted_path]
                fmt_res = subprocess.run(fmt_cmd, capture_output=True, text=True)
                if fmt_res.returncode != 0:
                    st.error(f"Format failed:\n{fmt_res.stderr}")
                    return
                st.session_state.yaml_type = "gb"
                st.session_state.input_path = formatted_path
                st.session_state.last_stem = stem
                go_to("yaml")

            col3.button("YAML", on_click=_go_to_yaml_gb)

            if pop_clicked:
                if not uploaded_gb:
                    st.error("Please upload a GenBank file first.")
                else:
                    upload_dir = "./uploads"
                    os.makedirs(upload_dir, exist_ok=True)
                    input_path = os.path.join(upload_dir, uploaded_gb.name)
                    with open(input_path, "wb") as f:
                        f.write(uploaded_gb.getbuffer())
                    stem = os.path.splitext(uploaded_gb.name)[0]
                    out_dir = os.path.join("out", stem)
                    if os.path.isdir(out_dir):
                        shutil.rmtree(out_dir)
                    formatted_path = os.path.join(upload_dir, f"{stem}_formatted.gb")
                    fmt_cmd = ["./iggypop.py", "format", "--i", input_path, "--o", formatted_path]
                    fmt_res = subprocess.run(fmt_cmd, capture_output=True, text=True)
                    if fmt_res.returncode != 0:
                        st.error(f"Format failed:\n{fmt_res.stderr}")
                    else:
                        cmd = ["./iggypop.py", "gb", "--i", formatted_path, "--o", stem] + common_args
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        if result.returncode == 0:
                            zip_path = shutil.make_archive(out_dir, "zip", out_dir)
                            st.session_state.result_zip = zip_path
                            st.session_state.last_stem = stem
                            st.session_state.last_page = "gb"
                            st.session_state.page = "results"
                            st.rerun()
                        else:
                            st.error(f"iggypop gb failed:\n{result.stderr}")

        # RESULTS PAGE
        elif st.session_state.page == "results":
            st.title("Run Complete!")
            stem = st.session_state.last_stem
            zip_path = st.session_state.result_zip
            st.write(f"Your results for **{stem}** are ready:")
            with open(zip_path, "rb") as fp:
                st.download_button("Download Results", data=fp,
                                   file_name=os.path.basename(zip_path),
                                   mime="application/zip")
            last = st.session_state.get("last_page", "cds")
            st.button("← Back to CDS" if last == "cds" else "← Back to GB",
                      on_click=go_to, args=(last,))

        # YAML PAGE
        elif st.session_state.page == "yaml":
            st.title("Select a YAML template")
            yaml_type = st.session_state.yaml_type
            if yaml_type == "cds":
                options = [
                    "domesticate_cds_hybrid.yml",
                    "domesticate_cds_mcu_gc_53.yml",
                    "domesticate_cds_mcu.yml",
                    "domesticate_cds_minimal.yml",
                    "domesticate_cds_ubc.yml",
                    "domesticate_cds.yml",
                    "domesticate_two_step_cds_hybrid.yml",
                    "domesticate_two_step_cds_mcu.yml",
                    "domesticate_two_step_cds_ubc.yml",
                    "domesticate_two_step_cds.yml",
                ]
            else:
                options = [
                    "domesticate_gb_hybrid.yml",
                    "domesticate_gb_mcu.yml",
                    "domesticate_gb_ubc.yml",
                    "domesticate_gb.yml",
                    "domesticate_two_step_gb_hybrid.yml",
                    "domesticate_two_step_gb_mcu.yml",
                    "domesticate_two_step_gb_ubc.yml",
                    "domesticate_two_step_gb.yml",
                ]
            choice = st.selectbox("Template:", options)
            col1, col2 = st.columns(2)
            pop_yaml = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=(yaml_type,))

            if pop_yaml:
                input_path = st.session_state.input_path
                stem = st.session_state.last_stem
                out_dir = os.path.join("out", stem)
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)

                cmd = [
                    "./iggypop.py", yaml_type,
                    "--i", input_path,
                    "--o", stem,
                    "--yml", f"yaml/{choice}"
                ]
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    zip_path = shutil.make_archive(out_dir, "zip", out_dir)
                    st.session_state.result_zip = zip_path
                    st.session_state.last_page = yaml_type
                    st.session_state.page = "results"
                    st.rerun()
                else:
                    st.error(f"iggypop failed:\n{result.stderr}")

    except Exception as e:
        st.session_state.page = "error"
        st.session_state.last_error = str(e)
        st.rerun()

    # ERROR PAGE
    if st.session_state.page == "error":
        st.title("An error occurred")
        st.error(f"**Error details:**\n\n```\n{st.session_state.last_error}\n```")
        st.button("← Back to main", on_click=go_to, args=("main",))

if __name__ == "__main__":
    run_streamlit()

