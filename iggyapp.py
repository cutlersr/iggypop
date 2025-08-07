import sys
import os
import subprocess
import shutil
import streamlit as st
import yaml  

class rsc:
    @staticmethod
    def get_image(filename):
        return f"./images/{filename}"

def setup_page():
    st.set_page_config(
        page_title="Iggypop",
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
            'https://github.com/cutlersr/iggypop/tree/main. The protocol for the gene assembly pipeline is available at: '
            'https://www.protocols.io/view/iggypop-rapid-and-large-scale-dna-assembly-method-eq2lyqyzqvx9/v1 '
            'Please, do not forget to cite our  paper if you found this useful for your research. '
        )
        st.markdown(
            "<div style='font-size:10px; margin-top:50px;'>"
            "G. Dvir, Z. Xing, I. Beldman, A. Rivera, I. Wheeldon, & S.R. Cutler, "
            "Synthesis of large single-transcript pathways from oligonucleotide pools: Design of STARBURST, "
            "an autobioluminescent reporter. Proc. Natl. Acad. Sci. U.S.A. 122 (31) e2508109122,"
            "https://doi.org/10.1073/pnas.2508109122 (2025)."
            "</div>",
            unsafe_allow_html=True,
        )

def go_to(page_name):
    st.session_state.page = page_name

def collect_common_options(defaults=None):
    """
    Render widgets for all common iggypop parameters,
    seeding defaults from the `defaults` dict.
    """
    if defaults is None:
        defaults = {}
    cmd_args = []

    st.markdown("### Sequence optimization parameters")
    modes = ["chisel", "no_mods", "no_hinge"]
    default_mode = defaults.get("mode", modes[0])
    idx_mode = modes.index(default_mode) if default_mode in modes else 0
    st.subheader("Mode")
    mode = st.radio("", modes, index=idx_mode, key="mode_opt", help = "Operation mode. chisel will modify the input sequence. no_mods just hinges and barcodes. no_hinge will chisel without hinging")
    cmd_args += ["--mode", mode]

    codon_opts = ["use_best_codon", "match_codon_usage", "hybrid", "none"]
    default_cod = defaults.get("codon_opt", codon_opts[-1])
    idx_cod = codon_opts.index(default_cod) if default_cod in codon_opts else len(codon_opts)-1
    st.subheader("Codon optimization method")
    codon_method = st.radio("", codon_opts, index=idx_cod, key="codon_opt", help = "Select your codon optimizitation strategy.")
    cmd_args += ["--codon_opt", codon_method]
    if codon_method != "none":
        def_sp = defaults.get("species", "arabidopsis")
        st.subheader("Species")
        species = st.text_input("Species", value=def_sp, key="species")
        cmd_args += ["--species", species]
    if codon_method == "hybrid":
        def_pct = defaults.get("pct", 20.0)
        st.subheader("Percent divergence")
        pct = st.number_input("Percent divergence", value=float(def_pct), min_value=0.0, max_value=100.0, key="pct_div")
        cmd_args += ["--pct", str(pct)]

    st.markdown("### Cloning parameters")
    def_ext = defaults.get("ext_overhangs", ["AATG","GCTT"])
    st.subheader("External overhangs")
    over5 = st.text_input("External overhang 5′", value=def_ext[0], key="ext5", help = "The 5´ overhang that will pair with your plasmid.")
    over3 = st.text_input("External overhang 3′", value=def_ext[1], key="ext3", help = "The 3´ overhang that will pair with your plasmid.")
    cmd_args += ["--ext_overhangs", over5, over3]

    def_b5 = defaults.get("base_5p_end", "AATGCGGTCTCTA")
    st.subheader("Base 5′ end")
    base5 = st.text_input("", value=def_b5, key="base5", help = "Constant sequence added to the 5´ of the final designed sequences.")
    cmd_args += ["--base_5p_end", base5]

    def_b3 = defaults.get("base_3p_end", "GCTTAGAGACCGCTT")
    st.subheader("Base 3′ end")
    base3 = st.text_input("", value=def_b3, key="base3", help = "Constant sequence added to the 3´ of the final designed sequences.")
    cmd_args += ["--base_3p_end", base3]

    def_p5 = defaults.get("pcr_5p_cut", "CGTCTCA")
    st.subheader("PCR 5′ CUT")
    p5 = st.text_input("", value=def_p5, key="pcr5", help = "Sequence appended to the 5' of each oligo for assembly.")
    cmd_args += ["--pcr_5p_cut", p5]

    def_p3 = defaults.get("pcr_3p_cut", "AGAGACG")
    st.subheader("PCR 3′ CUT")
    p3 = st.text_input("", value=def_p3, key="pcr3", help = "Sequence appended to the 3´ of each oligo for assembly.")
    cmd_args += ["--pcr_3p_cut", p3]

    st.markdown("### Two-step parameters")
    default_two = defaults.get("two_step", "off").lower()
    idx_two = 0 if default_two=="on" else 1
    st.subheader("Two-step assemblies")
    two = st.radio("", ["On","Off"], index=idx_two, key="two_step", help = "Enable a two-step assembly workflow.")
    if two=="On":
        cmd_args += ["--two_step","on"]
        def_len = defaults.get("two_step_length",1200)
        st.subheader("Two-step length")
        ts_len = st.number_input("", value=int(def_len), min_value=1, key="two_step_len", help = "Max length for two-step fragments.")
        cmd_args += ["--two_step_length",str(ts_len)]
        def_5p = defaults.get("two_step_5p_end","AATGCGTCTCA")
        st.subheader("Two-step 5′ end")
        ts5 = st.text_input("", value=def_5p, key="two_step_5", help = "5' end for two-step cloning.")
        cmd_args += ["--two_step_5p_end",ts5]
        def_3p = defaults.get("two_step_3p_end","AGAGACGGCTT")
        st.subheader("Two-step 3′ end")
        ts3 = st.text_input("",value=def_3p,key="two_step_3", help = "3' end for two-step cloning.")
        cmd_args += ["--two_step_3p_end",ts3]

    st.markdown("### Hinging parameters")
    def_ol = defaults.get("oligo_length",250)
    st.subheader("Oligo length")
    ol = st.number_input("", value=int(def_ol), min_value=1, key="oligo_len", help = "Max allowed oligo length.")
    cmd_args += ["--oligo_length",str(ol)]
    def_rad = defaults.get("radius",8)
    st.subheader("Radius")
    rad = st.number_input("",value=int(def_rad),min_value=1,key="radius", help = "Allowable distance from ideal cut sites for selecting overhangs.")
    cmd_args += ["--radius",str(rad)]
    def_nt = defaults.get("n_tries",10)
    st.subheader("Number of tries")
    nt = st.number_input("",value=int(def_nt),min_value=1,key="n_tries", help = "Number of potential overhang sets to consider.")
    cmd_args += ["--n_tries",str(nt)]
    def_pi = defaults.get("primer_index",1)
    st.subheader("Primer index")
    pi = st.number_input("", value=int(def_pi), min_value=1, key="primer_idx", help = "The first primer index of your current run.")
    cmd_args += ["--primer_index",str(pi)]
    def_rep = defaults.get("repeats",1)
    st.subheader("Repeats")
    rep = st.number_input("",value=int(def_rep),min_value=1,key="repeats", help = "Number of chiseled sequences to create per input sequence.")
    cmd_args += ["--repeats",str(rep)]
    def_mf = defaults.get("max_fragments",18)
    #st.subheader("Maximum fragments per PCR")
    #mf = st.number_input("",value=int(def_mf),min_value=1,key="max_frags", help = "Max number of fragments per PCR.")
    cmd_args += ["--max_fragments",str(def_mf)]

    def_repflag = defaults.get("no_reports",False)
    st.subheader("Reports")
    rflags = st.radio("",["On","Off"],index=(1 if def_repflag else 0),key="reports", help = "Enable creation of PDF DNA chisel reports.")
    if rflags=="Off":
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
            st.title("CDS Runtype")
            uploaded = st.file_uploader("Upload your FASTA file", type=["fasta","fa","txt"])
            if uploaded:
                st.write(f"Uploaded file: **{uploaded.name}**")
            with st.expander("Advanced mode"):
                common_args = collect_common_options()
            col1, col2, col3 = st.columns(3)
            pop = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=("main",))
            def to_yaml_cds():
                if not uploaded:
                    st.error("Please upload a FASTA file first.")
                    return
                upl="./uploads"; os.makedirs(upl,exist_ok=True)
                ip=os.path.join(upl,uploaded.name)
                with open(ip,"wb") as f: f.write(uploaded.getbuffer())
                st.session_state.yaml_type="cds"
                st.session_state.input_path=ip
                st.session_state.last_stem=os.path.splitext(uploaded.name)[0]
                go_to("yaml")
            col3.button("YAML", on_click=to_yaml_cds)
            if pop:
                if not uploaded:
                    st.error("Please upload a FASTA file first.")
                else:
                    upl="./uploads";os.makedirs(upl,exist_ok=True)
                    ip=os.path.join(upl,uploaded.name)
                    with open(ip,"wb") as f: f.write(uploaded.getbuffer())
                    stem=os.path.splitext(uploaded.name)[0]
                    outd=os.path.join("out",stem)
                    if os.path.isdir(outd): shutil.rmtree(outd)
                    cmd=["./iggypop.py","cds","--i",ip,"--o",stem] + common_args
                    res=subprocess.run(cmd,capture_output=True,text=True)
                    if res.returncode==0:
                        zp=shutil.make_archive(outd,"zip",outd)
                        st.session_state.result_zip=zp
                        st.session_state.last_stem=stem
                        st.session_state.last_page="cds"
                        st.session_state.page="results"
                        st.rerun()
                    else:
                        st.error(f"iggypop failed (exit {res.returncode}):\n{res.stderr}")

        elif st.session_state.page == "gb":
            st.title("GB Runtype")
            st.markdown("This runtype will format your Genbank file before processing it. "
                        "Please check that the contents of your formatted file are correct before ordering "
                        "an oligo pool. You can do so in the assets directory inside your zip folder.")
            uploaded = st.file_uploader("Upload your GenBank file", type=["gb","gbk"])
            if uploaded:
                st.write(f"Uploaded file: **{uploaded.name}**")
            with st.expander("Advanced mode"):
                common_args = collect_common_options()
            col1, col2, col3 = st.columns(3)
            pop = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=("main",))
            def to_yaml_gb():
                if not uploaded:
                    st.error("Please upload a GenBank file first.")
                    return
                upl="./uploads";os.makedirs(upl,exist_ok=True)
                ip=os.path.join(upl,uploaded.name)
                with open(ip,"wb") as f: f.write(uploaded.getbuffer())
                stem=os.path.splitext(uploaded.name)[0]
                fmt=os.path.join(upl,f"{stem}_formatted.gb")
                r=subprocess.run(["./iggypop.py","format","--i",ip,"--o",fmt],capture_output=True,text=True)
                if r.returncode!=0:
                    st.error(f"Format failed:\n{r.stderr}")
                    return
                st.session_state.yaml_type="gb"
                st.session_state.input_path=fmt
                st.session_state.last_stem=stem
                go_to("yaml")
            col3.button("YAML", on_click=to_yaml_gb)
            if pop:
                if not uploaded:
                    st.error("Please upload a GenBank file first.")
                else:
                    upl="./uploads";os.makedirs(upl,exist_ok=True)
                    ip=os.path.join(upl,uploaded.name)
                    with open(ip,"wb") as f: f.write(uploaded.getbuffer())
                    stem=os.path.splitext(uploaded.name)[0]
                    outd=os.path.join("out",stem)
                    if os.path.isdir(outd): shutil.rmtree(outd)
                    fmt=os.path.join(upl,f"{stem}_formatted.gb")
                    r=subprocess.run(["./iggypop.py","format","--i",ip,"--o",fmt],capture_output=True,text=True)
                    if r.returncode!=0:
                        st.error(f"Format failed:\n{r.stderr}")
                    else:
                        cmd=["./iggypop.py","gb","--i",fmt,"--o",stem] + common_args
                        res=subprocess.run(cmd,capture_output=True,text=True)
                        if res.returncode==0:
                            zp=shutil.make_archive(outd,"zip",outd)
                            st.session_state.result_zip=zp
                            st.session_state.last_stem=stem
                            st.session_state.last_page="gb"
                            st.session_state.page="results"
                            st.rerun()
                        else:
                            st.error(f"iggypop gb failed:\n{res.stderr}")

        elif st.session_state.page == "results":
            st.title("Run Complete!")
            stem = st.session_state.last_stem
            zp = st.session_state.result_zip
            st.write(f"Your results for **{stem}** are ready:")
            with open(zp,"rb") as fp:
                st.download_button("Download Results", data=fp,
                                   file_name=os.path.basename(zp),
                                   mime="application/zip")
            last = st.session_state.get("last_page","cds")
            st.button("← Back to CDS" if last=="cds" else "← Back to GB",
                      on_click=go_to, args=(last,))

        elif st.session_state.page == "yaml":
            st.title("Select a YAML template")
            yaml_type = st.session_state.yaml_type

            # dynamic listing
            yaml_dir = "./yaml"
            try:
                files = sorted(os.listdir(yaml_dir))
            except FileNotFoundError:
                st.error(f"Directory '{yaml_dir}' not found.")
                files = []
            opts = [f for f in files if f.endswith(".yml")]
            if not opts:
                st.error("No YAML templates found.")
                return

            choice = st.selectbox("Template:", opts, key="yaml_choice")

            # Confirm loads defaults
            if st.button("Confirm"):
                path = os.path.join(yaml_dir, choice)
                with open(path) as yf:
                    parsed = yaml.safe_load(yf) or {}
                st.session_state.yaml_defaults = parsed
                st.session_state.yaml_args = None
                st.success("Loaded defaults from YAML.")

            # Advanced mode: collect args using those defaults
            with st.expander("Advanced mode"):
                defaults = st.session_state.get("yaml_defaults", {})
                st.session_state.yaml_args = collect_common_options(defaults=defaults)

            # Retrieve collected args
            common_args = st.session_state.get("yaml_args", [])

            col1, col2 = st.columns(2)
            pop_yaml = col1.button("Pop!")
            col2.button("← Back", on_click=go_to, args=(yaml_type,))

            if pop_yaml:
                ip = st.session_state.input_path
                stem = st.session_state.last_stem
                outd = os.path.join("out", stem)
                if os.path.isdir(outd):
                    shutil.rmtree(outd)
                cmd = ["./iggypop.py", yaml_type, "--i", ip, "--o", stem] \
                      + common_args \
                      + ["--yml", os.path.join("yaml", choice)]
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    zp = shutil.make_archive(outd, "zip", outd)
                    st.session_state.result_zip = zp
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

