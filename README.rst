==========================
iggypop
==========================

**indexed golden gate assembly by PCRing barcoded oligo pools**
==========================


iggypop is a pipeline for creating synthetic genes at $3.00 - $7.00 per Kb in oligo costs. It uses the Edinburgh Genome Foundry's `dnachisel <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ to optimize sequences and `goldenhinges <https://github.com/Edinburgh-Genome-Foundry/GoldenHinges>`_ to fragment them into barcoded pieces that can be reassembled by golden gate assembly. High-fidelity assembly overhangs are rapidly selected from pre-computed sets, and the fragmented genes are amplified from oligo pools using experimentally validated barcode primer pairs.

iggypop enables end-to-end design, assembly, and validation of 100s of genetic parts in a single experiment.

.. image:: png/overview.png
   :alt: Overview



Installation
============



**Linux**

.. code:: bash

    git clone github.com/cutlersr/iggypop
    cd iggypop
    python -m venv .venv
    source .venv/bin/activate
    chmod +x setup.sh
    sudo ./setup.sh



**Docker**

.. code:: bash

    git clone github.com/cutlersr/iggypop
    cd iggypop
    docker build -t iggypop .
    docker run -v $(pwd):/app -it iggypop /bin/bash



Usage
=====

Coding Sequence Mode
---------------------



CDS mode is used to design oligos for reassembling coding sequences. It requires fasta formatted ORFs. The default cds mode `yaml <yaml/moclo_cds_mcu.yml>`_ design MoClo compatible ORFs that lack common gg sites (BsaI, BsmBI, BbsI, SapI, BtgZI), match Arabidopsis codon usage, minimize micro-homologies (10 bp repeats) and hairpins, have G/C content ≤ 0.60.

.. code:: bash

    # CDS in FASTA format
    ./iggypop.py cds --i in/test.fasta --o test_oligos



Sequence optimization parameters are set in YAML files using `dnachisel` `specifications <https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html>`_. The following YAML snippet sets dnachisel constraints to exclude BsaI sites, prevent changes to the first 200 base pairs of the sequence, and enforce synonymous changes:

.. code:: yaml

    constraints
          - type: AvoidPattern          # Pattern/RE site to avoid
            pattern: BsaI_site
            strand: both    
          - type: AvoidChanges          # Prevent changes in the first 200 base pairs
            location: [1, 200]
          - type: EnforceTranslation    # Force changes to be synonymous


You can also change settings on the command line, which overrides the parameters in the YAML file.

.. code:: bash

    ./iggypop.py cds  --i in/cds_test.fasta  --species s_cerevisiae    \    # use S. cerevisiae codon table
                      --pcr_5p_cut GGTCTCT  --pcr_3p_cut AGAGACC       \    # use BsaI ends instead of default BsmBI
                      --codon_opt use_best_codon  --oligo_length 300        # ubc w/ 300 bp oligos



Genbank File Mode
------------------

The parameters for optimizing GenBank files are set with annotations according to `dnachisel's genbank API <https://edinburgh-genome-foundry.github.io/DnaChisel/genbank/genbank_api.html>`_. `./iggypop.py format` is used to annotate the gb file using parameters specified in a YAML; once formatted, oligos can are then generated with `./iggypop.py gb`

.. code:: bash

    # Format a Genbank file 
    
    ./iggypop.py format --i in/test_unformatted.gb    \
                        --o in/test_formatted.gb      \
       # optionally override the defaults             \
                        --species b_subtilis          \
                        --codon_opt use_best_codon    
                        
    
    # Then run the formatted Genbank file to generate oligos
    ./iggypop.py gb  --i in/test_formatted.gb --o test_oligos


We recommend you check the formatting produced by `./iggypop.py format` in Snapgene, Geneious, Benchling, or your favorite viewer.



Design
=====

The `yaml/ <yaml/>`_ folder contains parameter files for some common design strategies. 

 


MoClo-compatible CDSs
-----------------------

The `moclo` YAML files have parameters to design reusable CDSs by adding a 5' BsaI/ATTG on the 5' end and a 3' GCTT/BsaI. The CDSs are first assembled with BsmBI, and the final clones are compatible with MoClo (i.e., BsaI digestion will release an ORF with AATG/GCTT overhangs). The figure below zooms in on the first and last oligonucleotides of an assembly to illustrate the mapping between parameters set in the YAML file and the designed sequence. This is the default design mode; see the YAML folder for other options. 

.. image:: png/molco.png
   :alt: MoClo Compatibility



Two-step assembly
-------------------

For target sequences longer than 3 Kb (~16-18 fragments encoded in 250 bp oligos), the frequency of proper assemblies is low enough that it can be more efficient to break the target sequences into smaller step one fragments that are cloned, sequence validated and then used for second step assemblies to yield the final target. The `two_step` YAML files have parameters to break a sequence into ~ 1 Kb chunks assembled with BbsI; once validated, they are assembled into the final sequence using BsmBI in the second step. The figure below shows the first and last oligos of a 2-step assembly. The cds mode defaults are for MoClo-compatibility; the gb versions not.

.. image:: png/two_step.png
   :alt: Two-step Assembly

.. code:: bash

    ./iggypop.py gb --i in/test.gb --two_step on  --o two_step



Versioning
---------------------

You may want to test different versions of the same coding sequence, as there can be substantial variation in expression between codon-optimized variants of the same amino acid sequence. The `--repeats` parameter allows you to generate multiple versions.

.. code:: bash

    ./iggypop.py gb --i in/test.gb --repeats 5 --o repeats



Deintronization
-----------------

Sequences ported from other organisms or newly designed sequences sometimes contain cryptic introns that reduce or kill expression in a eukaryotic host. With `./iggypop.py cds --deintronize on`, a chiseled CDS is generated, passed to a splicing model from the `Spliceator <https://link.springer.com/article/10.1007/s00438-016-1258-6>`_ project. Potential intron donor and acceptor sites are identified, if any, and fed back to `dnachisel` and eliminated using `@AvoidPattern`. The chiseled sequence is reanalyzed and continues up to 5 times or until a deintronized CDS is identified.

.. code:: bash

    ./iggypop.py cds  --i in/test.fasta --deintronize on --o deintronized



`Hybrid` codon optimization
-----------------

The two main methods of optimizing sequences are match_codon_usage (MCU), which randomly samples codons based on their usage frequency, and use_best_codon (UBC). MCU generates sequences that typically have `CAI <https://en.wikipedia.org/wiki/Codon_Adaptation_Index>`_ values of ~0.75, and UBC generates CAI values 1. In some cases, you may want CAI values between those ranges, for example, to create many versions of high CAI sequences (UBC usually generates only one sequence). The --codon_opt  hybrid parameter allows this with the `--pct` parameter determining the target sequence difference from the input sequence (the default values shoot for ~20% difference). You may need to tweak the pct parameter to hit your desired CAI value.

.. code:: bash

    ./iggypop.py cds --i in/test.fasta --codon_opt hybrid --pct 30 --o hybrid



Codon tables
=====

For cds mode, a condensed local version of the `cocoputs <https://pubmed.ncbi.nlm.nih.gov/31029701>`_ database is used for codon table lookups (but this can be changed using the `--codon_tbl kazusa` flag). Gb mode uses Kazusa codon tables. Based on our lab's most common use cases, **cds mode defaults to an arabidopsis codon table, and gb mode defaults to an *E.coli* codon table**. To change this, use the `--species flag`. To optimize for the monkeyflower *Erythranthe guttata*, these all work:

.. code:: bash

    # short name, species name, or taxid work cds mode
    ./iggypop.py cds  --i in/test.fasta --species e_guttata
    ./iggypop.py cds  --i in/test.fasta --species Erythranthe guttata   
    ./iggypop.py cds  --i in/test.fasta --species 4155


Reports & quiet
-----------------

You can generate a dnachisel report with `--reports`; if you want iggypop to print less to the screen use `--quiet`



Vectors
=======

We've developed a series of pPOP vectors for the one-step and two-step cloning modes; they are derivatives of pUPD2 and pCAMBIA. Sequences can be found `here <vectors/>`_.


Barcode primers
----------------

Our barcode primers were designed to have balanced Tms, lack commonly used restriction sites, not dimerize, and be as small as possible (to maximize the sequence per oligo dedicated to the target). In addition, we wanted to minimize off-target hybridization between barcodes and eliminate amplification of common contaminating DNAs (*E. coli*, T7, T4, others). To accomplish this, a large set of 18 bp primers was generated with `primer3`; `MFEprimer3 <https://academic.oup.com/nar/article/47/W1/W610/5486745>`_ was then used to minimize barcodes that cross-hybridized with one another or  were predicted to produce amplicons from contaminant DNA. 350 pairs of primers designed with this pipeline have been experimentally validated; 4 of 350 designs failed. We've purged the pairs we know are defective and will update the primer file as we get more validation. 

You probably don't need to start from scratch, but if you do, here's the pipeline:

.. code:: bash

    # example run settings; set `num_sequences` to more than you need to account for # bad primers eliminated after the MFEprimer steps.
    ./iggypop.py primers  --num_sequences 10    \
    --opt_tm 60  --opt_size 18 --gc_content 0.5 \
    --max_size 18 --max_size 18


Overhangs
-------------

We use `goldenhinges` to select fragment junctions and assembly overhangs. Given an input sequence, target fragment length, and list of candidate overhangs, `goldenhinges` searches for solutions within a given radius of the ideal target cut sites. To produce high-fidelity assemblies, we limit `goldenhinges` overhang choices to pre-computed high-fidelity overhang sets; `iggypop` searches through these sets to identify `n_tries` solutions, and returns the highest fidelity set obtained. The data below show the fidelities obtained for a run of 4,500 plant transcription factors using AATG/GCTT cloning overhang recombination with our overhang sets; in this run, the mean fragment number is 7 (~1.2 kB), and the mean assembly fidelity is predicted, to be 99.5%.

.. image:: png/fidelity_plot.png
   :alt: fidelity_plot

.

The overhang sets we use (`hf_oh_sets.xlsx`) were generated using both genetic algorithm and Monte Carlo optimizers. The sets were optimized with `AATG, GCTT` as the `fixed_overhangs` (i.e., external cloning overhangs); AATG and GCTT have near-perfect fidelity and are MoClo-compliant for CDSs, so it's easy to create high-fidelity sets using them. Fidelities are calculated using `Potapov et al. <https://pubs.acs.org/doi/10.1021/acssynbio.8b00333>`_ data for one-hour incubations at 25 ºC using T4 DNA ligase; you can change this with the `potapov_data` setting. You can specify whatever external overhangs you want but check with NEB's `fidelity calculator <https://ligasefidelity.neb.com/viewset/run.cgi>`_ to ensure they are a high-fidelity pair first. The following command will do a run with a target of a set of 20 overhangs. Due to how GAs work, sets with repeated sequences can arise; the `alpha` and `beta` parameters below control a penalty function that reduces repeated overhangs. For the overhang sets used, we ran a few thousand gaga runs on UCR's high-performance computing cluster and filtered the results to select the highest-scoring sets and maximally diverse subsets.


.. code:: bash

    # run a bunch of optimizations
    /iggypop.py gagga 
        --set_size=20            \
        --ngen=150               \
        --pop_size=1000          \
        --min_improve=.0005      \
        --alpha 2.4              \
        --beta 2.4               \
        --tournament_size 4 
    
    # then run this from the directory with all of your results
    Rscript scripts/process_gagga_runs.R --top_percent=2 --n_cliques=30

