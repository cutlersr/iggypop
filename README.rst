iggypop
========

Indexed Golden Gate Assembly of Fragments PCR Amplified from Oligo Pools
-------------------------------------------
.. image:: png/overview.png

**iggypop** is a pipeline for creating synthetic genes at $3.00 - $7.00 per kB in oligo costs. It uses the Edinburgh Genome Foundry's `dnachisel <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ to optimize sequences and `goldenhinges <https://github.com/Edinburgh-Genome-Foundry/GoldenHinges>`_ to fragment them into barcoded oligos that can be reassembled by Golden Gate cloning. The overhangs are selected using pre-computed, high-fidelity hingesets, and the fragmented genes are amplified from oligo pools using experimentally validated indexsets primers. Reactions are assembled and cloned into pPOP vectors and then nanopore sequenced using barcoded amplicons and the `iggyseq` pipeline.

Installation
------------

**Linux**

.. code-block:: bash

    # requires python 3.8
    git clone github.com/cutlersr/iggypop
    cd iggypop
    python3.8 -m venv .venv
    source .venv/bin/activate
    chmod +x setup.sh
    ./setup.sh

**Docker**

.. code-block:: bash

    git clone github.com/cutlersr/iggypop
    cd iggypop
    docker build -t iggypop .
    docker run -it -v $(pwd):/app iggypop

Fragmenting and Optimizing Coding Sequences
-------------------------------------------

Coding sequences are domesticated, fragmented, indexed, and appended with cut sites to yield oligonucleotides that can be amplified with gene-specific primers and then assembled. Sequence domestication and optimization is conducted using `dnachisel`; optimization parameters can be set in a yaml file using dnachisel's built-in `specifications <https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html>`_; several yaml files used in our common workflows are in the `/yaml` folder.

.. code-block:: bash

    ./iggypop.py cds --i "in/10_TFs.fasta" --o "10_TFs"


The default paramaters can be modified by creating new yaml files or from the command line. To codon optimize with an *E. coli* codon table, use BsaI sites for assembly, and synthesize 300 bp oligos:

.. code-block:: bash

    ./iggypop.py cds                                        \
        --i "test/10_TFs.fasta" --o "10_TFs_coli_mcu"       \
        --base_5p_end "GGTCTCA" --base_3p_end "AGAGACC"     \ # BsaI instead of BsmBI
        --codon_opt "match_codon_usage" --species "e_coli"  \
        --oligo_length 300  # default is 250

The default `iggypop.py cds` settings design ORFs that:

- Lack common Golden Gate cloning sites (BsaI, BsmBI, and BbsI)
- Enforce synonymous changes
- Assemble from oligos ≤ 250 bp with BsmBI
- Are GoldenBraid / MoClo compatible (inner 5'-BsaI-AATG...GCTT-BsaI-3')

Modify as needed.

Genbank File Mode
-----------------

The parameters for optimizing GenBank files differ and use annotations added to your GenBank file using `dnachisel's genbank API <https://edinburgh-genome-foundry.github.io/DnaChisel/genbank/genbank_api.html>`_. `iggypop format` allows easy parameter setting in a yaml file:

.. code-block:: bash

    # Format a genbank file using the default domesticate_gb.yml file
    ./iggypop.py format --i "in/sfGFP_unformatted.gb" --o "in/sfGFP_formatted.gb"

Default settings:

- Remove common GG Sites: BsaI, BsmBI, and BbsI with `@AvoidPattern` tags
- Protect annotated regulatory sites with `@AvoidChanges` tags
- Enforce synonymous changes to all annotated CDSs using `@EnforceTranslation` tag
- Assemble oligos ≤ 250 bp with BsmBI using AATG/GCTT overhangs

Review `./iggypop format` output in your viewer, then generate oligos:

.. code-block:: bash

    ./iggypop.py gb --i "test/sfGFP_formatted.gb" --o "sfGFP"

GoldenBraid / MoClo Compatible CDSs
-----------------------------------

Default `./iggypop cds` sequences are GoldenBraid/MoClo compatible with 5'-BsaI-AATG and GCTT-BsaI-3'. Adjust `base_5p_end` and `base_3p_end` as needed.

.. image:: png/goldenbraid.png

Two-Step Assembly
-----------------

For sequences >3 kb (~18 fragments with 250 bp oligos), use the two-step assembly mode.

.. image:: png/two_step.png

Use the provided two_step yaml files:

.. code-block:: bash

    ./iggypop.py cds --i "in/RUBY.fasta" --o "RUBY_two_step" --yml "yaml/two_step_cds.yml"

Changing Cloning Overhangs & Assembly Enzyme
--------------------------------------------

You can change the external overhangs and enzyme for cloning:

.. code-block:: bash

    ./iggypop.py cds --i "test/RUBY.fasta" \
        --pcr_5p_cut GGTCTCA  --pcr_3p_cut AGAGACC \ # BsaI
        --base_5p_end AAAA    --base_3p_end GCCG \ # new cloning ends
        --ext_overhangs AAAA GCCG

Combining Runs
--------------

Use `--primer_index` to specify the starting row of the indexset file for new runs.

.. code-block:: bash

    ./iggypop.py cds --i "test/edibles.fasta" --o "edibles"
    ./iggypop.py cds --i "test/juiceables.fasta" --o "juiceables" --primer_index 11

Combine files into one fasta file for ordering:

.. code-block:: bash

    cat out/juiceables/juiceables_oligo_pool.fasta \
        out/edibles/edibles_oligo_pool.fasta > oligo_order.fasta

Use `assemble_fragments.py` to simulate oligo assembly:

.. code-block:: bash

    python scripts/assemble_fragments.py --i "oligo_order.fasta" --o "assembled_ej_oligos.fasta"

Versioning
----------

Use the `repeat` option for multiple optimized versions:

.. code-block:: bash

    ./iggypop.py cds --i "test/RUBY.fasta" --o "five_RUBYs" --codon_opt "match_codon_usage" --repeats 5

Chisel Only
-----------

`--mode no_hinge` outputs only dnachisel'd sequences.

Reports
-------

`--reports` enables dnachisel's report function, adding a sub-folder with changes for each sequence.

Quiet Mode
----------

`--quiet on` suppresses most terminal output.

Reproducibility
---------------

Set `--seed 123` to force a specific seed.

pPOP-vectors
------------

The pPOP vectors support one-step and two-step cloning and plant transformation. Find sequences [here].

iggyseq
-------

`iggyseq` identifies error-free clones via nanopore sequencing of barcoded colony PCR amplicons. See more details in the documentation.

hingesets
---------

`iggypop` uses `goldenhinges` to identify overhang solutions using precomputed `hingesets`.

.. image:: png/fidelity_plot.png

custom hingesets
----------------

Use `iggypop gagga` to create new `hingesets`:

.. code-block:: bash

    iggypop gagga                        \
        --set_size=20 --pop_size=1000    \
        --min_improve=.0005 --alpha 2.4  \
        --beta 2.4 --tournament_size 4

Process multiple runs with `process_gagga_runs.R`:

.. code-block:: bash

    Rscript scripts/process_gagga_runs.R --top_percent=2 --n_cliques=30

indexsets
---------

`indexsets` primers are designed to minimize cross-hybridization and unwanted amplifications.

custom indexsets
----------------

Use the pipeline below for custom indexsets:

.. code-block:: bash

    ./iggypop primers                   \
        --num_sequences 10 --opt_tm 60  \
        --opt_size 18 --gc_content 0.5  \
        --max_size 18 --min_size 18

Example `MFEprimer3` output for scoring primers:

.. image:: png/MFEprimer3_output.png
