
# iggypop


![Overview](png/overview.png)


**iggypop** is a pipeline for designing and synthesizing genes at ~$3.00 per kb in oligo costs. It uses the Edinburgh Genome Foundry's [`dnachisel`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) to optimize sequences and [`goldenhinges`](https://github.com/Edinburgh-Genome-Foundry/GoldenHinges) to fragment them into indexed oligos that can be reassembled by Golden Gate cloning. Once assembled, sequence-verified constructs are identified by nanopore sequencing of barcoded amplicons.

*iggypop* enables economical end-to-end construction of 100s of genetic parts in days.

## Installation

### Linux

```bash
# requires python 3.8
git clone github.com/cutlersr/iggypop
cd iggypop
python3.8 -m venv .venv
source .venv/bin/activate
chmod +x setup.sh
./setup.sh
```

### Docker

```bash
git clone github.com/cutlersr/iggypop
cd iggypop
docker build -t iggypop .
docker run -it -v $(pwd):/app iggypop
```



## Working with Fasta Formatted Coding Sequences

Coding sequences are domesticated, fragmented, indexed, and appended with cut sites to yield oligonucleotides that can be amplified with gene-specific primers and then assembled. Sequence domestication and optimization is conducted using `dnachisel`; optimization parameters can be set in a YAML file using dnachisel's built-in [`specifications`](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html); several YAML files used in our common workflows are in the  [`yaml`](../yanml/) folder.

```bash
./iggypop.py cds --i "test/10_TFs.fasta" --o "10_TFs"
```


The default parameters can be modified by creating new YAML files or from the command line. To codon optimize with an *E. coli* codon table, use BsaI sites for assembly, and synthesize 300 bp oligos:

```bash
./iggypop.py cds                                        \
    --i "test/10_TFs.fasta" --o "10_TFs_coli_mcu"       \
    --base_5p_end "GGTCTCA" --base_3p_end "AGAGACC"     \ # BsaI instead of BsmBI
    --codon_opt "match_codon_usage" --species "e_coli"  \
    --oligo_length 300  # default is 250
```


The default CDS settings design ORFs that:

- Lack common Golden Gate cloning sites (BsaI, BsmBI, BbsI, SapI, BtgZI)
- Enforce synonymous changes
- Assemble from oligos ≤ 250 bp with BsmBI
- Lack hairpins or repeats >12 bp
- Are GoldenBraid / MoClo compatible (inner 5'-BsaI-AATG...GCTT-BsaI-3')


## Working with GenBank Formatted Sequences

The parameters for optimizing GenBank files differ and use annotations added to your GenBank file using dnachisel's [`genbank API`](https://edinburgh-genome-foundry.github.io/DnaChisel/genbank/genbank_api.html). `iggypop.py format` allows easy parameter setting in a YAML file:

```bash
# Format a genbank file using the default domesticate_gb.yml file
./iggypop.py format --i "test/sfGFP_unformatted.gb" --o "test/sfGFP_formatted.gb"
```

Default settings:

- Remove common GG Sites: BsaI, BsmBI, BbsI, SapI, and BtgZI with `@AvoidPattern` tags
- Protect annotated regulatory sites with `@AvoidChanges` tags
- Enforce synonymous changes to all annotated CDSs using `@EnforceTranslation` tag
- Assemble oligos ≤ 250 bp for BsmBI assembly using AATG/GCTT overhangs

Check the output in your favorite viewer, then generate your oligos:

```bash
./iggypop.py gb --i "test/sfGFP_formatted.gb" --o "sfGFP"
```



## Level 0 / Level-α CDS Modules

The default `iggypop` settings create GoldenBraid/MoClo-compatible level 0 coding sequences with 5'-BsaI-AATG and GCTT-BsaI-3'. ends Adjust the `base_5p_end` and `base_3p_end` parameters to modify this behavior. 

![GoldenBraid](png/goldenbraid.png)


## pPOP-vectors

The pPOP [`vectors`](../vectors/) support one-step and two-step cloning of level 0 / level α parts; the pPlantPOP-BsmBI vector supports `iggypop` assemblies of MoClo-compatible parts for direct testing *in planta* via Agrobacterium-mediated transformation.


## Two-Step Assemblies

For longer sequences >3 kb (~18 fragments with 250 bp oligos), we recommend that you use the two-step assembly mode, which breaks sequences into  "step one" blocks which are assembled from oligo pools using BbsI and the pPOP-BbsI vector. Once sequence validated step one clones are identified, the final genes are assembled using Golden Gate cloning using BsmBI into pPOP-BsmBI or pPlantPOP-BsmBI. 

![Two-Step Assembly](png/two_step.png)

To do this, use the provided `two_step` YAML files:

```bash
./iggypop.py cds --i "test/RUBY.fasta" --o "RUBY_two_step"   \
                --yml "yaml/domesticate_two_step_cds.yml"
```



## Changing Cloning Overhangs & Assembly Enzyme

You can change the external overhangs and enzyme for cloning:

```bash
./iggypop.py cds --i "test/RUBY.fasta"          \
    --pcr_5p_cut GGTCTCA  --pcr_3p_cut AGAGACC  \ # BsaI
    --base_5p_end AAAA    --base_3p_end GCCG    \ # new cloning ends
    --ext_overhangs AAAA GCCG
```


## Combining Runs

Use "--primer_index" to specify the starting row of the indexset file for new runs.

```bash
./iggypop.py cds --i "test/edibles.fasta" --o "edibles"
./iggypop.py cds --i "test/juiceables.fasta" --o "juiceables" --primer_index 11
```

Then combine files into one fasta file to create an oligo pool file for ordering:

```bash
cat out/juiceables/juiceables_oligo_pool.fasta \
    out/edibles/edibles_oligo_pool.fasta > oligo_order.fasta
```

You can use `assemble_fragments.py` to simulate golden gate assembly and confirm that none of your index primers are used on more than one gene:

```bash
python scripts/assemble_fragments.py --i "oligo_order.fasta"          \
                                     --o "assembled_ej_oligos.fasta"
```


## Hinging

Although our default settings work well, you can see marginal increases in the predicted fidelity of your assemblies by increasing the search `radius` around target cut sites and/or increasing the number of hinge solutions evaluated (`n_tries`). We only recommend changing these if your targets generate low fidelities (<95%) with defaults. Increasing the `radius` parameter can sometimes generate more oligos for the same target, so change this parameter cautiously.

```bash
./iggypop.py cds --i "test/RUBY.fasta" --o "five_RUBYs"  \
				 --radius 10                             \ #default is 8
				 --n_tries 50                            \ #default is 10
```

## Chiseling Only

`--mode no_hinge` will outputs only dnachisel'd sequences.


## Hinging Only

`--mode no_mods` will run the hinging process and output indexed oligo for input sequences without making any changes to your input sequences.


## Versioning

Use the `repeat` option if you'd like to generate multiple optimized versions of a sequence:

```bash
./iggypop.py cds --i "test/RUBY.fasta" --o "five_RUBYs"        \
                 --codon_opt "match_codon_usage" --repeats 5   
```



## Reports

`--reports` enables dnachisel's report function, adding a sub-folder with changes for each sequence.


## Quiet Mode

`--quiet on` suppresses most terminal output.


## Reproducibility

Set `--seed 123` to force a specific seed. The log files list the seeds used for each sequence.



## iggyseq

`iggyseq` identifies error-free clones via nanopore sequencing of barcoded colony PCR amplicons. See ... more details in the documentation.


## hingesets

iggypop uses *goldenhinges* to identify overhang solutions using precomputed `hingesets`. The data below show the typical fidelities using a set of 4500 plant transcription factors.

![Fidelity Plot](png/fidelity_plot.png)

### Custom hingesets

Use `gagga` to create new hingesets:

```bash
./iggypop.py gagga --set_size=20 --pop_size=1000      \
                   --min_improve=.0005 --alpha 2.4    \
                   --beta 2.4 --tournament_size 4
```

Process multiple runs with `process_gagga_runs.R`:

```bash
Rscript scripts/process_gagga_runs.R --top_percent=2 --n_cliques=30
```


## indexsets

Our primers used for amplifying fragments from pools were designed to minimize cross-hybridization and unwanted amplifications.

### Custom indexsets

Use the pipeline below for custom indexsets:

```bash
./iggypop.py primers                   \
    --num_sequences 10 --opt_tm 60     \
    --opt_size 18 --gc_content 0.5     \
    --max_size 18 --min_size 18
```
