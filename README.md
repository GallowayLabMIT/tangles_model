# Supercoiling-mediated feedback rapidly couples and tunes transcription
[<img src="https://img.shields.io/badge/code_license-MIT-green">](./LICENSE)
[<img src="https://img.shields.io/badge/text_license-CC--BY--4.0-green">](https://creativecommons.org/licenses/by/4.0/)

| Item          | DOI           |
| ------------- |:-------------:|
| Article (preprint)       | [![Article DOI](https://img.shields.io/badge/Article_DOI-10.1101/2022.04.20.488937-green)](https://doi.org/10.1101/2022.04.20.488937)     |
| Software (this repository)      | [![Software DOI](https://img.shields.io/badge/Software_DOI-10.5281/zenodo.7054395-blue)](https://doi.org/10.5281/zenodo.7054395)                                                              |
| Dataset for article figures     | [![Dataset DOI](https://img.shields.io/badge/Dataset_DOI-10.5281/zenodo.7041641-blue)](https://doi.org/10.5281/zenodo.7041641)     |

Included in this repo is complete code to:
1. Simulate gene expression under the control of supercoiling for any user-defined circuits.
2. Regenerate simulation data for the preprint.
3. Regenerate all figures from simulation data for the preprint.

If you'd like to simulate your own gene circuits, continue reading! If you just want to replicate the paper results, jump to [Replicating paper results](#replicating-paper-results). Instead of regenerating the data using the Julia scripts, you may want to [download](https://nextcloud.meson.us/index.php/s/nPeowqMAKLrETgB) the already-generated datasets used in generating all figures including in this paper.

## Simulating your own gene networks
The full model described in our work is mostly contained within [TanglesModel.jl](src/TanglesModel.jl), a Julia implementation.

### Install and code check
After installing Julia, try to install the needed dependencies. To do this, from the root folder of this repository, you can run:
```
julia --project=.
```
You will see a prompt like `julia>`. Enter package management mode by pressing `]` and type the following two commands (here, `pkg>` is the Julia prompt; do not type it!)
```
julia> ]
(TanglesModel) pkg> resolve
(TanglesModel) pkg> instantiate
(TanglesModel) pkg> precompile
```
You may see Julia output after typing these; resolve checks for any package updates, instantiate downloads all of the dependencies, and precompile precompiles all packages, including TanglesModel.
Precompiling also does some basic checks, so the model code is ready to go if you do not see an error! Successful output looks like:
```
Precompiling project...
  1 dependency successfully precompiled in 94 seconds (254 already precompiled)
```

For the Python part of the analysis, it is most convenient if you install the requirements into a virtual environment. To do this from the root, do:
```
python -m venv env       # on most OS's
python3 -m venv env      # on some modern MacOS's
.\env\scripts\activate   # on Windows
source env/bin/activate  # on Linux, Mac
pip install -r requirements.txt
```

### Running a simulation
To run a simulation, you need to decide on your system's _boundary conditions_ and the _type, location, and direction_ of genes. If you are looking for examples, you can look at the following:
- `src/simulate_fig_base_model_examples.jl` for an example of two uncoupled genes in various orientations.
- `src/simulate_fig_toggles.jl` for an example of two genes with relatively simple coupling.
- `src/simulate_fig_zinani.jl` for an example of two genes with complex, multi-species reactions occurring.

#### Specifying simulation parameters
All simulations need a `SimulationParameters` object to encode the various physical parameters used. To replicate the system used in a majority of the paper, use:

```
sim_params = SimulationParameters(
        DEFAULT_SIM_PARAMS.mRNA_params,
        DEFAULT_SIM_PARAMS.RNAP_params,
        DEFAULT_SIM_PARAMS.DNA_params,
        DEFAULT_SIM_PARAMS.temperature,
        DEFAULT_SIM_PARAMS.topoisomerase_rate,
        DEFAULT_SIM_PARAMS.mRNA_degradation_rate,
        true,    # Supercoiling dependence
        0.025    # alpha
    )
```
However, all other parameters could be changed as well

#### Boundary conditions
Two types of boundary conditions can be created. Circular boundary conditions are created via:
```
CircularBoundaryConditions(length_in_nm)
```
and represent a circular piece of DNA. Importantly, genes _cannot_ be defined spanning the origin, so select an origin that is not in a gene body.

Linear boundary conditions can be created via:
```
LinearBoundaryConditions(length_in_nm, left_is_free, right_is_free)
```
Here, the second and third parameters specify whether the left and right boundaries are fixed or free. Fixed boudary conditions act as walls that supercoiling can accumualte against, whereas free boundary conditions relax supercoiling away.


An example of a 10,000bp linear boundary conditions with a fixed left-end with a free right-end could be:
```
bcs = LinearBoundaryConditions(10000.0 * 0.34, false, true)
```

#### Specifying genes
Two types of genes can be specified. If your gene has some static expression rate, you can directly specify it as an UncoupledGene, assigning it a transcript index (1-indexed!):
```
UncoupledGene(basal_rate_in_inverse_seconds, transcript_idx, gene_start_in_nm, gene_end_in_nm)
```
Such a gene will still be supercoiling-dependent. Importantly, the gene orientation is inferred based on the start and end; e.g. a gene will be "pointing" to the left if the gene start is greater than the gene end.

An example of two uncoupled genes, convergently oriented, could be:
```
genes = [UncoupledGene(1/120, 1, 3000 * 0.34, 4000 * 0.34), UncoupledGene(1/120, 2, 7000 * 0.34, 6000 * 0.34)]
```

If you need to introduce additional coupling, then `CoupledGene`s are your ticket:
```
CoupledGene(basal_rate_in_inverse_seconds, transcript_idx, gene_start_in_nm, gene_end_in_nm, coupling_lambda)
```
What is `coupling_lambda`? This is a function that should take two arguments and return a number that is multiplied into the basal rate.
The first argumentis the current concentrations of the discrete simulation components, including transcripts.
The second argument is the current simulation time.

An example of a coupled gene, that only turns on after a 10,000 second delay is:
```
CoupledGene(1/120, 1, 3000 * 0.34, 4000 * 0.34, (_,t) -> t > 10000)
```

An example of a coupled gene that depends on the number of the first gene's transcripts is:
```
K = 10
CoupledGene(1/120, 2, 3000 * 0.34, 4000 * 0.34, (discrete_vals,_) -> K / (K + discrete_vals[1]))
```

Both `UncoupledGene`s and `CoupledGene`s create a single mRNA transcript when polymerases finish transcribing. If you'd like to create
multiple transcripts from the same "gene", you can use `MultiUncoupledGene` and `MultiCoupledGene`, which allow you to pass a list
of transcript indicies.
```
MultiUncoupledGene(basal_rate_in_inverse_seconds, transcript_idx_list, gene_start_in_nm, gene_end_in_nm)
MultiCoupledGene(basal_rate_in_inverse_seconds, transcript_idx_list, gene_start_in_nm, gene_end_in_nm, coupling_lambda)
```

#### Specifying a discrete configuration
With a list of genes, you can specify a DiscreteConfiguration. If you have no other stochastic reactions, then this is as simple as:
```
dconfig = DiscreteConfig([UncoupledGene(...), UncoupledGene(...)]) # or a list of other types of Genes
```

If you have other discrete reactions, then you need to pass:
```
DiscreteConfig(gene_list, n_other_discrete_species, discrete_reaction_list)
```
`n_other_discrete_species` is the number of other species that you want to simulate; they get (1-indexed) indicies after the transcripts.
For example, if your genes create three different transcripts, then the first other discrete component will have index 4.

Reactions are specified as pairs of the form:
```
(propensity_function(discrete, t) => [(species_index => change_in_species_number), (species_index => change_in_species_number), ...])
```

For an example, let's consider an irreversible dimerization protein reaction, that only "turns on" after one hour has elapsed. If
the reaction takes the monomer (index 2) and turns it into a dimer (index 3) with the reaction rate of `2 A -> B` occuring with rate
`r = 2 [A]^2`, then the reaction would be:
```
reaction = (discrete, t) -> 2 * discrete[2]^2 * (t > 3600) => [2 => -2, 3 => 1]
```

### Selecting the type of output
There are many types of output that can be saved. The most relevant ones are:

```
simulate_full_examples(filename, n_simulations, comment, sim_params, bcs, dconfig, end_time, n_timesteps)
```
which saves the complete simulation state at evenly-spaced timesteps. The complete simulation state includes the location of every polymerase
along with the excess twist at those locations, in addition to the number of transcripts.

```
simulate_discrete_runs(filename, n_simulations, comment, sim_params, bcs, dconfig, end_time, n_timesteps, discrete_ics, extra_metadata)
simulate_discrete_runs(filename, n_simulations, comment, sim_params, bcs, dconfig, end_time, n_timesteps, extra_metadata)
```
which just saves the number of transcripts (and other discrete species) as a function of simulation time.

```
simulate_sc_rnap_dynamics(filename, n_simulations, comment, sim_params, bcs, dconfig, n_supercoiling_steps, end_time, n_timesteps, extra_metadata)
```
which saves the supercoiling density across the simulated system, in addition to saving the polymerase initiation and termination times. This data
is used for supercoiling density and burst dynamics calculations.

The relevant arguments are:
- `filename`: the filename of the `.h5` file to write simulation results into. Results are _appended_, so you can redirect multiple runs into the same file.
- `n_simulations`: the number of simulations to perform.
- `comment`: a string that can help in identifying simulations written to the `.h5` file.
- `sim_params`: a SimulationParameters object defining important parameters (see above).
- `bcs`: a BoundaryConditions object encoding the boundaries of the simulation (see above).
- `dconfig`: A DiscreteConfig object encoding genes and other discrete reactions to simulate (see above).
- `end_time`: the length of the simulation, which will run from zero seconds to `end_time` seconds.
- `n_timesteps`: the number of timepoints, evenly spaced, to save to the output file.
- `discrete_ics`: if specified, sets the initial number of each of the discrete species (gene transcripts and any user-defined additional species). Otherwise all zeros.
- `extra_metadata`: a list of `(key_string => value_string)` pairs that gets written as additional `h5` metadata.



## Replicating paper results
There is a three-step pipeline for simulating results and generating the figures in the paper.

1. Simulate the raw simulation data, potentially spread across many workers in parallel.
2. Merge datasets together from multiple workers.
3. Pre-process and generate figures!

If you don't want to regenerate the data yourself, you can download the the raw, Julia-combined dataset (`unprocessed_datasets.zip`, 15.1GB)  and post-processed datasets (`preprocessed_datasets.zip`, 2.8GB) from [Zenodo](https://doi.org/10.5281/zenodo.7041641). You can unzip both the unprocessed and post-processed datasets into your data directory folder; they should not be in subdirectories.

After downloading the dataset or regenerating it, you can move to [setting a data directory](#setting-a-data-directory) and [replicating figures](#replicating-figures). In any case,
you will need to perform the Python dependency installation steps in the [code check](#install-and-code-check) section. If you are also
redoing the raw simulation data, you need to perform the Julia-specific steps.

### Simulating raw simulation data

Generating the raw simulation data is a time-intensive process, but can be
run in parallel (e.g. on a powerful computer or on a computing cluster). The main simulation scripts
are:

- `src/simulate_fig_base_model.jl, src/simulate_fig_base_model_spacing.jl` for Figure 3
- `src/simulate_base_model_examples.jl` for Figure 4
- `src/simulate_sc_bursting.jl` for all supercoiling density/bursting figures
- `src/simulate_fig_toggles.jl` for Figure 5
- `src/simulate_fig_zinani.jl` for Figure 6
- `src/simulate_fig_hyperparams.jl` for Figure S2.

All of these scripts are meant to be run in parallel, and take two command-line arguments. The first
is the number of parallel runs, and the second is the (zero-indexed) worker index. For example, to split
the base-model simulations across ten workers, you would run the following from the root directory:
```
julia --project=. src/simulate_fig_base_model.jl 10 0
julia --project=. src/simulate_fig_base_model.jl 10 1
...
julia --project=. src/simulate_fig_base_model.jl 10 9
```

If you are running this single-threaded, you would use `1 0` as the options. To reduce time requirements, you can
decrease the `n_repeats`, `n_examples_per_node`, or `n_full_examples_per_node` parameters, at the cost of
simulating a smaller ensemble.

If you happen to be using the [MIT Supercloud](https://supercloud.mit.edu), there are `slurm` launch scripts that
automate this process on the cluster.

### Sidebar: simulating the Zinani data
[Zinani et al.](https://www.nature.com/articles/s41586-020-03055-0) provides MATLAB code, that we use as a comparison in our work.
However, [the authors' code](https://github.com/ozbudak/zinani_genepairing) is on Github but is unlicensed, and thus is All Rights Reserved,
so we cannot legally include their code in our repository.

We made some small modifications to speed up execution of their code without affecting logic. We have provided a diff of our changes in
`src/zinani_paired.m.diff` and `src/zinani_unpaired.m.diff`.
By peforming the changes in those diff files to [genepaired_scenario1.m](https://github.com/ozbudak/zinani_genepairing/blob/main/genepaired_scenario1.m)
and [geneunpaired_scenario1.m](https://github.com/ozbudak/zinani_genepairing/blob/main/geneunpaired_scenario1.m), you can recreate the versions
of these MATLAB scripts we used here.

### Merging datasets
There are Python scripts that are responsible for taking worker-generated `.h5` files and turning them into merged datasets.
You can run the following scripts (some may have to be edited to update file paths, depending on where you generated the worker data) to perform this merge:
```
src/merge_model_bc_selection.py
src/merge_sc_bursting_sims.py
src/merge_spacing_sweep.py
src/merge_summaries.py
src/merge_zinani_matlab.py
src/merge_zinani_sims.py
src/merge_hyperparameters.py
```

### Setting a data directory
We use [rushd](https://github.com/GallowayLabMIT/rushd) to keep our data analysis tidy. Part of this is specifying
a *data directory*, potentially on another drive or outside of this repository, as the raw simulation data is several gigabytes.
After you decide where to put the data locally, create a file called `datadir.txt` that contains only the absolute path to the folder
in which you placed the data. For example, if you placed the data at `D:/tangles_data`, e.g. in this folder you have the unzipped `.h5` (unprocessed) or `.gzip` (processed) files , then your `datadir.txt` file should contain one line:

```
D:/tangles_data
```

## Replicating figures
With the datadir in place, you should be able to replicate our figures. Run the Jupyter notebook at `notebooks/modeling_paper_figures.ipynb` to replicate the figure results.
If you want to replicate the supplemental movies, run `notebooks/short_movies.ipynb`. The annotated movie (Supplemental Movie 1) was manually edited in Adobe Premiere, and is not fully reproducible in code.

There is a manual step in order to prepare figures that contain illustrations. For most of the main-text figures, running the notebook will create an SVG file containing
the graphics, under a filename like `fig_base_model_mpl.svg` (e.g. the base model figure, Figure 3, generated using `matplotlib`). This can be combined with the PDF-exported
version of the Adobe Illustrator files. to create the final version.

Don't have access to Adobe products? Just rename the `.ai`'s in `writeups/figures/modeling_paper` to `.pdf`'s and it should be importable.

## License and reuse
All figures and text (e.g. content within the `writeups` subfolder) are licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
All figure generating code and simulation code are licensed under the [MIT license](./LICENSE).

## Why is it called "TANGLES model"?
The original name for this project was Topologically-Affected Networks of Genes Linking Expression to Syntax, thus TANGLES :)
