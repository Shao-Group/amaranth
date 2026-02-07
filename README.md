# ![image](docs/image.jpeg) Amaranth: Single-Cell Transcript Assembler

<a href="http://bioconda.github.io/recipes/amaranth-assembler/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="install with bioconda"></a> <a href="https://github.com/Shao-Group/amaranth/releases"><img src="https://img.shields.io/github/v/tag/Shao-Group/Amaranth" alt="GitHub tag (latest SemVer)"></a> <a href="https://github.com/Shao-Group/amaranth/blob/master/LICENSE"><img src="https://img.shields.io/github/license/Shao-Group/Amaranth" alt="GitHub License"></a>

Amaranth is a reference-based transcript assembler specifically optimized for single-cell RNA-seq data. The development of Amaranth has been based on the [Scallop2](https://github.com/Shao-Group/scallop2)/[Scallop](https://github.com/Kingsford-Group/scallop) assembler series.

# Installation

Amaranth is available in bioconda channel. It is most recommended to install amaranth by using [Pixi](#installation-using-pixi) or [Mamba](#installation-using-mamba-micromamba). It can also be installed by Conda, but conda may be much slower in solving the environment, especially when user’s conda version is not up-to-date or when trying to install in a populated environment. Amaranth supports Linux and macOS. Windows users please consider using [wsl](https://learn.microsoft.com/en-us/windows/wsl/install) or other Unix-like system.

## Installation using Pixi

[Pixi](https://pixi.prefix.dev/latest/) is a new fast, modern, and reproducible package management tool and works cross-programming languages, including Python, C++ and Rust. It can install software from conda channels. It is the most recommended package management tool to install Amaranth. 

If not yet, users need to first install Pixi:

```bash
# install pixi
curl -fsSL https://pixi.sh/install.sh | sh
```

You may need to restart your terminal or run `source ~/.bashrc` (Linux) / `source ~/.zshrc` (macOS) to make the change effective. 

> If the installation is successful, you can type `pixi` in the terminal and pixi's help message should be printed.

Pixi organizes software and dependencies in projects (similar to virtual environments). To start, you need to initiate a pixi project in a directory using `pixi init <project_dir>`. This will create a `pixi.toml` config file and a hidden `.pixi` directory in the `<project_dir>`.  You will also need to specify `conda-forge` and `bioconda` channels. Pixi will search for software from those channels.

```bash
# Initiate pixi project in a directory (anywhere, e.g., your current dir)
pixi init . --channel conda-forge --channel bioconda
```

Once inside your `<project_dir>`, use the following command to install amaranth:

```bash
# to install
pixi add bioconda::amaranth-assembler
```

> Successful installation will print message like “Added bioconda::amaranth-assembler”.
>
> After installation, `amaranth` is installed in `.pixi` subdirectory, managed by Pixi. You won't see a new file in your current directory. To use, run `pixi run amaranth <arguments>`.

After installation, users can use `pixi run <tool_name> <tool_arguments>` to use the software. For example, the following command will print help messages.

```bash
# To test installation, you can try print the help message
pixi run amaranth --help
```

You can also download example data from our GitHub repo ([example data](https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam)) and test amaranth. After downloading, place the example file in your project directory:

```bash
# download example data
wget https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam
# to assemble example-input.bam and produce test_output.gtf file
pixi run amaranth -i example-input.bam -o test_output.gtf
```

Note: `pixi run` commands must be executed from within the project directory.

To directly use tools without `pixi run`, users can invoke `pixi shell`, similar to `conda activate env_name`. Then users can use `amaranth` anywhere, even outside the `<project_dir>`.

```bash
# activate pixi environment
pixi shell
# users can directly run tools without `pixi run`
amaranth -i example-input.bam -o test_output.gtf
```

> `pixi run <tool>` only works inside its corresponding pixi project_dir. After invoking `pixi shell` (from project_dir), users can use the software anywhere, including outside the project_dir.
>
> To leave pixi shell, run `exit`.  

## Installation using Mamba (Micromamba)

[Mamba](https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++. It is fully compatible with conda, and much faster. Micromamba is a tiny version of the mamba package manager. In this installation guide, we use `micromamba` as an example. If you already have conda or mamba installed, you can skip the installation step and replace `micromamba` with `conda` or `mamba` in the commands below. However, conda could be noticeably slower than mamba and micromamba, especially on older versions. If using conda, we recommend at least Version 24 or newer.

If not yet, users need to first install Micromamba:

```bash
# install Micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

You may need to restart your terminal or run `source ~/.bashrc` (Linux) / `source ~/.zshrc` (macOS) to make the change effective.

> If the installation is successful, you can type `micromamba` in the terminal and micromamba’s help message should be printed.

(Optional) It is recommended to create a new environment, to minimize conflicts and speed-up installation. You can also activate any existing environment that you want to install in.

```bash
# Optionally, create and activate a new environment
micromamba create -n amaranth_env -c conda-forge -c bioconda
micromamba activate amaranth_env
```

The following `micromamba` command will install amaranth:

```bash
# to install amaranth
micromamba install -c conda-forge bioconda::amaranth-assembler
```

> Successful installation will print message like "Executing transaction: done" or "Transaction finished". If Micromamba takes a long time to resolve the environment (e.g., 15+ minutes), consider using a fresh new environment.
>
> After installation, `amaranth` is available as a command in your terminal (i.e. from `$PATH`), but you won't see a new file in your current directory. The binary is installed inside the micromamba environment folder. Use `which amaranth` to see the exact location.

After installation, users can directly call `amaranth` to use the software. For example, to print help message or test on a small example data (download [data](https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam) from our GitHub repo):

```bash
# To test installation, you can try print the help message
amaranth --help
# or download and test on an example dataset
wget https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam
amaranth -i example-input.bam -o test_output.gtf
```

> If you created a new environment, remember to run `micromamba activate amaranth_env` each time you open a new terminal before using amaranth.
>
> To deactivate the environment, run `micromamba deactivate`.

## Use amaranth in Docker container

Alternatively, amaranth is also available in docker. Please ensure `Docker daemon` is running in the background before using docker. Root privilege may be required to run Docker daemon.

To pull (download) the image, users need to replace the `<tag>` with appropriate value from biocontainer ([biocontainers / amaranth-assembler](https://quay.io/repository/biocontainers/amaranth-assembler?tab=tags&tag=latest)). 

```bash
# replace the <tag> with valid value
docker pull quay.io/biocontainers/amaranth-assembler:<tag>
```

For example, the docker container tag for `v0.1.0`:

```bash
# actual tag for v0.1.0
docker pull quay.io/biocontainers/amaranth-assembler:0.1.0--h5ca1c30_0
```

To explore the container interactively and see what's inside:

```bash
docker run -it quay.io/biocontainers/amaranth-assembler:0.1.0--h5ca1c30_0 /bin/bash
# you can try to explore amaranth in the container, for example, print help message
amaranth --help
```

To run amaranth with your local files mounted (non-interactive):

```bash
docker run -v $(pwd):/data quay.io/biocontainers/amaranth-assembler:0.1.0--h5ca1c30_0 \
	amaranth -i /data/<your_input.bam> -o /data/<output_prefix>
```

The `-v $(pwd):/data` mounts your current directory to `/data` inside the container so you can access your files. That’s why `/data/` is also added before both input and output files’ local paths. Otherwise the input/ output files won’t be accessible to Docker.

## Installation from source code (with external libraries)

Download the source code that contains external libraries (htslib, boost, and zlib) from [releases](https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/amaranth-0.1.0-full.tar.gz). Use the following commands to uncompress and install:

```
tar xzvf amaranth-0.1.0-full.tar.gz
cd amaranth-0.1.0
./build.sh
```
The executable will appear as `src/amaranth`.


## Other installation methods and FAQs

For other installation methods, such as installing from source code (without external libraries), please read [INSTALL.md](./INSTALL.md). If you have questions on installation, please check [Installation FAQ](./INSTALL.md/#Installation-FAQ) or open an issue on github. 

# Usage

Assuming `amaranth` is available from command line (which may vary depending on the installation method), the general usage is:

```
amaranth -i <input.bam> -o <output>
```

> Users may need to replace `amaranth` with appropriate command with respect to their installation methods:
>
> - Pixi: go to corresponding `<project_dir>` and (1) run `pixi run amaranth -i <input.bam> -o <output>` or (2) run `pixi shell` and then run `amaranth -i <input.bam> -o <output>`.
> - Mamba/conda: please make sure the correct environment is activated before calling `amaranth -i <input.bam> -o <output>`.
> - Docker: `docker run -v $(pwd):/data quay.io/biocontainers/amaranth-assembler:<tag> amaranth -i /data/<input.bam> -o /data/<output>`. Remember to replace `<tag>` with the actual tag value.
> - Compiled from source, please use the path to amaranth executable binary. It is usually  `./src/amaranth`.

The `input.bam` is a read alignment file generated by an RNA-seq aligner. 

To correctly assemble single-cell RNA-seq reads (e.g. Smart-seq3), the bam file should have SAM optional field tag `BC` which stores cell barcodes and tag `UB` which stores UMI barcodes. See bam data in `example` directory as an example.

If you want to assemble each single cell independently, run one command for each cell. If you want to perform meta-assembly (which leverages information across all cells to improve individual cell assemblies), see [Meta-assembly Usage](#Meta-assembly-Usage).

Make sure that the bam file is sorted; otherwise run `samtools` to sort it:

```
samtools sort input.bam > input.sort.bam
```

The reconstructed transcripts shall be written as gtf format into `output.gtf`.

## Meta-assembly Usage

The usage of `amaranth` for meta-assembly is:

```bash
amaranth --meta -i <merged.bam> -o <output>
```

The argument `--meta` must be supplied to perform meta-assembly. 

The `merged.bam` is the **sorted** and **merged** read alignment file of all cells. The bam file should have SAM optional field tag `BC` which stores cell barcodes and tag `UB` which stores UMI barcodes.

If user has separate sorted read alignment files of each single cell, use `samtools` to merge them.

```bash
# merge n bam files with 32 threads
samtools merge -@32 -o merged.bam 1.bam 2.bam ... n.bam
```

The reconstructed transcripts for each cell will be written as gtf format into `output.<BC>.gtf`, where `<BC>` is the cell barcode (sam tag `BC`) of each cell. The meta-assembly (union of transcripts from all cells) will be written into `output.meta.gtf`.

To achieve the best performance, it is recommended to union a cell's transcripts from both amaranth's meta-assembly run and individual  (as a single cell) assembly run. Union of transcripts can be done using tools such as [TACO](https://tacorna.github.io/) or [gtfmerge](https://github.com/Shao-Group/rnaseqtools?tab=readme-ov-file#gtfmerge).

## Parameters

Here is a list of supported parameters. Please refer to additional explanations below the table.

 Parameters | Default Value | Description
 ------------------------- | ------------- | ----------
 --help  | | print usage of amaranth and exit
 --version | | print version of amaranth and exit
 --meta | not used | to perform meta-assembly. 
 --preview | | show the inferred `library_type` and exit
 --verbose | 1 | chosen from {0, 1, 2}
 -f/--transcript_fragments    | | file to which the assembled non-full-length transcripts will be written to
 --library_type               | empty | chosen from {empty, unstranded, first, second}; If empty, Amaranth will try to infer automatically. 
 --assemble_duplicates		  | 10 | the number of consensus runs of the decomposition
 --min_transcript_coverage    | 1.5 | the minimum coverage required to output a multi-exon transcript
 --min_single_exon_coverage   | 20 | the minimum coverage required to output a single-exon transcript
 --min_transcript_length_base      |150 | the minimum base length of a transcript
 --min_transcript_length_increase  | 50 | the minimum increased length of a transcript with each additional exon
 --min_mapping_quality        | 1 | ignore reads with mapping quality less than this value
 --max_num_cigar              | 1000 | ignore reads with CIGAR size larger than this value
 --min_bundle_gap             | 100 | the minimum distances required to start a new bundle
 --min_num_hits_in_bundle     | 5 | the minimum number of reads required in a bundle
 --min_flank_length           | 3 | the minimum match length required in each side for a spliced read
 --min_splice_boundary_hits    | 1 | the minimum number of spliced reads required to support a junction
 --min-umi-reads-bundle | 1 | (int) Bundle with less UMI reads than this threshold will be ignored
 --min-umi-ratio-bundle           | 0.0           | (float) Bundle with lower UMI reads ratio than this threshold will be ignored
 --both-umi-support               | not used      | If set a bundle need to satisfy both `min-umi-reads-bundle` and `min-umi-ratio-bundle`. Otherwise, satisfy either of them is ok
 --min-umi-reads-start-exon | 1 | (int) minimum number of UMI reads support of the first exon in a valid transcript
 --remove-retained-intron          | used          |
 --no-remove-retained-intron       | not used      |
 --max-ir-part-ratio-v | 0.5 | the ratio threshold of retained node to skip edge for partial introns (if greater than threshold, consider true transcript)
 --max-ir-part-ratio-e | 0.5 | the ratio threshold of retained node's edge to skip edge for partial introns (if greater than threshold, consider true transcript)
 --max-ir-full-ratio-v | 1.0 | the ratio threshold of retained node to skip edge for full introns (if greater than threshold, consider true transcript)
 --max-ir-full-ratio-e | 0.5 | the ratio threshold of retained node's edge to skip edge for full introns (if greater than threshold, consider true transcript)
 --max-ir-full-ratio-i | 10.0 | the ratio threshold of retained node to its own edge for full introns (if greater than threshold, consider true RETENTION)
 --max-ir-umi-support-full | 3 | (int) maximum number of UMI reads to support a partial exon rather than full intron retention
 --max-ir-umi-support-part | 5 | (int) maximum number of UMI reads to support a partial exon rather than partial intron retention
 --remove-pcr-duplicates          | 1             | 0 (not remove) or 1 (remove w.r.t alignment coordinates and CIGAR string)
 --min-cb-ratio | 0.3 | (float) for meta-assembly, minimum ratio of exons in a transcript supported by a cell's barcode 


1. For `--verbose`, 0: quiet; 1: one line for each splice graph; 2: details of graph decomposition.
2. `--min_transcript_coverage` is used to filter lowly expressed transcripts: amaranth will filter
  out transcripts whose (predicted) raw counts (number of moleculars) is less than this number.
3. `--min_transcript_length_base` and `--min_transcript_length_increase` is combined to filter
  short transcripts: the minimum length of a transcript is given by `--min_transcript_length_base`
    \+ `--min_transcript_length_increase` * num-of-exons-in-this-transcript. Transcripts that are less
    than this number will be filtered out.



# Example

Example test data is provided in `example/`. The example data is also available in [Release](https://github.com/Shao-Group/amaranth/releases/tag/v0.1.0).

Users can also use the following link to download from GitHub Release:

```bash
wget https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam
```

> ‼️ Note that `wget` works only with GitHub Release, but it does NOT work with GitHub repo files. It will be truncated if you use wget to download GitHub files directly from a repo,

You can use `sha256` to test the integrity of downloaded data.

```bash
sha256 example-input.bam  # should be `59e036720e6539336600409bb7a466dd82a51e8ad98a30016951aea42f21fba6`
```

Users can use the following command to do a basic example run (assuming you are in the same directory as `example-input.bam`)

```bash
amaranth -i example-input.bam -o test_output
```

The assembled transcripts will be in `test_output.gtf`.

# Questions and Bug Report
Please raise an issue on [GitHub](https://github.com/Shao-Group/amaranth/issues).

