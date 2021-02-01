# cohort_sv

A set of cohort plotting script for the SOPHIA SV caller.

This script parses output from multiple SOPHIA structural variation "minEventScore3" files, and create structures that can be plotted circlize

Whie this was developed for the SOPHIA SV caller, technically it should work with any set of SV calls in BED-PE format (where the first 3 cols are the start of the SV, and the next 3 columns and the end of the SV)

This script expects a folder sctructure created by the softlink project script created by @ishaque

The PID list and or file list can be forced using --pid_list and --file_list options

## Prerequisites

```
perl v5.20.0
bedtools 2.24.0
r v3.3.1
 - library(circlize)
 - xlibrary(dplyr)
```

## Conda

`conda env create cohort_sv perl=5.20.0 bedtools=2.24.0 r=3.3.1` 

## Modules (only on Eils and ODCF clusters)

```
module load perl/5.20.0
module load bedtools/2.24.0
module load R/3.3.1
```

## Usage

First the SV calls for the cohort have to be preprocessed

```
  perl process_SV_for_circos.pl
		--results_dir     -r   [path: file path to the results_per_pid folder. REQUIRED!]
		--output_prefix   -o   [string: output prefix. REQUIRED!]
		--sv_close        -s   [integer: distance to define close SV. Default 10000]
		--deldup_window   -d   [integer: size of del/dup windows for aggregation. Default 1000000]
		--pids            -p   [file or list: file containing a list of PIDs, or a space separated list of pids]
		--file_list       -f   [list of files to use... over rides PIDs]
		--help            -h   [print usage information]
```

Example call: `perl process_SV_for_circos.pl -r /path/tp/hipo_050/whole_genome_sequencing/results_per_pid -o test_H050 -p "H050-ABCD H050-EFGH H050-IJKL" `

This will create a TSV files which can be used by a subsequent script for plotting. An example call which shows all arcs targetting genes at least `2` times in grey, and at least `6` times in red.

```
R -f make_cohort_SV_circlize_plot.R --args [output from process_SV_for_circos] 6 2 "My title"
```

This will create a PDF of 3 pages, describing genes that are (i) hit directly, (ii) close to a gene, i.e defined by the `sv_close` parameter in `process_SV_for_circos.pl`, and (iii) to the nearest gene.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).
 - v1.0.0: first working version

## Authors

Naveed Ishaque

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgements
