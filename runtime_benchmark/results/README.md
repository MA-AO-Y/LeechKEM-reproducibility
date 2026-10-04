# Generated results

Running `../run_full_benchmark.sh` replaces or creates the following files:

- `leech_full.csv`
- `frodo640_full.csv`, `frodo976_full.csv`, `frodo1344_full.csv`
- `runtime_results.csv`: all raw summary rows in one file, including the
  independent case-pool size and full-cycle timing
- `runtime_results.md`: formatted table and paired ratios
- `runtime_results.tex`: ready-to-paste runtime table
- `runtime_ratios.tex`: ready-to-paste paired-ratio table
- `environment.txt`: CPU, OS/compiler, flags, and benchmark settings

Do not use any sample results from another computer in the paper.  Run the
default command on the reporting machine and retain both `runtime_results.csv`
and `environment.txt`.
