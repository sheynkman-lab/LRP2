# Support and Customization

LRP2 supports a range of customization:

- **HPC environment**: SLURM and LSF schedulers with configurable partition/queue and cluster options; Singularity, Docker, or Conda containers
- **Input flexibility**: RNA-only, protein-only, or paired RNA + protein samples; DDA or DIA mass spec data
- **Reference genomes**: GENCODE human and mouse across multiple versions
- **Proteomics**: alternative search engines (FragPipe or MetaMorpheus) and customizable FragPipe workflows

We welcome input from the community. Please open an Issue if you have a use case not covered by the defaults. 

## Getting Help

### Issues and Bug Reports

For bug reports, feature requests, or general issues, please open an **Issue**: [github.com/sheynkman-lab/LRP2/issues](https://github.com/sheynkman-lab/LRP2/issues)

When reporting an issue, please attach the following files to your Issue:
1. Your run logfile (`.nextflow.log`)
2. Your `samplesheet.csv` file for the run 

Additionally, please share the following details in the body of your Issue: 

1. Your command line and configuration
2. Nextflow version and executor (local, SLURM, LSF)
3. A brief overview of your data (number of samples, read depth / file sizes, etc.)
4. Your exact error message(s) and any other relevant log files (e.g. any other relevant logs for individual processes from the `work/` directory)

### Direct Contact

For direct inquiries or collaboration opportunities, please reach out directly:

- **Megan Schertzer**: cwp5au@virginia.edu
- **Julia Lewandowski**: jlewandowski@nygenome.org

## Community Support

Questions about:

- **General usage**: Check the [Usage Guide](usage/index.md) and [Parameters](usage/parameters.md)
- **Output interpretation**: See [Output Files](output/structure.md)
- **Installation issues**: Review [Installation](getting-started/installation.md)
- **Reference genomes**: See [Genomes](reference/genomes.md)

## Contributing

We welcome contributions from the community. We are currently building resources to faciliate this process. For now, if you would like to contribute, please open an issue directory and feel free to open a PR to develop your feature.

## Non-Academic Use

For non-academic or commercial use cases, particularly regarding FragPipe licensing, please contact the development team for guidance.
