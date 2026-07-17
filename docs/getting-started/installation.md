# Installation

## Prerequisites

- **Nextflow** ≥ 24.04.2 ([installation guide](https://www.nextflow.io/docs/latest/install.html))
- **Singularity/Apptainer** or **Docker** for containerized dependencies

## Containers

LRP2 uses containers to manage software dependencies. Containers allow for the packaging of libraries, code, and configurations such that each tool in the pipeline can be reliably run in any computing environment without compatibility issues.

To run LRP2, you **must have one** of the following installed:

### 1. Singularity/Apptainer (required for HPC)

Most HPC systems have Singularity or Apptainer pre-installed as a module. You can check this by running `module avail singularity` or `module avail apptainer`. If not available, contact your HPC administrator or see [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html).

### 2. Docker (may be used for local systems)

Installation guides: [Docker Desktop](https://docs.docker.com/get-docker/) (Mac/Windows) or [Docker Engine](https://docs.docker.com/engine/install/) (Linux)

!!! warning "Docker on HPC"
    Docker requires root/admin privileges, and is typically not permitted on shared HPC systems. We therefore strongly recommend the use of Singularity/Apptainer.

The pipeline will automatically pull and cache container images on first run. Note that Singularity images are cached in `work/singularity/` by default. 

### 3. Conda (not recommended)

While LRP2 technically supports Conda environments as a fallback, **we strongly discourage using Conda** due to significant reproducibility limitations. 

If you are interested in learning more about best practices for computational reproducibility, please refer to [this](https://pmc.ncbi.nlm.nih.gov/articles/PMC6263957/) paper. 

!!! danger "Why shouldn't I use Conda?"
    In short, **containers provide superior reproducibility** because they freeze the entire software stack including the operating system, tools, and all dependencies at a specific point in time. This ensures analyses run identically across different systems and over many years. **Conda environments**, by contrast, are resolved dynamically at runtime based on the package ecosystem at that moment. This means fresh Conda builds may pull newer versions of dependencies (including critical components like R or Python), which can cause results to differ when reruns are done months or years later. 

For these reasons, please always use Singularity/Apptainer or Docker for analyses unless you have a specific technical constraint that prevents container usage.

## Clone the Repository

```bash
git clone https://github.com/sheynkman-lab/LRP2.git
cd LRP2
```

## Next Steps

Once you have installed the prerequisites and cloned the repository, proceed to the [Quick Start](quickstart.md) guide to run the test dataset.
