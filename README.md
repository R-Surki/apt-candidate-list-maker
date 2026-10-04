# APT Candidate List Maker

A JupyterLab-based interactive tool for making and exploring candidate-ion lists for atom-probe tomography (APT) mass spectra.

The tool reads EPOS data through [APAV](https://gitlab.com/jesseds/apav), generates candidate atomic and molecular ions, displays their expected mass-to-charge positions, shows natural-isotope fingerprints, and exports the candidate list as a text file.

The software was developed primarily for analyzing Fe–C–Mn–Al–Si-based alloys. The current implementation also allows O, H, Ga, and Nb candidates, together with their supported atomic and molecular-ion families.

## Project relationship to APAV

This project is an independent JupyterLab-based candidate list maker built on top of APAV. It does not replace APAV and is not an official APAV interface.

APAV is intended to be used as a Python library. The APAV documentation explains that a GUI may eventually exist alongside APAV, but not as a replacement, and that APAV already includes interactive plotting tools for selected computations. This project explores that complementary GUI direction for candidate-ion selection and mass-spectrum visualization.

Please cite and respect the APAV project and its license when using or redistributing this work:

- APAV repository: <https://gitlab.com/jesseds/apav>
- APAV documentation: <https://apav.readthedocs.io/>
- APAV publication: Smith and Young, *APAV: An Open-Source Python Package for Mass Spectrum Analysis in Atom Probe Tomography*, JOSS 8(83), 4862. <https://doi.org/10.21105/joss.04862>

## Features

- Load an EPOS file once and reuse it while changing candidate selections.
- Select elements and maximum atomic charge states.
- Generate atomic ions and molecular-ion candidates.
- Generate hydrogen molecular ions, carbon clusters, oxygen clusters, carbon-hydrogen ions, hydrides, carbides, oxides, hydroxides, and mixed-metal ions.
- Display isotope fingerprints beside each molecular-ion checkbox.
- Display each selected ion family in a distinct color.
- Add manual candidate peaks.
- Export an interactive Plotly HTML spectrum.
- Export the generated candidate list as a tab-separated `.txt` file.

## Current element scope

The default GUI focuses on:

```text
Fe, C, Mn, Si, Al, O, H, Ga, Nb
```

N, P, S, Cr, Ni, and B are intentionally excluded from the current candidate lists.

## Installation

Create and activate a virtual environment, then install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

## Run in JupyterLab

Start JupyterLab:

```bash
jupyter lab
```

Run the two Python cells in this order:

```python
%run cell1_core.py
```

```python
%run cell2_gui.py
```

Enter the path to your own EPOS file in the GUI. EPOS files are intentionally excluded from this repository.

## Files

- `cell1_core.py` — isotope calculations, candidate generation, and plotting.
- `cell2_gui.py` — Jupyter widgets, EPOS caching, candidate export, and GUI callback.
- `requirements.txt` — Python dependencies.

## Status

This is an actively developing research and visualization tool. Candidate-ion generation is intended to support interpretation and exploration; candidate matches should be checked against the experimental spectrum and the relevant APT literature.

## Acknowledgment and citation

If this software contributes to a publication, presentation, thesis, or other research output, users are kindly requested to acknowledge the author and cite the project repository.

Suggested acknowledgment:

> The authors acknowledge R-Surki for developing the APT Candidate List Maker, an open-source JupyterLab tool for candidate-ion generation and mass-spectrum visualization based on APAV.

Suggested software citation:

> R-Surki. *APT Candidate List Maker*. GitHub repository: <https://github.com/R-Surki/apt-candidate-list-maker>. Accessed [DATE].

Please also cite APAV separately when it is used in the analysis.

## License and attribution

This project is an independent companion tool that uses APAV. Before publishing a final repository license, verify compatibility with APAV's license and decide which license you want for your own additions. Keep the APAV attribution and license notices intact.
