# Water Counting Script for UCSF Chimera

## Dependencies
* **Software:** UCSF Chimera 1.19 (Legacy Production Release, March 2025)
* **Python:** Python 2.7 (Bundled with Chimera)
* **Dependencies:** Standard Chimera modules (`chimera`, `Molecule`, `runCommand`)
  

## Usage
Run this script from the Chimera Command Line:
`open script-plane_water.py`

## Developer Notes (For AI Agents)
* **DO NOT upgrade to Python 3.** This script must remain compatible with Python 2.7.
* Strictly use `print` statements compatible with Py2.7.
* Do not use f-strings (e.g., `f"{var}"`); use `%` formatting or `.format()`.

## Citation

If you use **Water Counting** in your research, please cite it as follows:

> Peralta-Morales, M. F., Sciortino, G. & Maréchal, J.-D. Water Counting. *GitHub* https://github.com/insilichem/Water-count-per-plane (2026).

### BibTeX
```bibtex
@misc{PeraltaMorales2026WaterCounting,
  author       = {Peralta-Morales, Mar{\'i}a Fernanda and Sciortino, Giuseppe and Mar{\'e}chal, Jean-Didier},
  title        = {Water Counting: Counts water molecules per plane},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{[https://github.com/insilichem/Water-count-per-plane](https://github.com/insilichem/Water-count-per-plane)}}
}
