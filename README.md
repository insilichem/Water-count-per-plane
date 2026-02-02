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
