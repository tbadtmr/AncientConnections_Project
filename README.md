# 🧬 Ancient Connections

**Ancient Connections** is an interactive Streamlit app for exploring shared maternal lineages between modern users and ancient DNA samples. 
It visualizes genetic ancestry through interactive timelines, haplogroup trees, and a summary in an easy-to-use interface.

---

## Features

- Compare two users via mtDNA haplogroups (manual, database, or test inputs)
- Detect the last shared ancient maternal ancestor in the database
- Trace shared maternal lineage over time
- Visualize ancient samples on an interactive timeline and display lineages
---

## Installation

This app requires the following Python packages:

- `streamlit >= 1.43`
- `pandas >= 2.2`
- `numpy >= 1.24`
- `plotly >= 6.0`

To ensure full compatibility and reproducibility, it is recommended to set up a new Conda environment:

```bash
conda create -n ancient_connections python=3.9
conda activate ancient_connections
pip install streamlit pandas numpy plotly
```

## How to run the App

Once the environment is set, the app can be run with:
```bash
streamlit run AncientConnections.py
```
Then open the app in your browser at the local address provided by Streamlit.

## Data

This app uses the following datasets:
- Ancient mtDNA samples — curated from the AADR (Allen Ancient DNA Resource, 2025)
- mtDNA phylogenetic tree — to resolve maternal lineages
- Test users — 5 sample users with mtDNA haplogroups

The ancient DNA data was preprocessed using the script `AADR_dataParser.py` and manually cleaned for consistency.
All required data is included in the `Data/` folder.

## Project Repository

You can find the full project, source code, and documentation at:
  - GitHub: tbadtmr/AncientConnections_Project

--

**Have Fun Exploring!**
Whether you're tracing your own maternal line or diving into the deep ancestry of ancient populations —
enjoy exploring these ancient connections :)
