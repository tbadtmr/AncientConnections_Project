#%% 
# DOCUMENTATION SECTION - Ancient Connections App
#______________________________________________________________________________
#!/usr/bin/python3
"""
Description:

Procedure:

Usage:


Version: 1.00
Date: 2025-03-13
Name: Tabea Dittmar
"""

#%%
import streamlit as st
import pandas as pd
import numpy as np
import re
import plotly.graph_objects as go

#%%
# IMPORT DATA
#______________________________________________________________________________
# Set Page Configuration
st.set_page_config(page_title="Ancient Connections", layout="wide")

# Load Data
@st.cache_data
def load_data():
    # Ancient database samples
    mtDNA_data = pd.read_excel("Data/mtdb_metadata.xlsx")
    # Cleaned user selection database
    AADR_data = pd.read_excel("Data/AADR_clean.xlsx")
    
    # Read phylogenetic tree
    phylo_tree = {"mt-MRCA": None}
    with open("Data/mt_phyloTree_b17_Tree2.txt", "r") as f:
        next(f)  # Skip first line
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                child, parent = parts
                phylo_tree[child] = parent
            if parent not in phylo_tree and parent != "mt-MRCA":
                phylo_tree[parent] = "mt-MRCA"

    return mtDNA_data, AADR_data, phylo_tree

# Load Data
mtDNA_data, AADR_data, phylo_tree = load_data()

# Test Users & Haplogroups
test_users = {
    "User1": {"mtDNA": "X2c2", "Y-DNA": "R-S4634"},
    "User2": {"mtDNA": "K1c1b", "Y-DNA": "I-FT90548"},
    "User3": {"mtDNA": "U2", "Y-DNA": "R-M269"},
    "User4": {"mtDNA": "V7", "Y-DNA": "FGC11293"},
    "User5": {"mtDNA": "U5b2c", "Y-DNA": "R-CTS11962"},
}

#%%
# CHOOSE SAMPLE INPUT METHOD
#______________________________________________________________________________

#%%
# GET MOST RECENT HAPLOGROUP AND LINEAGE
#______________________________________________________________________________


#%%
# FIND ANCESTORS IN DATABANK
#______________________________________________________________________________


# %%
# PLOT
#______________________________________________________________________________


#%%
# DISPLAY EVERYTHING
#______________________________________________________________________________