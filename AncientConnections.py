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

# App Title
st.markdown("<h1 style='text-align: left;'>🧬 Ancient Connections</h1>", unsafe_allow_html=True)
st.markdown("""**Trace your maternal ancestry through DNA connections & find ancient shared ancestors in an interactive timeline.**""")

# App Layout (2 columns)
col1, col2 = st.columns([1, 2])

# Order Epochs for selection
epoch_order = [
    "Modern Age", "Middle Ages", "Iron Age", "Bronze Age",
    "Neolithic", "Mesolithic", "Paleolithic"
]

# User Selection Methods
with col1:
    # First Person
    user1_method = st.radio("**Select Person 1**", ["Test User", "Manual Input", "Database Sample"], key="user1_method")
    
    # Select one of the test users
    if user1_method == "Test User":
        user1 = st.selectbox("Select a test user", list(test_users.keys()), key="user1")
        user1_haplo = test_users[user1]["mtDNA"]
        user1_year = 2025
        
    # Manually input user
    elif user1_method == "Manual Input":
        user1_haplo = st.text_input("Enter mtDNA haplogroup", key="user1_haplo")
        user1_year = st.number_input("Enter estimated year (negative for BCE)", value=1000, key="user1_year")
        
    # Choose Database sample
    else:
        # Filter dataset to exclude N/A values in mt_hg column
        valid_data = AADR_data[AADR_data["mt_hg"].notna() & (AADR_data["mt_hg"] != "N/A")]
        continent1 = st.selectbox("Select continent", valid_data["Continent"].unique(), key="continent1")
        country1 = st.selectbox("Select country", valid_data[valid_data["Continent"] == continent1]["Country"].unique(), key="country1")

        # Filter valid epochs for the selected country
        epoch1_list = valid_data[(valid_data["Continent"] == continent1) & (valid_data["Country"] == country1)]["Epoch"].dropna().unique()
        epoch1_sorted = sorted(epoch1_list, key=lambda x: epoch_order.index(x) if x in epoch_order else len(epoch_order))
        epoch1 = st.selectbox("Select epoch", epoch1_sorted, key="epoch1")

        # Select sample only from valid data
        sample1 = st.selectbox("Select sample", valid_data[
            (valid_data["Continent"] == continent1) &
            (valid_data["Country"] == country1) &
            (valid_data["Epoch"] == epoch1)
        ]["Identifier"], key="sample1")

        # Retrieve the mt_hg and Year values
        user1_haplo = valid_data[valid_data["Identifier"] == sample1]["mt_hg"].values[0]
        user1_year = valid_data[valid_data["Identifier"] == sample1]["Year"].values[0]

    # Repeat for Person 2
    user2_method = st.radio("**Select Person 2**", ["Test User", "Manual Input", "Database Sample"], key="user2_method")

    if user2_method == "Test User":
        user2 = st.selectbox("Select a test user", list(test_users.keys()), key="user2")
        user2_haplo = test_users[user2]["mtDNA"]
        user2_year = 2025
    elif user2_method == "Manual Input":
        user2_haplo = st.text_input("Enter mtDNA haplogroup", key="user2_haplo")
        user2_year = st.number_input("Enter estimated year (negative for BCE)", value=-1000, key="user2_year")
    else:
        valid_data = AADR_data[AADR_data["mt_hg"].notna() & (AADR_data["mt_hg"] != "N/A")]
        continent2 = st.selectbox("Select continent", valid_data["Continent"].unique(), key="continent2")
        country2 = st.selectbox("Select country", valid_data[valid_data["Continent"] == continent2]["Country"].unique(), key="country2")
        epoch2_list = valid_data[(valid_data["Continent"] == continent2) & (valid_data["Country"] == country2)]["Epoch"].dropna().unique()
        epoch2_sorted = sorted(epoch2_list, key=lambda x: epoch_order.index(x) if x in epoch_order else len(epoch_order))
        epoch2 = st.selectbox("Select epoch", epoch2_sorted, key="epoch2")

        sample2 = st.selectbox("Select sample", valid_data[
            (valid_data["Continent"] == continent2) &
            (valid_data["Country"] == country2) &
            (valid_data["Epoch"] == epoch2)
        ]["Identifier"], key="sample2")

        user2_haplo = valid_data[valid_data["Identifier"] == sample2]["mt_hg"].values[0]
        user2_year = valid_data[valid_data["Identifier"] == sample2]["Year"].values[0]

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