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


def get_lineage(haplogroup, tree):
    """
    Returns the ancestral path of a haplogroup up to the root (mt-MRCA).
    If a haplogroup is missing, it finds the closest known parent.
    """
    lineage = []

    # If haplogroup is missing, try to infer its closest known ancestor
    if haplogroup not in tree:
        haplogroup = infer_parent(haplogroup, tree) or "mt-MRCA"

    # Traverse the tree upwards
    while haplogroup in tree and tree[haplogroup] is not None:
        lineage.append(haplogroup)
        haplogroup = tree[haplogroup]

    lineage.append("mt-MRCA")  # Ensure root is included
    return lineage[::-1]  # Reverse to start from the root

def infer_parent(haplogroup, tree):
    """
    Tries to infer the most likely parent haplogroup if it's missing in the tree.
    It progressively removes numbers from the haplogroup name (e.g., R1a2 → R1a).
    """
    while haplogroup:
        haplogroup = re.sub(r'\d+$', '', haplogroup)  # Remove trailing numbers
        if haplogroup in tree:
            return haplogroup  # Found a valid parent in the tree
        if len(haplogroup) <= 1:  # Prevent guessing too far back
            break
    return None  # No valid parent found

def find_common_mtDNA(haplo1, haplo2, tree):
    """
    Finds the Most Recent Common Ancestor (MRCA) of two haplogroups.
    """
    lineage1 = get_lineage(haplo1, tree)
    lineage2 = get_lineage(haplo2, tree)

    # Find the last common ancestor
    common_ancestor = "mt-MRCA"
    for h1, h2 in zip(lineage1, lineage2):
        if h1 == h2:
            common_ancestor = h1
        else:
            break
    return common_ancestor


# Find Most Recent Common Ancestor Haplogroup
#----------------------------------------------
common_mtDNA = find_common_mtDNA(user1_haplo, user2_haplo, phylo_tree)

# Get lineage
#---------------------------------------------
# Compute Shared Lineage (From mt-MRCA to Shared Haplogroup)**
shared_lineage = get_lineage(common_mtDNA, phylo_tree)
shared_lineage_set = set(shared_lineage)  # ✅ Convert to set for faster lookups

# Compute Individual Lineages
user1_lineage = get_lineage(user1_haplo, phylo_tree)
user2_lineage = get_lineage(user2_haplo, phylo_tree)

# Extract the Part of the Lineage That Differs
if common_mtDNA in user1_lineage:
    lineage_to_user1 = user1_lineage[user1_lineage.index(common_mtDNA):]  # From MRCA → User 1
else:
    lineage_to_user1 = ["Not found"]
if common_mtDNA in user2_lineage:
    lineage_to_user2 = user2_lineage[user2_lineage.index(common_mtDNA):]  # From MRCA → User 2
else:
    lineage_to_user2 = ["Not found"]


#%%
# FIND ANCESTORS IN DATABANK
#______________________________________________________________________________

most_recent_identifier, most_recent_haplogroup, most_recent_year, most_recent_year2, most_recent_country = None, None, None, None, None
most_ancient_identifier, most_ancient_haplogroup, most_ancient_year, most_ancient_year2, most_ancient_country = None, None, None, None, None
traceback_samples = []
related_samples = []
show_recent = False
show_ancient = False
show_traceback = False
show_relatives = False
show_common_mtDNA = False

# Get all ancient samples from the lineage
lineage_samples = mtDNA_data[
    # Must be in shared lineage
    (mtDNA_data["mt_hg"].isin(shared_lineage_set)) &
    # Must be older than both users
    (mtDNA_data["year_from"] <= min(user1_year, user2_year))
]

# Get most recent connection (EXACT haplogroup match)
if not lineage_samples.empty:
    recent_ancestor = lineage_samples[lineage_samples["mt_hg"] == common_mtDNA]
    if not recent_ancestor.empty:
        recent_ancestor = recent_ancestor.sort_values(by="year_from", ascending=False).iloc[0]  # Youngest match
        most_recent_identifier = recent_ancestor["identifier"]
        most_recent_haplogroup = recent_ancestor["mt_hg"]
        most_recent_year = recent_ancestor["year_from"]
        most_recent_year2 = recent_ancestor["year_to"]
        most_recent_country = recent_ancestor["country"]
        most_recent_sex = recent_ancestor["sex"]
        most_recent_region = recent_ancestor["region"]
        most_recent_epoch = recent_ancestor["epoch"]
        show_recent = True

# Get most ancient connection (EXACT haplogroup match)
if not lineage_samples.empty:
    ancient_ancestor = lineage_samples[lineage_samples["mt_hg"] == common_mtDNA]
    if not ancient_ancestor.empty:
        ancient_ancestor = ancient_ancestor.sort_values(by="year_from", ascending=True).iloc[0]  # Oldest match
        most_ancient_identifier = ancient_ancestor["identifier"]
        most_ancient_haplogroup = ancient_ancestor["mt_hg"]
        most_ancient_year = ancient_ancestor["year_from"]
        most_ancient_year2 = ancient_ancestor["year_to"]
        most_ancient_country = ancient_ancestor["country"]
        most_ancient_sex = ancient_ancestor["sex"]
        most_ancient_region = ancient_ancestor["region"]
        most_ancient_epoch = ancient_ancestor["epoch"]
        show_ancient = True

# Check if both are the same sample, if so hide one
if (show_recent and show_ancient and 
    most_recent_identifier == most_ancient_identifier):
    show_recent = False

if show_ancient == False:
    show_common_mtDNA = True
    
# Check if there more ancient samples from the lineage
if show_ancient:
    # Only look for ancestors older than the most ancient ancestor
    lineage_samples["year_from"] = pd.to_numeric(lineage_samples["year_from"], errors="coerce")
    additional_ancestors = lineage_samples[
        lineage_samples["year_from"] < most_ancient_year
        ]
else:
    # If no ancient ancestor was found, include all lineage samples
    additional_ancestors = lineage_samples.copy()

# Oder additional ancestors
additional_ancestors["lineage_depth"] = additional_ancestors["mt_hg"].apply(
    lambda hg: shared_lineage.index(hg) if hg in shared_lineage else len(shared_lineage)
)
# Sort first by lineage (following haplogroup tree), then by year
additional_ancestors = additional_ancestors.sort_values(by=["lineage_depth", "year_from"], ascending=[True, True])
# Convert to list format
additional_ancestors_list = additional_ancestors.to_dict(orient="records")

# %%
# PLOT
#______________________________________________________________________________


#%%
# DISPLAY EVERYTHING
#______________________________________________________________________________