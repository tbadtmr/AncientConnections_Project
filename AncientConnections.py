#!/usr/bin/python3
#%% 
# DOCUMENTATION SECTION - Ancient Connections App
#______________________________________________________________________________

"""
Description:
Ancient Connections is an interactive Streamlit application for tracing shared maternal ancestry 
between modern users and ancient individuals based on mitochondrial DNA (mtDNA) haplogroups. 
It identifies the most recent common maternal ancestor and visualizes related ancient samples 
on a timeline with lineage paths and historical context.

Procedure:
1. Load user input (test user, manual entry, or ancient sample).
2. Determine the most recent common haplogroup using a phylogenetic tree.
3. Search for ancient samples with matching or ancestral haplogroups.
4. Visualize connections on an interactive timeline with lineage paths.

Usage:
Run the app locally with:
    streamlit run AncientConnections.py

User-defined functions:
- get_lineage(haplogroup, tree): 
    Returns the full maternal lineage for a given haplogroup using the mtDNA phylogenetic tree.
- infer_parent(haplogroup, tree): 
    Heuristically infers the closest ancestral haplogroup if a specific node is missing from the tree.
- find_common_mtDNA(haplo1, haplo2, tree): 
    Identifies the most recent common maternal haplogroup shared between two users.
- bezier_curve(x_start, x_end, y_start, y_end, ...): 
    Generates smooth Bezier curves to visually connect samples along the timeline.
- format_year(year):
    restructures continous timescale back to historical annotation (BCE, CE)

Version: 2.00  
Date: 2025-03-24  
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
# DEFINE FUNCTIONS
#______________________________________________________________________________

def get_lineage(haplogroup, tree):
    """
    Returns the ancestral path of a haplogroup up to the root (mt-MRCA).
    If a haplogroup is missing, it finds the closest known parent.
    """
    lineage = []

    # If haplogroup is missing in the tree, try to infer its closest known ancestor
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
    Finds the most recent shared haplogroup of two haplogroups.
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
        user1_haplo = None
        user1_haplo_input = st.text_input("Enter mtDNA haplogroup", key="user1_haplo")
        user1_year = st.number_input("Enter estimated year (negative for BCE)", value=1000, key="user1_year")

        if user1_haplo_input:
            haplo = user1_haplo_input.strip().upper()
            # Validate haplogroup entry (has to be in the tree)
            if haplo in phylo_tree:
                user1_haplo = haplo
            # Otherwise print warning
            else:
                st.warning("Please enter a valid Haplogroup.")

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

        # Retrieve the mt_hg and year values
        user1_haplo = valid_data[valid_data["Identifier"] == sample1]["mt_hg"].values[0]
        user1_year = valid_data[valid_data["Identifier"] == sample1]["Year"].values[0]

    # Repeat for Person 2
    user2_method = st.radio("**Select Person 2**", ["Test User", "Manual Input", "Database Sample"], key="user2_method")

    if user2_method == "Test User":
        user2 = st.selectbox("Select a test user", list(test_users.keys()), key="user2")
        user2_haplo = test_users[user2]["mtDNA"]
        user2_year = 2025
    elif user2_method == "Manual Input":
        user2_haplo = None
        user2_haplo_input = st.text_input("Enter mtDNA haplogroup", key="user2_haplo")
        user2_year = st.number_input("Enter estimated year (negative for BCE)", value=1000, key="user2_year")

        if user2_haplo_input:
            haplo = user2_haplo_input.strip().upper()
            if haplo in phylo_tree:
                user2_haplo = haplo
            else:
                st.warning("Please enter a valid Haplogroup.")
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
        
# CHECK: Make sure two different samples are selected
# ---------------------------------------------------
# Only proceed if both users have valid haplogroups
if user1_haplo and user2_haplo:
    # Case: Identical haplogroup and year — likely the same person
    if (user1_haplo == user2_haplo) and (user1_year == user2_year):
        skip_comparison = True
        identity_warning = "Person 1 and Person 2 are identical. Please select different samples."
    else:
        # Valid comparison between two distinct users
        skip_comparison = False
        identity_warning = ""
else:
    # One or both haplogroups missing — skip comparison
    skip_comparison = True
    identity_warning = ""

#%%
# GET MOST RECENT HAPLOGROUP AND LINEAGE
#______________________________________________________________________________

# Find Most Recent COMMON HAPLOGROUP
common_mtDNA = find_common_mtDNA(user1_haplo, user2_haplo, phylo_tree)

# GET LINEAGE
#------------
# Compute Shared Lineage (From mt-MRCA to Shared Haplogroup)
shared_lineage = get_lineage(common_mtDNA, phylo_tree)
shared_lineage_set = set(shared_lineage)

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

# Get all ancient samples from the lineage
lineage_samples = mtDNA_data[
    (mtDNA_data["mt_hg"].isin(shared_lineage_set)) &
    # have to be younger than users
    (mtDNA_data["year_from"] <= min(user1_year, user2_year))
].copy()

# Assign lineage depth
lineage_samples["lineage_depth"] = lineage_samples["mt_hg"].apply(
    lambda hg: shared_lineage.index(hg) if hg in shared_lineage else len(shared_lineage)
)

# Sort by lineage depth and year
lineage_samples = lineage_samples.sort_values(by=["lineage_depth", "year_from"])

# Convert to list for plotting
lineage_samples_list = lineage_samples.to_dict(orient="records")


# Filter lineage samples to preserve lineage + time order
# -------------------------------------------------------
# Goal: Select a clean subset of ancient samples along the shared maternal lineage,
# keeping only one representative per haplogroup (or closest in time),
# while ensuring samples get progressively more ancient as you move up the lineage.

filtered_lineage = []
last_depth = -1  # Track position in the shared haplogroup lineage
last_year = float('-inf')  # Track most recent sample added
last_hg = None  # Track last haplogroup added

for sample in lineage_samples_list:
    hg = sample["mt_hg"]
    year = sample["year_from"]

    # Get the haplogroup's depth in the shared lineage
    depth = shared_lineage.index(hg) if hg in shared_lineage else None

    # Skip if not in the defined shared lineage path
    if depth is None:
        continue

    # Always accept the first valid sample
    if not filtered_lineage:
        filtered_lineage.append(sample)
        last_depth = depth
        last_year = year
        last_hg = hg
        continue

    # Allow:
    # - same haplogroup as previous, but more recent
    # - new haplogroup deeper in lineage and more recent
    if (hg == last_hg and year > last_year) or (depth > last_depth and year > last_year):
        filtered_lineage.append(sample)
        last_depth = depth
        last_year = year
        last_hg = hg

# Replace lineage_samples_list with the filtered version
lineage_samples_list = filtered_lineage

# Define whether to show the common haplogroup as a node
show_common_mtDNA = not any(ancestor["mt_hg"] == common_mtDNA for ancestor in lineage_samples_list)


# %%
# PLOT
#______________________________________________________________________________

# Initialize timeline data
labels = ["User 1", "User 2"]
years = [user1_year, user2_year]
colors = ["orange", "purple"]
hover_texts = [f"User 1: {user1_haplo}", f"User 2: {user2_haplo}"]

# Dynamic Y spacing
ancestor_step = 0.3  # space between ancestors
user_gap = 0.6       # space between youngest ancestor and users

# How many ancestors?
num_ancestors = len(lineage_samples_list)

# Y positions (ancestors from top down)
y_positions = [(num_ancestors - i) * ancestor_step for i in range(num_ancestors)]

# Timeline lists
labels = []
years = []
colors = []
hover_texts = []

# ADD ALL ANCESTORS TO THE PLOT
# ------------------------------
# Loop through each ancient sample in the shared lineage and add it to the plot
for i, ancestor in enumerate(lineage_samples_list):
    labels.append(ancestor["identifier"])                # Sample ID
    years.append(ancestor["year_from"])                  # Date for X-axis
    # Highlight the one that exactly matches the shared haplogroup
    colors.append("darkblue" if ancestor["mt_hg"] == common_mtDNA else "lightblue")
    
    # Detailed hover info per sample
    hover_texts.append(
        f"<b>Ancestor</b>: {ancestor['identifier']}<br>"
        f"<b>Haplogroup</b>: {ancestor['mt_hg']}<br>"
        f"<b>Country</b>: {ancestor['country']}<br>"
        f"<b>Region</b>: {ancestor['region']}<br>"
        f"<b>Lifespan</b>: {ancestor['year_from']} to {ancestor['year_to']}"
    )

# Add a synthetic node for the shared haplogroup if no direct ancient sample exists
if not any(ancestor["mt_hg"] == common_mtDNA for ancestor in lineage_samples_list):
    # Estimate position further back in time
    shared_haplogroup_year = min(user1_year, user2_year) - 500
    # Estimate Y position just above users
    shared_haplogroup_pos = (y_positions[0] + ancestor_step) if y_positions else 1.0

    # Insert shared haplogroup node at the beginning of the timeline
    labels.insert(0, common_mtDNA)
    years.insert(0, shared_haplogroup_year)
    y_positions.insert(0, shared_haplogroup_pos)
    colors.insert(0, "black")
    hover_texts.insert(0, (
        f"<b>Shared Haplogroup</b>: {common_mtDNA}<br>"
        "This is the most recent common haplogroup shared by both users. No ancient database sample found."
    ))
    
# ADD USERS AT THE BOTTOM
# -------------------------
# Position users at the bottom of the plot, below all ancestors
user1_y = -user_gap                    # Base Y position for User 1
user2_y = user1_y - 0.4               # Slight vertical offset to separate User 2 visually

# Add users to the plot data
labels.extend(["User 1", "User 2"])   # Display names
years.extend([user1_year, user2_year])  # X-axis values (time)
y_positions.extend([user1_y, user2_y])  # Y-axis placement
colors.extend(["orange", "purple"])     # User-specific colors

# Add basic hover text with haplogroup info
hover_texts.extend([
    f"User 1: {user1_haplo}",
    f"User 2: {user2_haplo}"
])

# CREATE FIGURE
#--------------
fig = go.Figure()

# Add Scatter Points
fig.add_trace(go.Scatter(
    x=years, 
    y=y_positions, 
    mode="markers+text", 
    marker=dict(size=15, color=colors, line=dict(width=2, color="black")), 
    text=labels,
    textposition="top center",
    hovertext=hover_texts,
    hoverinfo="text"
))

# ADD CONNECTION LINES
#----------------------
# Bezier Curve Function for Smooth Connections
def bezier_curve(x_start, x_end, y_start, y_end, curve_direction="up", num_points=100, stop_offset=5):
    """ Generate curved Bezier path for smooth connections. """
    t = np.linspace(0, 1, num_points)
    # Dynamic Control Point:
    control_y = (y_start + y_end) / 2  # Default middle control point
    control_x = (x_start + x_end) / 2  # Midpoint in X
    # Push Control Point Higher or Lower to Strengthen the Curve
    vertical_distance = abs(y_end - y_start)
    curve_strength = vertical_distance * 0.5  # Increase strength based on spacing
    if curve_direction == "up":
        control_y += curve_strength  # Push control point UP
    else:
        control_y -= curve_strength  # Push control point DOWN
    # Bezier Curve Calculation
    bezier_x = (1 - t) ** 2 * x_start + 2 * (1 - t) * t * control_x + t ** 2 * x_end
    bezier_y = (1 - t) ** 2 * y_start + 2 * (1 - t) * t * control_y + t ** 2 * y_end
    # Apply Offset to Avoid Sharp Stops
    bezier_x[-1] = bezier_x[-1] - stop_offset  

    return bezier_x, bezier_y

# ADD SHARED HAPLOGROUP NODE (IF NOT ALREADY IN DATA)
# -----------------------------------------------------
if show_common_mtDNA:
    # Place shared haplogroup before all samples (a synthetic anchor point)
    shared_haplogroup_year = min(user1_year, user2_year) - 500  # Push back in time
    shared_haplogroup_pos = (y_positions[0] + y_positions[-1]) / 2  # Center between users

    labels.insert(0, common_mtDNA)
    years.insert(0, shared_haplogroup_year)
    y_positions.insert(0, shared_haplogroup_pos)
    colors.insert(0, "black")
    hover_texts.insert(0, (
        f"<b>Shared Haplogroup</b>: {common_mtDNA}<br>"
        "This is the most recent common haplogroup shared by both users."
    ))

# Build a lookup for each label's Y position
y_position_map = {label: y for label, y in zip(labels, y_positions)}

# 1. Connect Users → Most Recent Ancestor
#      (only if shared haplogroup dot not shown)
# ------------------------------------------------------------------
if lineage_samples_list and not show_common_mtDNA:
    youngest = lineage_samples_list[-1]  # last in list is most recent
    if youngest["identifier"] in y_position_map:
        for user_label, user_year in zip(["User 1", "User 2"], [user1_year, user2_year]):
            x_vals, y_vals = bezier_curve(
                user_year, youngest["year_from"],
                y_position_map[user_label], y_position_map[youngest["identifier"]],
                "up"
            )
            fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="lightblue", width=2)))

# 2. Connect Ancestors to Each Other (in order)
# ------------------------------------------------------------------
for i in range(len(lineage_samples_list) - 1):
    older = lineage_samples_list[i]
    younger = lineage_samples_list[i + 1]

    if older["identifier"] in y_position_map and younger["identifier"] in y_position_map:
        x_vals, y_vals = bezier_curve(
            older["year_from"], younger["year_from"],
            y_position_map[older["identifier"]], y_position_map[younger["identifier"]],
            "down"
        )
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="lightblue", width=2)))

# 3. If shared haplogroup was added (i.e., not found in data)
# ------------------------------------------------------------------
if show_common_mtDNA and common_mtDNA in y_position_map:

    # Connect users → shared haplogroup
    for user_label, user_year in zip(["User 1", "User 2"], [user1_year, user2_year]):
        x_vals, y_vals = bezier_curve(
            user_year, shared_haplogroup_year,
            y_position_map[user_label], y_position_map[common_mtDNA],
            "up"
        )
        fig.add_trace(go.Scatter(
            x=x_vals, y=y_vals, mode="lines",
            line=dict(color="grey", width=2, dash="dash")
        ))

    # Connect shared haplogroup → oldest ancestor
    if lineage_samples_list:
        oldest = lineage_samples_list[0]  # first in list
        if oldest["identifier"] in y_position_map:
            x_vals, y_vals = bezier_curve(
                shared_haplogroup_year, oldest["year_from"],
                y_position_map[common_mtDNA], y_position_map[oldest["identifier"]],
                "down"
            )
            fig.add_trace(go.Scatter(
                x=x_vals, y=y_vals, mode="lines",
                line=dict(color="grey", width=2, dash="dash")
            ))
            
#%%
# PLOT VISUALS
#______________________________________________________________________________           

# Epoch scale:
epoch_definitions = {
    "Paleolithic": (-50000, -12000),
    "Mesolithic": (-12000, -8000),
    "Neolithic": (-8000, -3000),
    "Bronze Age": (-3000, -1200),
    "Iron Age": (-1200, 500),
    "Middle Ages": (500, 1500),
    "Modern Age": (1500, 2025)
}

epoch_colors = {
    "Paleolithic": "#6E6E6E",
    "Mesolithic": "#2E8B57",
    "Neolithic": "#BDB76B",
    "Bronze Age": "#8B4513",
    "Iron Age": "#CD5C5C",
    "Middle Ages": "#800080",
    "Modern Age": "#4682B4"
}
# Dynamically Filter Epochs in Range
# Prepare Epoch Data
filtered_epochs = []
epoch_positions = []
epoch_tick_labels = []
min_year, max_year = min(years), max(years)

for epoch, (start, end) in epoch_definitions.items():
    if start <= max_year and end >= min_year:
        filtered_epochs.append(epoch)
        epoch_positions.append((start + end) / 2)
        epoch_tick_labels.append(epoch)

# Define Time Scale (X-Axis) Tick Formatting
tick_positions = np.linspace(min_year, max_year, num=10).astype(int)
tick_labels = [f"{abs(t)} BCE" if t < 0 else f"{t} CE" for t in tick_positions]

# Y-Axis Range — Add dynamic space below lowest Y for epoch bar
y_buffer_bottom = 0.6  # space below users
y_buffer_top = 0.5     # optional space above ancestors
y_min = min(y_positions) - y_buffer_bottom
y_max = max(y_positions) + y_buffer_top

# Update Layout
fig.update_layout(
    xaxis=dict(
        title="Time (Years BCE/CE)",
        showgrid=False,
        tickmode="array",
        tickvals=tick_positions,
        ticktext=tick_labels,
        side="top",
    ),
    xaxis2=dict(
        title="Epochs",
        showgrid=False,
        tickmode="array",
        tickvals=epoch_positions,
        ticktext=epoch_tick_labels,
        anchor="y",
        overlaying="x",
        side="bottom",
    ),
    yaxis=dict(
        visible=False,
        range=[y_min, y_max],
    ),
    showlegend=False,
    margin=dict(l=60, r=60, t=60, b=120),
)

# Draw Epoch Bars and Labels
for epoch, (start, end) in epoch_definitions.items():
    if start <= max_year and end >= min_year:
        x0 = max(start, min_year)
        x1 = min(end, max_year)

        # Background colored bar for each epoch
        fig.add_shape(
            type="rect",
            x0=x0, x1=x1,
            y0=y_min + 0.05, y1=y_min + 0.25,  # below user nodes
            fillcolor=epoch_colors[epoch],
            opacity=0.4,
            layer="below",
            line=dict(width=0)
        )

        # Text label for epoch
        text_x = (x0 + x1) / 2
        fig.add_annotation(
            x=text_x,
            y=y_min + 0.15,
            text=epoch,
            showarrow=False,
            font=dict(size=11, color="black"),
            textangle=-25,
            xanchor="center",
            yanchor="middle"
        )

        
#%%
# GENERATE SUMMARY TEXR
#______________________________________________________________________________

# Find the most recent (youngest) ancient individual in the shared lineage
if lineage_samples_list:
    best_ancestor = max(lineage_samples_list, key=lambda x: x["year_from"])
else:
    best_ancestor = None

# Extract User Info
user1_name = (
    "Test User 1" if user1_method == "Test User"
    else sample1 if user1_method == "Database Sample"
    else f"Manual Entry ({user1_haplo})"
)

user2_name = (
    "Test User 2" if user2_method == "Test User"
    else sample2 if user2_method == "Database Sample"
    else f"Manual Entry ({user2_haplo})"
)
# Time formatting: negative = BCE, positive = CE
def format_year(year):
    return f"{abs(int(year))} BCE" if year < 0 else f"{int(year)} CE"

user1_time = "Present" if user1_method == "Test User" else format_year(user1_year)
user2_time = "Present" if user2_method == "Test User" else format_year(user2_year)

user1_country = "Unknown" if user1_method == "Test User" else (country1 if user1_method == "Database Sample" else "Unknown")
user2_country = "Unknown" if user2_method == "Test User" else (country2 if user2_method == "Database Sample" else "Unknown")

# Generate Shared Connection Text
if best_ancestor:
    recent_connection_text = f"""
    **Shared Ancient Connection**  
    {user1_name} (Country: {user1_country}, Time: {user1_time}) and {user2_name} (Country: {user2_country}, Time: {user2_time}) share a common maternal ancestor: **{best_ancestor['identifier']}**.
    This individual likely lived in {best_ancestor['region']}, {best_ancestor['country']}, during the {best_ancestor['epoch']} period, between {format_year(best_ancestor['year_from'])} and {format_year(best_ancestor['year_to'])}.

    For other ancient connections, explore the **timeline below**.
    """
else:
    recent_connection_text = f"""
    <b>No ancient ancestor found for {user1_name} and {user2_name}.</b><br>
    These individuals have no ancient ancestor in the Database. Try selecting different samples or input methods.
    """

#%%
# DISPLAY IN STREAMLIT
#______________________________________________________________________________

# Display warning or summary box
with col2:
    # Show a warning if the same individual was selected
    if identity_warning:
        st.warning(identity_warning)

    # Show summary box if valid comparison
    if not skip_comparison:
        st.markdown(
            f"""
            <div style="
                padding: 15px; 
                border: 2px solid black; 
                background-color: #edf7fe; 
                border-radius: 8px; 
                font-size: 16px;
            ">
                {recent_connection_text}
            """, 
            unsafe_allow_html=True
        )

# Timeline and Lineage Visuals
    if not skip_comparison:
        # Timeline plot section
        st.subheader("Timeline")
        st.plotly_chart(fig, use_container_width=True)

        # Lineage section for both users and shared ancestry
        st.subheader("Individual and Shared Lineages")

        # CSS for layout styling
        st.markdown(
            f"""
            <style>
                .haplo-container {{
                    display: flex;
                    justify-content: center;
                    align-items: flex-start;
                    gap: 30px;
                    margin-bottom: 20px;
                    flex-wrap: wrap;
                }}
                .haplo-left {{
                    display: flex;
                    flex-direction: column;
                    gap: 15px;
                }}
                .haplo-box {{
                    text-align: center;
                    padding: 12px;
                    font-size: 14px;
                    font-weight: bold;
                    border-radius: 10px;
                    border: 2px solid #f8f8f8;
                    background-color: #f8f8f8;
                    min-width: 180px;
                }}
                .haplo-shared {{
                    background-color: #f8f8f8;
                    border: 2px solid #f8f8f8;
                    padding: 15px;
                }}
                .lineage-box {{
                    font-size: 12px;
                    margin-top: 5px;
                    font-weight: normal;
                    color: #333;
                }}
            </style>

            <div class="haplo-container">
                <div class="haplo-left">
                    <div class="haplo-box">
                        User 1 Haplogroup: {user1_haplo}
                        <div class="lineage-box">
                            { " → ".join(lineage_to_user1) if lineage_to_user1 and lineage_to_user1 != ["Not found"] else "No lineage found" }
                        </div>
                    </div>
                    <div class="haplo-box">
                        User 2 Haplogroup: {user2_haplo}
                        <div class="lineage-box">
                            { " → ".join(lineage_to_user2) if lineage_to_user2 and lineage_to_user2 != ["Not found"] else "No lineage found" }
                        </div>
                    </div>
                </div>
                <div class="haplo-box haplo-shared">
                    Shared Haplogroup: {common_mtDNA}
                    <div class="lineage-box">
                        { " → ".join(shared_lineage) if shared_lineage else "No shared lineage" }
                    </div>
                </div>
            </div>
            """,
            unsafe_allow_html=True
        )