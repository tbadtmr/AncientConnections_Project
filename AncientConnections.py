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
#FUNCTIONS

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
        # Validate haplogroup entry
            if haplo in phylo_tree:
                user1_haplo = haplo
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
        user2_haplo = None
        user2_haplo_input = st.text_input("Enter mtDNA haplogroup", key="user2_haplo")
        user2_year = st.number_input("Enter estimated year (negative for BCE)", value=1000, key="user2_year")

        if user2_haplo_input:
            haplo = user2_haplo_input.strip().upper()
        # Validate haplogroup entry
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

#%%
# GET MOST RECENT HAPLOGROUP AND LINEAGE
#______________________________________________________________________________





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

# DEFINE SCALES
#------------------------------------------------------------------------------
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

# **Create Timeline Data**
labels = ["User 1", "User 2"]
years = [user1_year, user2_year]
colors = ["orange", "purple"]
hover_texts = [f"User 1: {user1_haplo}", f"User 2: {user2_haplo}"]

# **Define Starting Y Positions**
y_base = 2.0  # ✅ Start high for oldest ancestors
y_positions = [1.0, 1.4]  # Users at bottom

# **Add Additional Ancestors (from Oldest to Newest)**
if additional_ancestors_list:
    for i, ancestor in enumerate(additional_ancestors_list):
        labels.insert(i, ancestor["identifier"])
        years.insert(i, ancestor["year_from"])
        y_positions.insert(i, y_base - (i * 0.2))  # ✅ Stepwise descending Y
        colors.insert(i, "lightblue")

        # **Enhanced Hover Text with More Details**
        hover_texts.insert(i, (
            f"<b>Shared Ancestor</b>: {ancestor['identifier']}<br>"
            f"<b>Haplogroup</b>: {ancestor['mt_hg']}<br>"
            f"<b>Country</b>: {ancestor['country']}<br>"
            f"<b>Region</b>: {ancestor['region']}<br>"
            f"<b>Lifespan</b>: {ancestor['year_from']} to {ancestor['year_to']}"
        ))

# **Add Most Ancient Shared Ancestor**
if show_ancient:
    labels.insert(len(additional_ancestors_list), most_ancient_identifier)
    years.insert(len(additional_ancestors_list), most_ancient_year)
    y_positions.insert(len(additional_ancestors_list), y_base - (len(additional_ancestors_list) * 0.2))
    colors.insert(len(additional_ancestors_list), "darkblue")

    # **Enhanced Hover Text**
    hover_texts.insert(len(additional_ancestors_list), (
        f"<b>Most Ancient Shared Ancestor</b>:{most_ancient_identifier}<br>"
        f"<b>Haplogroup</b>: {most_ancient_haplogroup}<br>"
        f"<b>Country</b>: {most_ancient_country}<br>"
        f"<b>Region</b>: {most_ancient_region}<br>"
        f"<b>Lifespan</b>: {most_ancient_year} to {most_ancient_year2}"
    ))

# **Add Most Recent Shared Ancestor**
if show_recent:
    labels.insert(len(additional_ancestors_list) + 1, most_recent_identifier)
    years.insert(len(additional_ancestors_list) + 1, most_recent_year)
    y_positions.insert(len(additional_ancestors_list) + 1, y_positions[-1] + 0.2)  # ✅ Just above users
    colors.insert(len(additional_ancestors_list) + 1, "darkblue")

    # **Enhanced Hover Text**
    hover_texts.insert(len(additional_ancestors_list) + 1, (
        f"<b>Most Recent Shared Ancestor</b>: {most_recent_identifier}<br>"
        f"<b>Haplogroup</b>: {most_recent_haplogroup}<br>"
        f"<b>Country</b>: {most_recent_country}<br>"
        f"<b>Region</b>: {most_recent_region}<br>"
        f"<b>Lifespan</b>: {most_recent_year} to {most_recent_year2}"
    ))

# **Add Shared Haplogroup**
if show_common_mtDNA:
    shared_haplogroup_year = min(user1_year, user2_year) - 500  # ✅ Place it further back
    shared_haplogroup_pos = (y_positions[0] + y_positions[-1]) / 2  # ✅ Place it between users & ancestors
    labels.insert(len(additional_ancestors_list), common_mtDNA)
    years.insert(len(additional_ancestors_list), shared_haplogroup_year)
    y_positions.insert(len(additional_ancestors_list), shared_haplogroup_pos)
    colors.insert(len(additional_ancestors_list), "black")  # ✅ Neutral color for shared haplogroup

    # **Enhanced Hover Text**
    hover_texts.insert(len(additional_ancestors_list), (
        f"<b>Shared Haplogroup</b>: {common_mtDNA}<br>"
        "This is the most recent common haplogroup shared by both users."
    ))


# **Create Plotly Figure**
fig = go.Figure()

# **Add Scatter Points**
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
#------------------------------------------------------------------------------
# **Bezier Curve Function for Smooth Connections**
def bezier_curve(x_start, x_end, y_start, y_end, curve_direction="up", num_points=100, stop_offset=5):
    """ Generate curved Bezier path for smooth connections. """
    t = np.linspace(0, 1, num_points)
    # **Dynamic Control Point:**
    control_y = (y_start + y_end) / 2  # Default middle control point
    control_x = (x_start + x_end) / 2  # Midpoint in X
    # **Push Control Point Higher or Lower to Strengthen the Curve**
    vertical_distance = abs(y_end - y_start)
    curve_strength = vertical_distance * 0.5  # Increase strength based on spacing
    if curve_direction == "up":
        control_y += curve_strength  # Push control point UP
    else:
        control_y -= curve_strength  # Push control point DOWN
    # **Bezier Curve Calculation**
    bezier_x = (1 - t) ** 2 * x_start + 2 * (1 - t) * t * control_x + t ** 2 * x_end
    bezier_y = (1 - t) ** 2 * y_start + 2 * (1 - t) * t * control_y + t ** 2 * y_end
    # **Apply Offset to Avoid Sharp Stops**
    bezier_x[-1] = bezier_x[-1] - stop_offset  

    return bezier_x, bezier_y

# Curved Connections
y_position_map = {label: y for label, y in zip(labels, y_positions)}

# Connect most recent
if show_recent and show_ancient:
    # Connect users to Most Recent
    if most_recent_identifier in y_position_map:
        x_vals, y_vals = bezier_curve(user1_year, most_recent_year, 
                                      y_position_map["User 1"], y_position_map[most_recent_identifier], 
                                      (y_position_map["User 1"] + y_position_map[most_recent_identifier]) / 2)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="blue", width=2)))

        x_vals, y_vals = bezier_curve(user2_year, most_recent_year, 
                                      y_position_map["User 2"], y_position_map[most_recent_identifier], 
                                      (y_position_map["User 2"] + y_position_map[most_recent_identifier]) / 2)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="blue", width=2)))

    # Connect Most Recent to Most Ancient
    if most_ancient_identifier in y_position_map:
        x_vals, y_vals = bezier_curve(most_recent_year, most_ancient_year, 
                                      y_position_map[most_recent_identifier], y_position_map[most_ancient_identifier], 
                                      (y_position_map[most_recent_identifier] + y_position_map[most_ancient_identifier]) / 2)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="blue", width=2)))

elif show_ancient:
    # If no Most Recent, connect Users directly to Most Ancient
    if most_ancient_identifier in y_position_map:
        x_vals, y_vals = bezier_curve(user1_year, most_ancient_year, 
                                      y_position_map["User 1"], y_position_map[most_ancient_identifier], 
                                      (y_position_map["User 1"] + y_position_map[most_ancient_identifier]) / 2)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="blue", width=2)))

        x_vals, y_vals = bezier_curve(user2_year, most_ancient_year, 
                                      y_position_map["User 2"], y_position_map[most_ancient_identifier], 
                                      (y_position_map["User 2"] + y_position_map[most_ancient_identifier]) / 2)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="blue", width=2)))

# Connect Users & Additional Ancestors to Shared Haplogroup
if show_common_mtDNA:
    # Connect Users
    x_vals, y_vals = bezier_curve(user1_year, shared_haplogroup_year, y_position_map["User 1"], shared_haplogroup_pos, "up")
    fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="black", width=2, dash="dash")))

    x_vals, y_vals = bezier_curve(user2_year, shared_haplogroup_year, y_position_map["User 2"], shared_haplogroup_pos, "up")
    fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="black", width=2, dash="dash")))

    # Connect Additional Ancestors to Shared Haplogroup
    if additional_ancestors_list:
        x_vals, y_vals = bezier_curve(
            additional_ancestors_list[-1]["year_from"], shared_haplogroup_year,
            y_position_map[additional_ancestors_list[-1]["identifier"]], shared_haplogroup_pos,
            (y_position_map[additional_ancestors_list[-1]["identifier"]] + shared_haplogroup_pos) / 2
        )
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="black", width=2, dash="dash")))


# Connections to Additional ANcestors
if additional_ancestors_list:
    for i in range(len(additional_ancestors_list) - 1):
        older_ancestor = additional_ancestors_list[i]
        younger_ancestor = additional_ancestors_list[i + 1]

        if older_ancestor["identifier"] in y_position_map and younger_ancestor["identifier"] in y_position_map:
            x_vals, y_vals = bezier_curve(
                older_ancestor["year_from"], younger_ancestor["year_from"],
                y_position_map[older_ancestor["identifier"]], y_position_map[younger_ancestor["identifier"]],
                (y_position_map[older_ancestor["identifier"]] + y_position_map[younger_ancestor["identifier"]]) / 2
            )
            fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="lightblue", width=2)))

    # **Connect last additional ancestor to Most Ancient OR Shared Haplogroup**
    if show_ancient and additional_ancestors_list[-1]["identifier"] in y_position_map:
        x_vals, y_vals = bezier_curve(
            additional_ancestors_list[-1]["year_from"], most_ancient_year,
            y_position_map[additional_ancestors_list[-1]["identifier"]], y_position_map[most_ancient_identifier],
            (y_position_map[additional_ancestors_list[-1]["identifier"]] + y_position_map[most_ancient_identifier]) / 2
        )
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="lightblue", width=2)))

    elif show_common_mtDNA and additional_ancestors_list[-1]["identifier"] in y_position_map:
        # **Connect last additional ancestor to Shared Haplogroup instead of Users**
        x_vals, y_vals = bezier_curve(
            additional_ancestors_list[-1]["year_from"], shared_haplogroup_year,
            y_position_map[additional_ancestors_list[-1]["identifier"]], shared_haplogroup_pos,
            (y_position_map[additional_ancestors_list[-1]["identifier"]] + shared_haplogroup_pos) / 2
        )
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", line=dict(color="black", width=2, dash="dash")))

# **Step 2: Dynamically Filter Epochs in Range**
filtered_epochs = []
epoch_positions = []
epoch_tick_labels = []
min_year, max_year = min(years), max(years)

for epoch, (start, end) in epoch_definitions.items():
    if start <= max_year and end >= min_year:  # Only include epochs within range
        filtered_epochs.append(epoch)
        epoch_positions.append((start + end) / 2)  # Center label in its period
        epoch_tick_labels.append(epoch)  # Add text label

# **Step 3: Define Time Scale (BCE/CE Formatting)**
tick_positions = np.linspace(min_year, max_year, num=10).astype(int)
tick_labels = [f"{abs(t)} BCE" if t < 0 else f"{t} CE" for t in tick_positions]

# **Step 4: Update Plot Layout**
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
    yaxis=dict(visible=False),  # Hide y-axis
    showlegend=False,
    margin=dict(l=60, r=60, t=60, b=120),  
)

# **Step 5: Add Colored Epoch Bars & Labels**
for epoch, (start, end) in epoch_definitions.items():
    if start <= max_year and end >= min_year:
        x0 = max(start, min_year)
        x1 = min(end, max_year)

        #  Add the epoch background color bar
        fig.add_shape(
            type="rect",
            x0=x0, x1=x1,
            y0=-0.25, y1=-0.05,  
            fillcolor=epoch_colors[epoch],
            opacity=0.4, layer="below", line=dict(width=0)
        )

        text_x = (x0 + x1) / 2  # Center of epoch range
        fig.add_annotation(
            x=text_x,
            y=-0.15,  # Positioned just below the bar
            text=epoch,
            showarrow=False,
            font=dict(size=11, color="black"),  # Adjusted for readability
            textangle=-25,  # Slightly tilted for clarity
            xanchor="center",
            yanchor="middle"
        )
        
# GENERATE TEXT ABOUT RECENT CONNECTION

# **Determine the Most Recent Connection**
if show_recent:
    best_ancestor = {
        "identifier": most_recent_identifier,
        "region": most_recent_region,
        "country": most_recent_country,
        "epoch": most_recent_epoch,
        "year_from": most_recent_year,
        "year_to": most_recent_year2
    }
elif show_ancient:
    best_ancestor = {
        "identifier": most_ancient_identifier,
        "region": most_ancient_region,
        "country": most_ancient_country,
        "epoch": most_ancient_epoch,
        "year_from": most_ancient_year,
        "year_to": most_ancient_year2
    }
elif additional_ancestors_list:
    # **Find the youngest additional ancestor**
    best_ancestor = additional_ancestors_list[-1]  # Last one is the most recent
else:
    best_ancestor = None  # No ancestor found

# **Extract User Information**
user1_name = "Test User 1" if user1_method == "Test User" else (sample1 if user1_method == "Database Sample" else "Manual Entry")
user2_name = "Test User 2" if user2_method == "Test User" else (sample2 if user2_method == "Database Sample" else "Manual Entry")

user1_time = "Present" if user1_method == "Test User" else f"{user1_year} BCE"
user2_time = "Present" if user2_method == "Test User" else f"{user2_year} BCE"

user1_country = "Unknown" if user1_method == "Test User" else (country1 if user1_method == "Database Sample" else "Unknown")
user2_country = "Unknown" if user2_method == "Test User" else (country2 if user2_method == "Database Sample" else "Unknown")

# **Generate Shared Connection Text**
if best_ancestor:
    recent_connection_text = f"""
    **Shared Ancient Connection**  
    {user1_name} (Country: {user1_country}, Time: {user1_time}) and {user2_name} ({user2_country}, {user2_time}) share a common maternal ancestor: **{best_ancestor['identifier']}**.
    This individual likely lived in {best_ancestor['region']}, {best_ancestor['country']}, during the {best_ancestor['epoch']} period, between {best_ancestor['year_from']} and {best_ancestor['year_to']} .

    For other ancient connections, explore the **timeline below**.
    """
else:
    if show_ancient == False and additional_ancestors.empty:
        recent_connection_text = f"""
        <b> No ancient ancestor found for {user1_name} and {user2_name}.</b>  
        This lineage has no direct ancestors in the database, please explore a different connection."""

# **Display in Streamlit**
with col2:
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
    

#%%
# DISPLAY EVERYTHING IN COLUMN 2
#______________________________________________________________________________
with col2:
    
    

    # **Show Plot**
    st.subheader("Timeline")
    st.plotly_chart(fig, use_container_width=True)
    
    # **Custom CSS for New Layout**
    st.subheader("Individual and Shared Lineages")
    st.markdown(
        f"""
        <style>
            .haplo-container {{
                display: flex;
                justify-content: center;
                align-items: center;
                gap: 30px;
                margin-bottom: 20px;
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

