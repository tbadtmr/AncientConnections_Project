#%% 
# DOCUMENTATION SECTION - AADR Dataset Parser
#______________________________________________________________________________
#!/usr/bin/python3
"""
Description:
    This script processes an ancient DNA dataset in excel format (AADR) and extracts relevant information 
    for downstream analysis (Ancient Connections App). It dynamically identifies and selects key columns
    such as haplogroups, sex, location etc. It performs data cleaning, standardizes haplogroup naming, assigns
    moder countries and continents, extracts the earliest year from date strings and categorizes samples into 
    broad historical epochs.

User-defined functions:
    find_closest_column(): Matches approximate column names in the dataset
    is_invalid(): Identifies invalid or missing haplogroup entries
    clean_haplogroup(): Extracts base-level haplogroup notation
    assign_continent(): Maps country names to continents
    extract_earliest_year(): Parses date strings to find earliest associated year
    assign_epoch(): Assigns historical periods based on extracted year

Procedure:
    1. Load Excel file
    2. Identify key columns, select and rename columns
    3. Clean haplogroup data
    4. Standardize geographic information, map country and continent
    5. Parse time data and classify into epochs
    6. Output cleaned Excel file

Usage:
    This script is designed to be executed directly in a Python environment.
    File paths are defined in the script and can be adapted as needed.

Version: 1.00
Date: 2025-03-13
Name: Tabea Dittmar
"""
#%%
# IMPORT MODULES AND FILES
#______________________________________________________________________________
import pandas as pd
import re

# Define the input and output file paths
input_file = r"\\Data\AADR_Annotation.xlsx"
output_file = r"\\Data\AADR_cleaned.xlsx"

# Load the Excel file
xls = pd.ExcelFile(input_file)

# Select the first sheet
df = xls.parse(xls.sheet_names[0])

#%%
# Find and select columns
#______________________________________________________________________________

# Function to find the closest column name
def find_closest_column(target_name, df_columns):
    for col in df_columns:
        if target_name.lower().strip() in col.lower().strip():
            return col
    return None

# Find the correct column names dynamically
correct_mtdna_col = find_closest_column("mtDNA haplogroup if >2x or published", df.columns)
correct_yhaplo_col = find_closest_column("Y haplogroup (manual curation in ISOGG format)", df.columns)
correct_political_entity_col = find_closest_column("Political Entity", df.columns)
correct_sex_col = find_closest_column("Molecular Sex", df.columns)
correct_identifier_col = find_closest_column("Master ID", df.columns)

if not correct_mtdna_col or not correct_yhaplo_col or not correct_political_entity_col or not correct_sex_col or not correct_identifier_col:
    raise ValueError("Error: Could not find required columns in the dataset.")

# Define the required columns
selected_columns = [
    "#",
    "Genetic ID",
    correct_identifier_col,
    "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE",
    "Group ID",
    "Locality",
    correct_political_entity_col,  
    correct_sex_col,  
    correct_yhaplo_col,  
    correct_mtdna_col
]

# Filter the dataframe with selected columns
df_selected = df[selected_columns].copy()

# Rename columns
df_selected.rename(columns={
    correct_political_entity_col: "Political Entity",
    correct_sex_col: "Sex",
    correct_identifier_col: "Identifier",
    correct_yhaplo_col: "ychr_hg",
    correct_mtdna_col: "mt_hg",
    "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE": "year_from"
}, inplace=True)


#%%

# HAPLOGROUP COLUMN
#-------------------

# List of invalid values
invalid_values = ["..", "n/a", "n/a (>2x)", "n/a (<2x)", "n/a (sex unknown)", "n/a (female)"]
# Function to check if a value is invalid
def is_invalid(value):
    return pd.isna(value) or any(invalid in str(value).lower() for invalid in invalid_values)

# Remove rows where both Y haplogroup and mtDNA haplogroup contain invalid values
df_selected = df_selected[
    ~((df_selected["ychr_hg"].apply(is_invalid)) &
      (df_selected["mt_hg"].apply(is_invalid)))
].copy()

# Function to clean haplogroup data by extracting only the base haplogroup
def clean_haplogroup(haplo):
    if pd.isna(haplo) or is_invalid(haplo):
        return "N/A"
    return re.match(r"^[A-Z][\d\w]*", haplo).group(0) if re.match(r"^[A-Z][\d\w]*", haplo) else haplo

# Apply cleaning function to haplogroup columns
df_selected.loc[:, "ychr_hg"] = df_selected["ychr_hg"].apply(clean_haplogroup)
df_selected.loc[:, "mt_hg"] = df_selected["mt_hg"].apply(clean_haplogroup)

# COUNTRX COLUMN
#----------------

# Mapping of Political Entity to Country
political_to_country = {
    "Crimea": "Ukraine / Russia (disputed)",
    "Canary Islands": "Spain",
    "Channel Islands": "United Kingdom",
    "Curacao": "Netherlands",
    "Czechia": "Czech Republic",
    "DR Congo": "Democratic Republic of Congo",
    "Faroes": "Denmark",
    "French Polynesia": "France",
    "Gibraltar": "United Kingdom",
    "Greenland": "Denmark",
    "Guadeloupe": "France",
    "Guam": "USA",
    "Isle of Man": "United Kingdom",
    "Micronesia": "Federated States of Micronesia",
    "Puerto Rico": "USA"
}

# Assign country based on Political Entity
df_selected["Country"] = df_selected["Political Entity"].replace(political_to_country).fillna(df_selected["Political Entity"])

# CONTINENT COLUMN
#-----------------

# Define continent mapping
continent_map = {
    "Europe": ["Albania", "Austria", "Belgium", "Bosnia-Herzegovina", "Bulgaria", "Croatia", "Czech Republic", "Cyprus", "Denmark",
               "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia",
               "Lithuania", "Luxembourg", "Montenegro", "Moldova", "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal",
               "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", "Ukraine", "United Kingdom"],
    "Asia": ["Afghanistan", "Armenia", "Azerbaijan", "China", "India", "Indonesia", "Iran", "Iraq", "Israel", "Japan",
             "Jordan", "Kazakhstan", "Kyrgyzstan", "Lebanon", "Laos", "Malaysia", "Mongolia", "Nepal", "Pakistan", "Philippines",
             "Syria", "Taiwan", "Tajikistan", "Thailand", "Turkey", "Turkmenistan", "Uzbekistan", "Vietnam", "Russia"],
    "Africa": ["Botswana", "Cameroon", "Democratic Republic of Congo", "Egypt", "Ethiopia", "Kenya", "Malawi", "Morocco",
               "South Africa", "Sudan", "Tanzania", "Zambia"],
    "North America": ["Bahamas", "Canada", "Cuba", "Dominican Republic", "Haiti", "St. Lucia", "USA", "Mexico", "Greenland", "Guam", "Panama", "Puerto Rico"],
    "South America": ["Argentina", "Belize", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "Guyana", "Paraguay",
                      "Peru", "Uruguay", "Venezuela"],
    "Oceania": ["Australia", "Fiji", "French Polynesia", "Micronesia", "New Zealand", "Papua New Guinea", "Solomon Islands",
                "Tonga", "Vanuatu"]
}


def assign_continent(country):
    for continent, countries in continent_map.items():
        if country in countries:
            return continent
    return "Unknown"

df_selected["Continent"] = df_selected["Country"].apply(assign_continent)

# YEAR COLUMN
#------------

# Extract earliest year from date column
def extract_earliest_year(date_str):
    if pd.isna(date_str):
        return None
    if "present" in date_str.lower():
        return 2024

    # Normalize dashes
    date_str = date_str.replace("–", "-").replace("−", "-")
    date_str = re.sub(r"\s*-\s*", "-", date_str.strip())

    # 1. Try to match full date+era range (e.g. 149 calBCE - 65 calCE)
    match = re.search(r"(\d{1,6})\s*(calBCE|BCE|calCE|CE)\s*-\s*(\d{1,6})\s*(calBCE|BCE|calCE|CE)", date_str)
    if match:
        y1, e1, y2, e2 = int(match.group(1)), match.group(2), int(match.group(3)), match.group(4)
    else:
        # 2. Try to match two numbers and a shared era after (e.g. 1090-900 BCE)
        match = re.search(r"(\d{1,6})-(\d{1,6})\s*(calBCE|BCE|calCE|CE)", date_str)
        if match:
            y1, y2, e1 = int(match.group(1)), int(match.group(2)), match.group(3)
            e2 = e1  # Same era assumed for both
        else:
            return None  # No valid match

    # Convert to actual year
    def convert(y, era):
        return -y if "BCE" in era else y
    year1 = convert(y1, e1)
    year2 = convert(y2, e2)

    return min(year1, year2)

df_selected["Year"] = df_selected["year_from"].apply(extract_earliest_year)

# EPOCH COLUMN
#--------------

# Define epoch classification
def assign_epoch(year):
    if year is None:
        return "Unknown"
    if year <= -10000:
        return "Paleolithic"
    elif year <= -8000:
        return "Mesolithic"
    elif year <= -3000:
        return "Neolithic"
    elif year <= -1200:
        return "Bronze Age"
    elif year <= 500:
        return "Iron Age"
    elif year <= 1500:
        return "Middle Ages"
    elif year <= 2025:
        return "Modern Age"
    elif year == "present":
        return "Modern Age"
    return "Unknown"

df_selected["Epoch"] = df_selected["Year"].apply(assign_epoch)

#%%
# WRITE OUTPUT FILE
#______________________________________________________________________________

# Save the filtered data to a new Excel file
df_selected.to_excel(output_file, index=False)

print(f"Filtered data with Political Entity, matched Country, Continent, and Epoch saved to {output_file}")
