import streamlit as st
import pandas as pd
import base64
import pickle
import numpy as np
from rdkit import Chem
import pubchempy as pcp

class PubChemFingerprint:
    def __init__(self):
        self.get_pubchem_compounds = pcp.get_compounds

    def calculate(self, smiles):
        pubchem_compound = pcp.get_compounds(smiles, 'smiles')[0]
        feature = [int(bit) for bit in pubchem_compound.cactvs_fingerprint]
        return np.asarray(feature)


# Custom CSS styles with more accurate targeting
st.markdown("""
    <style>
    /* Sidebar background color */
    [data-testid="stSidebar"] {
        background-color: #969696;
    }

    /* Main app background color */
    [data-testid="stAppViewContainer"] {
        background-color: #000000;
    }

    /* Button styling in the sidebar and app */
    button[kind="primary"] {
        background-color: #006EFF;
        border: none;
        padding: 10px;
        font-size: 16px;
        border-radius: 32px;
        cursor: pointer;
    }



    /* Sidebar text color */
    [data-testid="stSidebar"] .css-hxt7ib { 
        color: black;
    }

    /* Ensuring all main content text is white */
    [data-testid="stAppViewContainer"] .css-10trblm {
        color: white;
    }
    </style>
    """, unsafe_allow_html=True)


# Molecular descriptor calculator
def desc_calc(smiles_list, ids):
    fingerprints = []
    pubchem_fp = PubChemFingerprint()
    for smiles in smiles_list:
        fp = pubchem_fp.calculate(smiles)
        fingerprints.append(fp)
       # Create a DataFrame with molecule IDs as the first column
    df = pd.DataFrame(fingerprints)
    
    # Rename columns to "PubchemFP0", "PubchemFP1", ..., "PubchemFP880"
    df.columns = [f'PubchemFP{i}' for i in range(df.shape[1])]
    
    # Insert the molecule IDs as the first column
    df.insert(0, 'ChEMBL_ID', ids)
    
    return df

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href


# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)


# Page title
st.markdown("""
# Protein Bioactivity Prediction App

This app predicts the modularity of drugs targeting the inhibition of the `Acetylcholinesterase` enzyme. It takes in the atomic details of molecules in forms of canonical smiles notations, calculates their Pubchem fingerprints (using PaDEL & Java,) eliminates the salts, organic acids, & standardizes their nitro groups. After preprocessing, it uses the high variance functional groups (subset of descriptors previously selected in training) in a multi-layer ensemble model & predicts their corresponding pIC50. Further explanation:

Firstly, a dataset of drugs with the target protein of Acetylcholinesteras was picked from Chembl database. To describe pharmacokinetic properties of the molecules, Lipinski descriptors were calculated (using rdkit.)
Next, IC50 values were converted to pIC50 scale (negative logarithmic transformation) & Mann-Whitney U test was performed for pIC50 & Lipinski discriptors: with threshold values (IC50 < 1,000 nM = Actives while IC50 › 10,000 nM = Inactives, pIC50 › 6 = Actives & pIC50 < 5 = Inactives) pIC50 of actives & inactives displayed statistically significant difference. Of the 4 Lipinski's descriptors (MW, LogP, NumHDonors & NumHAcceptors), except LogP, all other 3 showed significant difference as well.

Lastly, to describe the physico-chemical features & see which functional groups are essential to design the model (& subsequently a high potency drug), Pubchem fingerprints were calculated. After a 10- fold cross validatinon, permutative sensitivity analysis, & gridsearch hyper-parameter optimization, over 40 different models were trained & fine tuned. Model comparison showed LGBM, HistGradient, XGBR, RF & GPR had the best scores, but XGBR & Gaussian Process Regressor were picked due to their calculation time to RMSE/ R2 ratio. Ultimately a hybrid, meta-learner ensemble model blends the XGBR & GPR's predictions as the final output.

---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    # Upload functionality remains as is
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])

    # Adding the Example Input File Download Button
    st.sidebar.markdown("### Download Example Input File")
    
    # Read the local example file and provide a download button for users
    with open('example_acetylcholinesterase.txt', 'rb') as example_file:
        st.sidebar.download_button(
            label="Download Example File",
            data=example_file,
            file_name='example_acetylcholinesterase.txt',
            mime='text/plain'
        )

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    smiles_list = load_data[0].tolist()  # Assuming SMILES are in the first column
    ids = load_data[1].tolist()
    
    with st.spinner("Calculating descriptors..."):
        desc = desc_calc(smiles_list, ids)

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')
