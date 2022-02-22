import urllib.parse

import pubchempy
import py3Dmol
import requests
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from sqlalchemy.orm import Session
from st_jsme import st_jsme
from stmol import showmol

import crud, models, schemas
from database import SessionLocal, engine


def smiles_to_mol_block(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)
    return mol_block


def render_3d_mol(smiles):
    mol_block = smiles_to_mol_block(smiles)
    xyzview = py3Dmol.view()  # (width=400,height=400)
    xyzview.addModel(mol_block, "mol")
    xyzview.setStyle({"stick": {}})
    xyzview.setBackgroundColor("white")
    xyzview.zoomTo()
    showmol(xyzview, height=500, width=500)


def render_2d_mol(smiles, highlight_atoms_list=[]):
    mol = Chem.MolFromSmiles(smiles)
    mol_image = Draw.MolToImage(mol, highlightAtoms=highlight_atoms_list)
    st.image(mol_image)


def identify_functional_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)

    group_to_smarts_map = {
        "Carboxylic acid ": "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]",
        "Anhydride": "[CX3;$([H0][#6]),$([H1])](=[OX1])[#8X2][CX3;$([H0][#6]),$([H1])](=[OX1])",
        "Ester": "[#6][CX3](=O)[OX2H0][#6]",
        "Acyl chloride": "[C,$([C]([#6])(=[O]))](=O)[Cl]",
        "Amide": "C(=O)-N",
        "Aldehyde": "[CX3H1](=O)[#6]",
        "Ketone": "[#6][CX3](=O)[#6]",
        "Benzene": "c1ccccc1",
        "Hemiacetal": "[OX2H][CX4H1,!$(C(O)(O)[!#6])][OX2][#6;!$(C=[O,S,N])]",
        "Primary alcohol": "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]",
        "Secondary alcohol": "[CH1][OH1]",
        "Tertiary alcohol": "[OX2H][CX4;$([H0])]",
        # "Alcohol": "[OX2H][CX4;!$(C([OX2H])[O,S,#7,#15])]",
        "Ether": "[OD2;!$(OC~[!#1!#6])]([#6])[#6]",
        "Primary amine": "[NX3H2,NX3H2+0,NX4H3+;!$([N]~[#7,#8,#15,#16])]",
        "Secondary amine": "[NX3H1,NX3H1+0,NX4H2+;!$([N]~[#7,#8,#15,#16])]",
        "Tertiary amine": "[#6][NX3;H0;!$(NC=O);!$(N=O)]([#6])[#6]",
        "Alkene": "[C]=[C]",
        "Alkyne": "[C]#[C]",
        "AlkylHalide": "[CX4][FX1,ClX1,BrX1,IX1]",
        "Primary alkyl halide": "[CH2][F,Cl,Br,I]",
        "Secondary alkyl halide": "[CH1][F,Cl,Br,I]",
        "Tertiary alkyl halide": "[C][F,Cl,Br,I]",
        "Alkane": "[CX4;$([H3][#6]),$([H2]([#6])[#6]),$([H1]([#6])([#6])[#6]),$([#6]([#6])([#6])([#6])[#6])]",
    }

    functional_groups = []
    for group, smarts in group_to_smarts_map.items():
        group_to_find = Chem.MolFromSmarts(smarts)
        res = mol.GetSubstructMatches(group_to_find)
        if len(res) != 0:
            functional_groups.append(group)
            st.write(group)
            st.write(res)
            render_2d_mol(smiles, highlight_atoms_list=res[0])

    return functional_groups


def fetch_and_display_iupac_name(smiles):
    if not smiles:
        return

    encoded_smiles = urllib.parse.quote_plus(smiles)
    url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_smiles}/iupac_name"
    response = requests.get(url)
    if response.status_code == 200:
        st.write(response.text)
    else:
        try:
            compounds = pubchempy.get_compounds(smiles, namespace="smiles")
            name = compounds[0].iupac_name
            if not name:
                raise pubchempy.NotFoundError
            else:
                st.warning(
                    "Trying alternative IUPAC source, results may be less accurate"
                )
                st.write(name)

        except (pubchempy.BadRequestError, pubchempy.NotFoundError):
            st.error("Failed getting IUPAC name!")


def present_possible_reactions(functional_groups):
    if not functional_groups:
        return

    db = SessionLocal()
    possible_reagents = set()
    try:
        for fg in functional_groups:
            raw_reagents = crud.get_reaction_reagents_by_substance(db, fg)
            reagents = [item[0] for item in raw_reagents]
            possible_reagents.update(reagents)
    finally:
        db.close()

    st.write(possible_reagents)

    with st.form("reactions_form"):
        reagents_selectbox = st.selectbox("Choose a reagent:", possible_reagents)
        submitted = st.form_submit_button("Submit")

        if submitted:
            possible_reactions = set()
            try:
                for fg in functional_groups:
                    raw_reactions = crud.get_reactions_by_substance_and_reagent(
                        db, substance=fg, reagent=reagents_selectbox
                    )
                    # reactions = [item[0] for item in raw_reactions]
                    possible_reactions.update(raw_reactions)

            finally:
                db.close()

            st.subheader("Possible reactions: ")

            col1, col2, col3 = st.columns(3)

            with col1:
                st.subheader("Name")
                for reaction in possible_reactions:
                    st.write(reaction.name)

            with col2:
                st.subheader("Environment")
                for reaction in possible_reactions:
                    st.write(reaction.environment)

            with col3:
                st.subheader("Products")
                for reaction in possible_reactions:
                    st.write(reaction.product)
                    st.empty()


def main():
    st.set_page_config(
        page_title="OrgChem Helper",
        page_icon="💊", #"🧪",
        initial_sidebar_state="expanded",
        menu_items={
            "Get Help": "https://www.extremelycoolapp.com/help",
            "Report a bug": "https://www.extremelycoolapp.com/bug",
            "About": "# This is a header. This is an *extremely* cool app!",
        },
    )
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", "CCC")
    st.subheader("SMILES (for debugging)")
    st.write(smiles)

    st.subheader("IUPAC Name")
    fetch_and_display_iupac_name(smiles)

    st.subheader("Functional Groups")
    functional_groups = []
    try:
        functional_groups = identify_functional_groups(smiles)
    except TypeError:
        pass
    st.write(functional_groups)

    st.subheader("Reactions")
    present_possible_reactions(functional_groups)


if __name__ == "__main__":
    main()
