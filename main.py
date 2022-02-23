import urllib.parse

import pandas as pd
import pubchempy
import py3Dmol
import requests
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from sqlalchemy.orm import Session
from st_jsme import st_jsme
from stmol import showmol

import crud
import models
from database import SessionLocal


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


def get_functional_groups(smiles):
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
        # "AlkylHalide": "[CX4][FX1,ClX1,BrX1,IX1]",
        "Primary alkyl halide": "[CH2][F,Cl,Br,I]",
        "Secondary alkyl halide": "[CH1][F,Cl,Br,I]",
        "Tertiary alkyl halide": "[X4&H0][F,Cl,Br,I]",
        # "Alkane": "[CX4;$([H3][#6]),$([H2]([#6])[#6]),$([H1]([#6])([#6])[#6]),$([#6]([#6])([#6])([#6])[#6])]",
        # "AlkanePrimary": "[CX4H3][#6][X]",
        # "AlkaneSecondary": "[CX4H2]([#6])[#6][X]",
        # "AlkaneTertiary": "[CX4H1]([#6])([#6])[#6][X]",
    }

    functional_groups = {}
    for group, smarts in group_to_smarts_map.items():
        group_to_find = Chem.MolFromSmarts(smarts)
        res = mol.GetSubstructMatches(group_to_find)
        if len(res) != 0:
            functional_groups[group] = res

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

    with st.form("reactions_form"):
        reagents_selectbox = st.selectbox("Choose a reagent:", possible_reagents)
        show_all_checkbox = st.checkbox("Show all relevant reactions")
        submitted = st.form_submit_button("Submit")

        if submitted:
            st.subheader("Possible reactions: ")
            possible_reactions = []
            try:
                for fg in functional_groups:
                    df = pd.read_sql(
                        db.query(models.Reaction)
                        .filter(models.Reaction.substance.contains(fg))
                        .filter(models.Reaction.reagent.contains(reagents_selectbox))
                        .statement,
                        db.bind,
                    )
                    possible_reactions.append(df)

            finally:
                db.close()

            if len(possible_reactions) > 1 and show_all_checkbox:
                hide_table_row_index = """
                <style>
                tbody th {display:none}
                .blank {display:none}
                </style>
                """
                st.markdown(hide_table_row_index, unsafe_allow_html=True)
                st.table(pd.concat(possible_reactions).drop(columns="id"))
            else:
                pass


def show_formatted_table(possible_reactions):
    hide_table_row_index = """
    <style>
    tbody th {display:none}
    .blank {display:none}
    </style>
    """
    st.markdown(hide_table_row_index, unsafe_allow_html=True)
    st.table(pd.concat(possible_reactions).drop(columns="id"))


def struct_to_iupac_mode():
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", "CCC")

    st.subheader("SMILES (for debugging)")
    st.write(smiles)

    st.subheader("IUPAC Name")
    fetch_and_display_iupac_name(smiles)


def struct_show_all_reactions_mode():
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", "CCC")

    st.subheader("SMILES (for debugging)")
    st.write(smiles)

    functional_groups = dict()
    try:
        functional_groups = get_functional_groups(smiles)
    except TypeError:
        pass

    st.subheader("Reactions")
    if not functional_groups:
        return

    db = SessionLocal()
    possible_reactions = []
    try:
        for fg in functional_groups:
            df = pd.read_sql(
                db.query(models.Reaction)
                .filter(models.Reaction.substance.contains(fg))
                .statement,
                db.bind,
            )
            possible_reactions.append(df)
    finally:
        db.close()

    show_formatted_table(possible_reactions)


def find_all_functional_groups_mode():
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", "CCC")

    st.subheader("SMILES (for debugging)")
    st.write(smiles)

    st.subheader("Functional Groups")
    functional_groups = dict()
    try:
        functional_groups = get_functional_groups(smiles)
    except TypeError:
        pass

    for group in functional_groups.items():
        st.write(group[0])
        render_2d_mol(smiles, highlight_atoms_list=group[1][0])


if __name__ == "__main__":
    st.set_page_config(
        page_title="OrgChem Helper",
        page_icon="💊",
        initial_sidebar_state="expanded",
        menu_items={
            "About": "I hope this app will help with some of the basic needs you might come across in your organic chemistry course.",
        },
    )

    st.sidebar.title("What to do")
    app_mode = st.sidebar.selectbox(
        "Choose the app mode",
        [
            "Structure to IUPAC Name",
            "Show All Reactions From Structure",
            "Find All Funcional Groups",
        ],
    )

    match app_mode:
        case "Show All Reactions From Structure":
            struct_show_all_reactions_mode()
        case "Structure to IUPAC Name":
            struct_to_iupac_mode()
        case "Find All Funcional Groups":
            find_all_functional_groups_mode()
