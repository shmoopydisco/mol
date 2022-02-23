import py3Dmol
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from stmol import showmol
from streamlit_option_menu import option_menu

from apps.functional_groups import find_all_functional_groups_mode
from apps.iupac import iupac_to_struct_mode, struct_to_iupac_mode
from apps.reactions import (
    match_reaction_by_reagent_mode,
    show_all_reactions_from_struct_mode,
)


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


if __name__ == "__main__":
    st.set_page_config(
        page_title="OrgChem Helper",
        page_icon="💊",
        initial_sidebar_state="expanded",
        menu_items={
            "About": "I hope this app will help with some of the basic needs you might come across in your organic chemistry course.",
        },
    )

    with st.sidebar.title("What to do"):
        app_mode = option_menu(
            "Choose the app mode",
            [
                "Structure to IUPAC Name",
                "IUPAC Name to Structure",
                "Find All Funcional Groups",
                "Show All Reactions From Structure",
                "Match Reaction By Reagent",
            ],
            default_index=1,
        )

    match app_mode:
        case "Structure to IUPAC Name":
            struct_to_iupac_mode()
        case "IUPAC Name to Structure":
            iupac_to_struct_mode()
        case "Find All Funcional Groups":
            find_all_functional_groups_mode()
        case "Show All Reactions From Structure":
            show_all_reactions_from_struct_mode()
        case "Match Reaction By Reagent":
            match_reaction_by_reagent_mode()
