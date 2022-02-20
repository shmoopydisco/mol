import py3Dmol
import requests
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from stmol import showmol
from st_jsme import st_jsme


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
        "Ester": "[#6][CX3](=O)[OX2H0][#6]",
        "Acyl chloride": "[C,$([C]([#6])(=[O]))](=O)[Cl]",
        "Amide": "C(=O)-N",
        "Aldehyde": "[CX3H1](=O)[#6]",
        "Ketone": "[#6][CX3](=O)[#6]",
        "Primary alcohol": "[CH2][OH1]",
        "Secondary alcohol": "[CH1][OH1]",
        "Tertiary alcohol": "[OX2H][CX4;$([H0])]",
        "Primary amine": "[#6][NX3;H2;!$(NC=O)]([H])[H]",
        "Secondary amine": "[#6][NX3;H;!$(NC=O)]([#6])[H]",
        "Tertiary amine": "[#6][NX3;H0;!$(NC=O);!$(N=O)]([#6])[#6]",
        "Alkene": "[C]=[C]",
        "Alkyne": "[C]#[C]",
        "AlkylHalide": "[CX4][FX1,ClX1,BrX1,IX1]",
        "Primary alkyl halide": "[CH2][X]",
        "Secondary alkyl halide": "[CH1][X]",
        "Tertiary alkyl halide": "[C][X]",
    }

    for group, smarts in group_to_smarts_map.items():
        group_to_find = Chem.MolFromSmarts(smarts)
        res = mol.GetSubstructMatches(group_to_find)
        if len(res) != 0:
            st.write(group)
            st.write(res)
            render_2d_mol(smiles, highlight_atoms_list=res[0])


def get_iupac_name(smiles):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return f"Failed to get IUPAC name! err: {response.status_code}"


def main():
    smiles = st_jsme('500x', '350px', 'C')
    st.write(smiles)
    if smiles:
        st.header(get_iupac_name(smiles))
    identify_functional_groups(smiles)


if __name__ == "__main__":
    main()
