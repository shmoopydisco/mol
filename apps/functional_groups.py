import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from st_jsme import st_jsme


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


def render_2d_mol(smiles, highlight_atoms_list=[]):
    mol = Chem.MolFromSmiles(smiles)
    mol_image = Draw.MolToImage(mol, highlightAtoms=highlight_atoms_list)
    st.image(mol_image)


def find_all_functional_groups_mode():
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", st.session_state.smiles)

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

    return smiles
