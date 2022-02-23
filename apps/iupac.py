import urllib.parse

import pubchempy
import requests
import streamlit as st
from st_jsme import st_jsme


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


def iupac_to_struct_mode():
    st.title("Organic Chemistry Helper")

    with st.form(key="iupac_to_struct_form"):
        iupac = st.text_input("Please enter IUPAC name to search")
        submitted = st.form_submit_button("Submit")

        if submitted:
            # url = f"https://cactus.nci.nih.gov/chemical/structure/{iupac}/image?width=1000&height=1000"
            url = f"https://opsin.ch.cam.ac.uk/opsin/{iupac}.png"
            response = requests.get(url)
            if response.status_code == 200:
                st.image(response.content)
            else:
                st.error("Failed getting a structure!")


def struct_to_iupac_mode():
    st.title("Organic Chemistry Helper")

    smiles = st_jsme("500x", "350px", "CCC")

    st.subheader("SMILES (for debugging)")
    st.write(smiles)

    st.subheader("IUPAC Name")
    fetch_and_display_iupac_name(smiles)
