import crud
import models
import pandas as pd
import streamlit as st
from awesome_table import AwesomeTable, Column
from database import SessionLocal
from st_jsme import st_jsme

from .functional_groups import get_functional_groups


def match_reaction_by_reagent_mode():
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
        st.error("No functional groups to react with!")
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

            show_formatted_table(possible_reactions)


def present_possible_reactions(functional_groups):
    if not functional_groups:
        st.error("No functional groups to react with!")
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
    if len(possible_reactions) < 1:
        return

    df = pd.concat(possible_reactions)
    AwesomeTable(
        df.drop(columns="id"),
        columns=[
            Column(name="name", label="Name"),
            Column(name="substance", label="Main Substance"),
            Column(name="reagent", label="Reagent"),
            Column(name="environment", label="Environment"),
            Column(name="product", label="Products"),
        ],
        show_search=True,
    )


def show_all_reactions_from_struct_mode():
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
        st.error("No functional groups to react with!")
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
