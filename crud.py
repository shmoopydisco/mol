from statistics import mode
from sqlalchemy.orm import Session

import models


def get_reaction(db: Session, reaction_id: int):
    return db.query(models.Reaction).filter(models.Reaction.id == reaction_id).first()


def get_reaction_by_name(db: Session, name: str):
    return db.query(models.Reaction).filter(models.Reaction.name == name).first()


def get_reactions_by_substance(db: Session, substance: str):
    return (
        db.query(models.Reaction)
        .filter(models.Reaction.substance.contains(substance))
        .all()
    )


def get_reaction_reagents_by_substance(db: Session, substance: str):
    return (
        db.query(models.Reaction.reagent)
        .filter(models.Reaction.substance.contains(substance))
        .all()
    )


def get_reactions_by_substance_and_reagent(db: Session, substance: str, reagent: str):
    return (
        db.query(models.Reaction)
        .filter(models.Reaction.substance.contains(substance))
        .filter(models.Reaction.reagent.contains(reagent))
        .all()
    )


def get_reactions(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Reaction).offset(skip).limit(limit).all()


def get_reagents(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Reaction.reagent).offset(skip).limit(limit).all()
