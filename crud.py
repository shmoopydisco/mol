from sqlalchemy.orm import Session

import database


def get_reaction(db: Session, reaction_id: int):
    return db.query(database.Reaction).filter(database.Reaction.id == reaction_id).first()


def get_reaction_by_name(db: Session, name: str):
    return db.query(database.Reaction).filter(database.Reaction.name == name).first()


def get_reactions_by_substance(db: Session, substance: str):
    return (
        db.query(database.Reaction)
            .filter(database.Reaction.substance.contains(substance))
            .all()
    )


def get_reaction_reagents_by_substance(db: Session, substance: str):
    return (
        db.query(database.Reaction.reagent)
            .filter(database.Reaction.substance.contains(substance))
            .all()
    )


def get_reactions_by_substance_and_reagent(db: Session, substance: str, reagent: str):
    return (
        db.query(database.Reaction)
            .filter(database.Reaction.substance.contains(substance))
            .filter(database.Reaction.reagent.contains(reagent))
            .all()
    )


def get_reactions(db: Session, skip: int = 0, limit: int = 100):
    return db.query(database.Reaction).offset(skip).limit(limit).all()


def get_reagents(db: Session, skip: int = 0, limit: int = 100):
    return db.query(database.Reaction.reagent).offset(skip).limit(limit).all()
