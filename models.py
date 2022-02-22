from sqlalchemy import Column, Integer, String

from database import Base


class Reaction(Base):
    __tablename__ = "reactions"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String)
    substance = Column(String)
    reagent = Column(String)
    environment = Column(String)
    product = Column(String)
