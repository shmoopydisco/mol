from pydantic import BaseModel


class Reaction(BaseModel):
    id: int
    name: str
    substance: str
    reagent: str
    environment: str
    product: str

    class Config:
        orm_mode = True
