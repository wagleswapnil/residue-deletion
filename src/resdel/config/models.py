from dataclasses import dataclass
from typing import Optional

@dataclass
class SystemConfig:
    name: str
    sequence: str
    residue_to_delete: str

@dataclass
class IOConfig:
    working_dir: Optional[str]
    input_mdp_dir: Optional[str]
    output_dir: Optional[str]
    output_prefix: Optional[str]

@dataclass
class StateConfig:
    topology: str
    structrue: str
    moleucle_name: str

@dataclass
class TransformConfig:
    molecule_name: Optional[str]
    number_of_lambdas: int
    lambda_vector: Optional[list]
    vdw_lambda_vector: Optional[list]
    coul_lambda_vector: Optional[list]
    states: dict

@dataclass
class Config:
    system: SystemConfig
    io: IOConfig
    transform: TransformConfig

