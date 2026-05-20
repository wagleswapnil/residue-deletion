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
    log_file: Optional[str]
    error_file: Optional[str]

@dataclass
class PrepareConfig:
    GMX_executable: Optional[str]
    forcefield: Optional[str]
    water_model: Optional[str]

@dataclass
class StateConfig:
    topology: Optional[str]
    structure: Optional[str]
    molecule_name: Optional[str]

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
    prepare: PrepareConfig
    transform: TransformConfig

