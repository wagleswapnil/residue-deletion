import yaml

from resdel.config.models import (
    Config,
    IOConfig,
    StateConfig,
    SystemConfig,
    TransformConfig
)


def load_config(config_file: str) -> Config:
    with open(config_file) as f:
        raw = yaml.safe_load(f)

    system = SystemConfig(**raw["system"])
    io = IOConfig(**raw["io"])
    states = {
        name: StateConfig(**state)
        for name, state in raw["transform"]["states"].items()
    }
    transform_data = raw["transform"].copy()
    transform_data["states"] = states
    transform = TransformConfig(**transform_data)
    return Config(
        system=system,
        io=io,
        transform=transform
    )
