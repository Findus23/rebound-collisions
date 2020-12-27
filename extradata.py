import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

from rebound import Particle


class ExtraData:

    def __init__(self):
        self.tree = CollisionTree()
        self.pdata: Dict[int, ParticleData] = {}
        self.meta = Meta()
        self.energy = EnergyConsercation()

    def save(self, filename: Path):
        pdata = {}
        for k, v in self.pdata.items():
            pdata[k] = v.__dict__

        with filename.open("w") as f:
            json.dump({
                "meta": self.meta.save(),
                "pdata": pdata,
                "tree": self.tree.save(),
                "dEs": self.energy._dEs
            }, f, indent=2)

    @classmethod
    def load(cls, filename: Path):
        with filename.open() as f:
            data = json.load(f)
        self = cls()
        self.meta = Meta(**data["meta"])

        self.tree.load(data["tree"])
        # self.tree._dEs = data["dEs"]

        for k, v in data["pdata"].items():
            self.pdata[int(k)] = ParticleData(**v)

        return self


@dataclass
class ParticleData:
    water_mass_fraction: float


@dataclass
class Meta:
    tmax: float = None
    savesteps: int = None
    max_n: int = None
    walltime: int = None  # seconds
    cputime: int = None  # seconds
    current_time: float = None
    current_steps: float = None

    def save(self):
        return self.__dict__


class CollisionTree:
    def __init__(self):
        self._tree = {}

    def add(self, source1: Particle, source2: Particle, to: Particle, metadata: Dict):
        self._tree[to.hash.value] = {"parents": [source1.hash.value, source2.hash.value], "meta": metadata}

    def save(self):
        return self._tree

    def load(self, tree):
        self._tree = tree


class EnergyConsercation:
    def __init__(self):
        self._dEs = []
        self.initial_energy = None

    def set_initial_energy(self, initial_energy: float) -> None:
        self.initial_energy = initial_energy

    def add_energy_value(self, energy: float) -> None:
        dE = abs((energy - self.initial_energy) / self.initial_energy)
        self._dEs.append(dE)
