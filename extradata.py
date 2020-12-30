import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

from rebound import Particle


@dataclass
class ParticleData:
    water_mass_fraction: float
    type: str
    active: bool = True


@dataclass
class Meta:
    tmax: float = None
    num_savesteps: int = None
    per_savestep: float = None
    max_n: int = None
    walltime: int = None  # seconds
    cputime: int = None  # seconds
    current_time: float = None
    perfect_merging: bool = None

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


class EnergyConservation:
    def __init__(self):
        self._dEs = []
        self.initial_energy = None

    def set_initial_energy(self, initial_energy: float) -> None:
        self.initial_energy = initial_energy

    def add_energy_value(self, energy: float) -> None:
        dE = abs((energy - self.initial_energy) / self.initial_energy)
        self._dEs.append(dE)

    def save(self):
        return {
            "initial_energy": self.initial_energy,
            "dEs": self._dEs
        }

    def load(self, data):
        self.initial_energy = data["initial_energy"]
        self._dEs = data["dEs"]


class ExtraData:

    def __init__(self):
        self.tree = CollisionTree()
        self.pdata: Dict[int, ParticleData] = {}
        self.meta = Meta()
        self.energy = EnergyConservation()

    def save(self, filename: Path):
        pdata = {}
        for k, v in self.pdata.items():
            pdata[k] = v.__dict__

        with filename.open("w") as f:
            json.dump({
                "meta": self.meta.save(),
                "pdata": pdata,
                "tree": self.tree.save(),
                "energy": self.energy.save()
            }, f, indent=2)

    @classmethod
    def load(cls, filename: Path):
        with filename.open() as f:
            data = json.load(f)
        self = cls()
        self.meta = Meta(**data["meta"])

        self.tree.load(data["tree"])
        self.energy.load(data["energy"])
        # self.tree._dEs = data["dEs"]

        for k, v in data["pdata"].items():
            self.pdata[int(k)] = ParticleData(**v)

        return self

    def pd(self, particle: Particle) -> ParticleData:
        return self.pdata[particle.hash]
