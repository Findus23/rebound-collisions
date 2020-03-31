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

    def save(self, filename: Path):
        pdata = {}
        for k, v in self.pdata.items():
            pdata[k] = v.__dict__

        with filename.open("w") as f:
            json.dump({
                "meta": self.meta.save(),
                "pdata": pdata,
                "tree": self.tree.save()
            }, f, indent=2)

    @classmethod
    def load(cls, filename: Path):
        with filename.open() as f:
            data = json.load(f)
        self = cls()
        self.meta = Meta(**data["meta"])

        self.tree.load(data["tree"])

        for k, v in data["pdata"].items():
            self.pdata[k] = ParticleData(**v)

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

    def add(self, source1: Particle, source2: Particle, to: Particle):
        self._tree[to.hash.value] = [source1.hash.value, source2.hash.value]

    def save(self):
        return self._tree

    def load(self, tree):
        self._tree = tree
