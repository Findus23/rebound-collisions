import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List

from rebound import Particle
from scipy.constants import astronomical_unit, year


@dataclass
class ParticleData:
    water_mass_fraction: float
    type: str
    escaped: float = None


@dataclass
class Input:
    alpha: float
    velocity_original: float
    escape_velocity: float
    gamma: float
    projectile_mass: float
    target_water_fraction: float
    projectile_water_fraction: float
    velocity_esc: float = None

    def __post_init__(self):
        self.velocity_esc = self.velocity_si / self.escape_velocity

    @property
    def velocity_si(self):
        return self.velocity_original * astronomical_unit / year


@dataclass
class CollisionMeta:
    collision_velocities: Tuple[List[float], List[float]] = None
    interpolation_input: List[float] = None
    raw_water_retention: float = None
    raw_mass_retention: float = None
    water_retention: float = None
    mass_retention: float = None
    total_mass: float = None
    final_wmf: float = None
    final_radius: float = None
    target_wmf: float = None
    projectile_wmf: float = None
    time: float = None
    input: Input = None
    adjusted_input: Input = None


@dataclass
class Meta:
    tmax: float = None
    num_savesteps: int = None
    per_savestep: float = None
    initial_N: int = None
    initial_N_planetesimal: int = None
    initial_N_embryo: int = None
    walltime: int = None  # seconds
    cputime: int = None  # seconds
    current_time: float = None
    perfect_merging: bool = None

    def save(self):
        return self.__dict__


class CollisionTree:
    def __init__(self):
        self._tree = {}

    def add(self, source1: Particle, source2: Particle, to: Particle, metadata: CollisionMeta):
        data = {"parents": [source1.hash.value, source2.hash.value], "meta": metadata}
        self._tree[to.hash.value] = data

    def save(self):
        savetree = {}
        tmpcopy = deepcopy(self._tree)
        for key, data in tmpcopy.items():
            metadata = data["meta"]
            metadata = metadata.__dict__
            metadata["input"] = metadata["input"].__dict__
            metadata["adjusted_input"] = metadata["adjusted_input"].__dict__
            data["meta"] = metadata
            savetree[key] = data
        return savetree

    def load(self, tree):
        self._tree = {}
        for key, data in tree.items():
            metadata = data["meta"]
            metadata["input"] = Input(**metadata["input"])
            metadata["adjusted_input"] = Input(**metadata["adjusted_input"])

            metadata = CollisionMeta(**metadata)
            data["meta"] = metadata
            self._tree[key] = data

    def get_tree(self) -> Dict:
        return self._tree


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
        return self.pdata[particle.hash.value]
