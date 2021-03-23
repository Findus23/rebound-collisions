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
    core_mass_fraction: float
    type: str
    escaped: float = None
    collided_with_sun: float = None
    total_mass: float = None

    @property
    def water_mass(self) -> float:
        return self.total_mass * self.water_mass_fraction

    @property
    def mantle_mass_fraction(self) -> float:
        return 1 - self.core_mass_fraction - self.water_mass_fraction


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
    collision_positions: Tuple[List[float], List[float]] = None
    collision_radii: Tuple[float, float] = None
    interpolation_input: List[float] = None
    raw_water_retention: float = None
    raw_mantle_retention: float = None
    raw_core_retention: float = None
    water_retention: float = None
    mantle_retention: float = None
    core_retention: float = None
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
    initcon_file: str = None
    tmax: float = None
    num_savesteps: int = None
    per_savestep: float = None
    initial_N: int = None
    initial_N_planetesimal: int = None
    initial_N_embryo: int = None
    walltime: int = None  # seconds
    cputime: int = None  # seconds
    current_time: float = None
    hash_counter: int = 0
    git_hash: str = None
    rebound_hash: str = None
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
            if metadata["input"]:
                metadata["input"] = metadata["input"].__dict__
                metadata["adjusted_input"] = metadata["adjusted_input"].__dict__
            data["meta"] = metadata
            savetree[key] = data
        return savetree

    def load(self, tree):
        self._tree = {}
        for key, data in tree.items():
            metadata = data["meta"]
            if metadata["input"]:
                metadata["input"] = Input(**metadata["input"])
                metadata["adjusted_input"] = Input(**metadata["adjusted_input"])

            metadata = CollisionMeta(**metadata)
            data["meta"] = metadata
            self._tree[int(key)] = data

    def get_tree(self) -> Dict:
        return self._tree

    def get(self, particle: Particle):
        return self._tree[particle.hash.value]


class History:
    def __init__(self):
        self.energy = []
        self.momentum = []
        self.total_mass = []
        self.time = []
        self.N = []
        self.N_active = []

    def append(self, energy: float, momentum: float, total_mass: float, time: float, N: int, N_active: int):
        self.energy.append(energy)
        self.momentum.append(momentum)
        self.total_mass.append(total_mass)
        self.time.append(time)
        self.N.append(N)
        self.N_active.append(N_active)

    def save(self):
        return self.__dict__

    def load(self, data):
        self.__dict__ = data


class ExtraData:

    def __init__(self):
        self.tree = CollisionTree()
        self.pdata: Dict[int, ParticleData] = {}
        self.meta = Meta()
        self.history = History()

    def save(self, base_filename: Path):
        pdata = {}
        for k, v in self.pdata.items():
            pdata[k] = v.__dict__

        with base_filename.with_suffix(".extra.json").open("w") as f:
            json.dump({
                "meta": self.meta.save(),
                "pdata": pdata,
                "tree": self.tree.save(),
            }, f, indent=2)
        with base_filename.with_suffix(".history.json").open("w") as f:
            json.dump(self.history.save(), f, indent=2)

    @classmethod
    def load(cls, base_filename: Path):
        with base_filename.with_suffix(".extra.json").open() as f:
            data = json.load(f)
        with base_filename.with_suffix(".history.json").open() as f:
            history = json.load(f)
        self = cls()
        self.meta = Meta(**data["meta"])
        self.history = History()
        self.history.load(history)
        self.tree.load(data["tree"])

        for k, v in data["pdata"].items():
            self.pdata[int(k)] = ParticleData(**v)

        return self

    def pd(self, particle: Particle) -> ParticleData:
        return self.pdata[particle.hash.value]
