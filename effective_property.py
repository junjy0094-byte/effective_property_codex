"""Compute effective composite properties using PyMAPDL.

This script builds a 1x1x1 mm representative volume element (RVE) with a
central cylindrical fiber sized to achieve a 10% volume fraction. The
model is solved in multiple load cases to extract effective elastic
constants and coefficients of thermal expansion (CTE).

The script is intentionally verbose so it can serve as a template for
further customization in Ansys MAPDL / Material Designer workflows.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Tuple

from ansys.mapdl.core import launch_mapdl

# Basic geometric defaults (millimeters)
CUBE_EDGE = 1.0
VOLUME_FRACTION = 0.10


def fiber_radius_from_volume_fraction(
    vf: float, cube_edge: float, fiber_height: float
) -> float:
    """Return the fiber radius (mm) needed for a target volume fraction.

    Volume fraction is defined as ``fiber_volume / cube_volume``. For a
    vertical cylinder, ``fiber_volume = pi * r**2 * height``.
    """

    cube_volume = cube_edge**3
    fiber_volume = vf * cube_volume
    return math.sqrt(fiber_volume / (math.pi * fiber_height))


@dataclass
class IsotropicMaterial:
    name: str
    ex: float
    nu: float
    rho: float = 1.0
    cte: float = 0.0
    density_unit: str = ""
    _mapdl_id: int = field(init=False, default=-1)

    def register(self, mapdl) -> int:
        self._mapdl_id = mapdl.get("MID", "", "NEXT") or 1
        mapdl.mp("EX", self._mapdl_id, self.ex)
        mapdl.mp("PRXY", self._mapdl_id, self.nu)
        mapdl.mp("ALPX", self._mapdl_id, self.cte)
        if self.density_unit:
            mapdl.mp("DENS", self._mapdl_id, self.rho)
        return self._mapdl_id


@dataclass
class OrthotropicMaterial:
    name: str
    ex: float
    ey: float
    ez: float
    nuxy: float
    nuyz: float
    nuzx: float
    gxy: float
    gyz: float
    gzx: float
    ctex: float = 0.0
    ctey: float = 0.0
    ctez: float = 0.0
    rho: float = 1.0
    density_unit: str = ""
    _mapdl_id: int = field(init=False, default=-1)

    def register(self, mapdl) -> int:
        self._mapdl_id = mapdl.get("MID", "", "NEXT") or 1
        mapdl.mp("EX", self._mapdl_id, self.ex)
        mapdl.mp("EY", self._mapdl_id, self.ey)
        mapdl.mp("EZ", self._mapdl_id, self.ez)
        mapdl.mp("PRXY", self._mapdl_id, self.nuxy)
        mapdl.mp("PRYZ", self._mapdl_id, self.nuyz)
        mapdl.mp("PRZX", self._mapdl_id, self.nuzx)
        mapdl.mp("GXY", self._mapdl_id, self.gxy)
        mapdl.mp("GYZ", self._mapdl_id, self.gyz)
        mapdl.mp("GZX", self._mapdl_id, self.gzx)
        mapdl.mp("ALPX", self._mapdl_id, self.ctex)
        mapdl.mp("ALPY", self._mapdl_id, self.ctey)
        mapdl.mp("ALPZ", self._mapdl_id, self.ctez)
        if self.density_unit:
            mapdl.mp("DENS", self._mapdl_id, self.rho)
        return self._mapdl_id


def build_rve(mapdl, vf: float = VOLUME_FRACTION, cube_edge: float = CUBE_EDGE) -> Tuple[int, int]:
    """Create the matrix and fiber volumes, returning their volume numbers."""

    mapdl.clear()
    mapdl.prep7()

    # Create the cube that hosts the inclusion
    cube_vnum = int(mapdl.blc4(0, 0, cube_edge, cube_edge, cube_edge))

    # Create the cylindrical fiber aligned to the Z-axis
    radius = fiber_radius_from_volume_fraction(vf, cube_edge, cube_edge)
    # CYL4 expects positional arguments; provide radius for both ends and full depth
    fiber_vnum = int(
        mapdl.cyl4(cube_edge / 2, cube_edge / 2, radius, radius, 0, 360, cube_edge)
    )

    # Subtract fiber from cube to obtain the matrix volume
    mapdl.vsbv(cube_vnum, fiber_vnum)
    return cube_vnum, fiber_vnum


def assign_materials(
    mapdl,
    matrix_volume: int,
    fiber_volume: int,
    matrix: IsotropicMaterial,
    fiber: OrthotropicMaterial,
    element_size: float = 0.1,
) -> None:
    """Mesh the RVE and assign materials to volumes."""

    mid_matrix = matrix.register(mapdl)
    mid_fiber = fiber.register(mapdl)

    mapdl.esize(element_size)

    mapdl.vsel("S", "VOLU", vmin=matrix_volume, vmax=matrix_volume)
    mapdl.vatt(mid_matrix, 1, 1, 0)
    mapdl.allsel()

    mapdl.vsel("S", "VOLU", vmin=fiber_volume, vmax=fiber_volume)
    mapdl.vatt(mid_fiber, 1, 1, 0)
    mapdl.allsel()

    mapdl.vmesh("ALL")


def _begin_static_solution(mapdl) -> None:
    """Reset boundary conditions and start a static solution block."""

    mapdl.finish()
    mapdl.run("/SOLU")
    mapdl.antype("STATIC")
    mapdl.allsel()
    mapdl.ddele("ALL", "ALL")
    mapdl.fdele("ALL", "ALL")


def _solve_uniaxial(mapdl, axis: str, displacement: float, cube_edge: float) -> float:
    comp = {"X": "UX", "Y": "UY", "Z": "UZ"}[axis]

    _begin_static_solution(mapdl)

    mapdl.allsel()
    mapdl.nsel("S", "LOC", axis, 0.0)
    mapdl.d("ALL", "UX", 0.0)
    mapdl.d("ALL", "UY", 0.0)
    mapdl.d("ALL", "UZ", 0.0)

    mapdl.nsel("S", "LOC", axis, cube_edge)
    mapdl.d("ALL", comp, displacement)

    mapdl.solve()
    mapdl.finish()

    mapdl.post1()
    mapdl.set(1)
    mapdl.nsel("S", "LOC", axis, cube_edge)
    reaction = float(mapdl.fsum(return_value=True)[comp])
    avg_stress = reaction / (cube_edge**2)
    avg_strain = displacement / cube_edge
    return avg_stress / avg_strain


def _solve_shear(mapdl, primary_axis: str, secondary_axis: str, shear_disp: float, cube_edge: float) -> float:
    """Apply a tangential displacement on the face normal to ``primary_axis``."""

    comp = {"XY": "UY", "YZ": "UZ", "ZX": "UX"}[f"{primary_axis}{secondary_axis}"]
    normal_comp = {"X": "UX", "Y": "UY", "Z": "UZ"}[primary_axis]

    _begin_static_solution(mapdl)

    mapdl.allsel()
    mapdl.nsel("S", "LOC", primary_axis, 0.0)
    mapdl.d("ALL", "UX", 0.0)
    mapdl.d("ALL", "UY", 0.0)
    mapdl.d("ALL", "UZ", 0.0)

    mapdl.nsel("S", "LOC", primary_axis, cube_edge)
    mapdl.d("ALL", comp, shear_disp)
    mapdl.d("ALL", normal_comp, 0.0)

    mapdl.solve()
    mapdl.finish()

    mapdl.post1()
    mapdl.set(1)
    mapdl.nsel("S", "LOC", primary_axis, cube_edge)
    reaction = float(mapdl.fsum(return_value=True)[comp])
    avg_shear_stress = reaction / (cube_edge**2)
    avg_shear_strain = shear_disp / cube_edge
    return avg_shear_stress / avg_shear_strain


def _solve_cte(mapdl, cube_edge: float, delta_t: float = 1.0) -> Dict[str, float]:
    _begin_static_solution(mapdl)
    mapdl.bf("ALL", "TEMP", delta_t)

    # Fix one corner to remove rigid body motions
    mapdl.nsel("S", "LOC", "X", 0.0)
    mapdl.nsel("R", "LOC", "Y", 0.0)
    mapdl.nsel("R", "LOC", "Z", 0.0)
    mapdl.d("ALL", "ALL", 0)
    mapdl.allsel()

    mapdl.solve()
    mapdl.finish()

    mapdl.post1()
    mapdl.set(1)

    ctes = {}
    for axis, comp in ("X", "UX"), ("Y", "UY"), ("Z", "UZ"):
        mapdl.nsel("S", "LOC", axis, cube_edge)
        disp = float(mapdl.get("_CTE", "NODE", 0, "U", comp[-1]))
        ctes[f"CTE{axis}"] = disp / (cube_edge * delta_t)
    mapdl.allsel()
    return ctes


def compute_effective_properties() -> Dict[str, float]:
    mapdl = launch_mapdl()
    try:
        cube_v, fiber_v = build_rve(mapdl)

        # Example material data (edit as needed)
        epoxy = IsotropicMaterial(name="Epoxy", ex=3.5e3, nu=0.35, cte=6.5e-6)
        carbon = OrthotropicMaterial(
            name="Carbon Fiber",
            ex=240e3,
            ey=15e3,
            ez=15e3,
            nuxy=0.2,
            nuyz=0.25,
            nuzx=0.25,
            gxy=5e3,
            gyz=5e3,
            gzx=5e3,
            ctex=-0.5e-6,
            ctey=15e-6,
            ctez=15e-6,
        )

        assign_materials(mapdl, cube_v, fiber_v, epoxy, carbon)

        # Compute elastic modulus
        ex = _solve_uniaxial(mapdl, "X", displacement=1e-3, cube_edge=CUBE_EDGE)
        ey = _solve_uniaxial(mapdl, "Y", displacement=1e-3, cube_edge=CUBE_EDGE)
        ez = _solve_uniaxial(mapdl, "Z", displacement=1e-3, cube_edge=CUBE_EDGE)

        # Compute shear modulus
        gxy = _solve_shear(mapdl, "X", "Y", shear_disp=1e-3, cube_edge=CUBE_EDGE)
        gyz = _solve_shear(mapdl, "Y", "Z", shear_disp=1e-3, cube_edge=CUBE_EDGE)
        gzx = _solve_shear(mapdl, "Z", "X", shear_disp=1e-3, cube_edge=CUBE_EDGE)

        # Poisson's ratios estimated from lateral reactions in uniaxial tests
        def _poisson(primary: str, lateral: str, disp: float) -> float:
            comp = {"X": "UX", "Y": "UY", "Z": "UZ"}[lateral]
            mapdl.post1()
            mapdl.set(1)
            mapdl.nsel("S", "LOC", primary, CUBE_EDGE)
            lateral_disp = float(mapdl.get("_POIS", "NODE", 0, "U", comp[-1]))
            return -(lateral_disp / CUBE_EDGE) / (disp / CUBE_EDGE)

        vxy = _poisson("X", "Y", 1e-3)
        vyz = _poisson("Y", "Z", 1e-3)
        vzx = _poisson("Z", "X", 1e-3)

        # Thermal expansion
        cte = _solve_cte(mapdl, cube_edge=CUBE_EDGE, delta_t=1.0)

        return {
            "Ex": ex,
            "Ey": ey,
            "Ez": ez,
            "Gxy": gxy,
            "Gyz": gyz,
            "Gzx": gzx,
            "vxy": vxy,
            "vyz": vyz,
            "vzx": vzx,
            "CTEx": cte["CTEX"],
            "CTEy": cte["CTEY"],
            "CTEz": cte["CTEZ"],
        }
    finally:
        mapdl.exit()


if __name__ == "__main__":
    properties = compute_effective_properties()
    for name, value in properties.items():
        print(f"{name}: {value:.4e}")
