"""Sweep a velocity-driven rigid duct with ElmerSolver and compare inlet pressure."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
LENGTH, HEIGHT = 10.0, 0.1
DENSITY, SOUND_SPEED, VELOCITY = 1.2, 343.0, 1.0
NY = 2


def write_mesh(folder, nx):
    folder.mkdir()
    nodes, elements, edges = [], [], []

    def node(i, j):
        return 1 + i + (nx + 1) * j

    def cell(i, j):
        return 1 + i + nx * j

    def boundary(tag, parent, a, b):
        edges.append(f"{len(edges)+1} {tag} {parent} 0 202 {a} {b}")

    for j in range(NY + 1):
        for i in range(nx + 1):
            nodes.append(f"{node(i,j)} -1 {LENGTH*i/nx:.12g} {HEIGHT*j/NY:.12g} 0")
    for j in range(NY):
        for i in range(nx):
            a, b = node(i, j), node(i+1, j)
            c, d = node(i+1, j+1), node(i, j+1)
            parent = cell(i, j)
            elements.append(f"{parent} 1 404 {a} {b} {c} {d}")
            if i == 0:
                boundary(1, parent, d, a)
            if i == nx-1:
                boundary(2, parent, b, c)
            if j == 0:
                boundary(3, parent, a, b)
            if j == NY-1:
                boundary(4, parent, c, d)
    header = [f"{len(nodes)} {len(elements)} {len(edges)}", "2",
              f"404 {len(elements)}", f"202 {len(edges)}"]
    for suffix, lines in (("nodes", nodes), ("elements", elements),
                          ("boundary", edges), ("header", header)):
        (folder / f"mesh.{suffix}").write_text("\n".join(lines) + "\n", encoding="ascii")


def inlet_pressure(vtu):
    piece = ET.parse(vtu).find(".//Piece")
    coords = list(map(float, piece.find("Points/DataArray").text.split()))
    inlet = sorted((coords[3*i+1], i) for i in range(len(coords)//3)
                   if abs(coords[3*i]) < 1e-10)
    if len(inlet) != NY + 1:
        raise ValueError(f"Unexpected inlet nodes: {inlet}")
    arrays = {a.get("Name").lower(): a for a in piece.findall("PointData/DataArray")}
    def average(key):
        values = list(map(float, arrays[key].text.split()))
        return sum((y2-y1)*(values[i]+values[k])/2
                   for (y1,i),(y2,k) in zip(inlet, inlet[1:])) / HEIGHT
    return complex(average("pressure wave 1"), average("pressure wave 2"))


def run(mesh_sizes):
    solver = shutil.which("ElmerSolver")
    if not solver:
        raise SystemExit("Add ElmerSolver to PATH first.")
    template = (ROOT / "duct.sif.in").read_text(encoding="ascii")
    with tempfile.TemporaryDirectory(dir=ROOT) as temp:
        work = Path(temp)
        for nx in mesh_sizes:
            case = work / str(nx)
            case.mkdir()
            write_mesh(case / "mesh", nx)
            rows = []
            for frequency in range(5, 101):
                omega = 2 * math.pi * frequency
                sif = template
                for key, val in (("OMEGA", omega), ("INLET_FLUX", omega*DENSITY*VELOCITY),
                                 ("SOUND_SPEED", SOUND_SPEED), ("DENSITY", DENSITY)):
                    sif = sif.replace(f"@{key}@", format(val, ".17g"))
                (case / "case.sif").write_text(sif, encoding="ascii")
                output = subprocess.run([solver, "case.sif"], cwd=case, capture_output=True,
                                        env={**os.environ, "OMP_NUM_THREADS":"1",
                                             "OPENBLAS_NUM_THREADS":"1"})
                log = (output.stdout + output.stderr).decode(errors="replace")
                if output.returncode or "ELMER SOLVER FINISHED" not in log:
                    raise RuntimeError(log[-6000:])
                vtu, = list(case.rglob("*.vtu"))
                p = inlet_pressure(vtu)
                ref = -DENSITY*SOUND_SPEED*VELOCITY/math.tan(omega*LENGTH/SOUND_SPEED)
                error = abs(p - complex(0, ref))
                rows.append((frequency, p.real, p.imag, abs(p), ref, abs(ref),
                             error, 100*error/abs(ref)))
                if nx == 1600 and frequency == 69:
                    shutil.copyfile(case / "case.sif", ROOT / "case.sif")
                    mesh_out = ROOT / "mesh"
                    mesh_out.mkdir(exist_ok=True)
                    for mesh_file in (case / "mesh").glob("mesh.*"):
                        shutil.copyfile(mesh_file, mesh_out / mesh_file.name)
                    shutil.copyfile(vtu, ROOT / "pressure_69Hz.vtu")
                    (ROOT / "solver_69Hz.log").write_text("\n".join(log.splitlines()) + "\n",
                                                         encoding="utf-8")
                vtu.unlink()
            csv_path = ROOT / f"response_{nx}.csv"
            with csv_path.open("w", newline="", encoding="ascii") as file:
                writer = csv.writer(file)
                writer.writerow(("frequency_Hz", "elmer_real_Pa", "elmer_imag_Pa",
                                 "elmer_abs_Pa", "analytic_imag_Pa", "analytic_abs_Pa",
                                 "complex_error_Pa", "relative_error_pct"))
                writer.writerows(rows)
            print(f"Wrote {csv_path.name}, {len(rows)} frequencies", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", nargs="+", type=int, choices=(200, 400, 1600),
                        default=(200, 400, 1600), help="axial mesh sizes to sweep")
    run(parser.parse_args().nx)
