"""Rebuild a case's initial condition from its predecessor's local solution."""
import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess


def read_field(folder, name):
    path = folder / name
    if path.is_file():
        return path.read_text(encoding="ascii")
    with gzip.open(str(path) + ".gz", "rt", encoding="ascii") as stream:
        return stream.read()


def values(text):
    import numpy as np
    match = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if match is None:
        raise ValueError("Expected a nonuniform ASCII field")
    array = np.fromstring(match[3].replace("(", " ").replace(")", " "), sep=" ")
    width = 3 if match[1] == "vector" else 1
    return array.reshape(int(match[2]), width)


def centres(case, time):
    folder = case / time
    if not any((folder / n).is_file() for n in ("C", "C.gz")):
        subprocess.run(["postProcess", "-case", str(case), "-func", "writeCellCentres", "-time", time], check=True)
    return values(read_field(folder, "C"))


def internal_block(a, digits):
    vector = a.shape[1] == 3
    fmt = "%." + str(digits) + "g"
    rows = ["(" + " ".join(fmt % v for v in row) + ")" for row in a] if vector else [fmt % v for v in a[:, 0]]
    kind = "vector" if vector else "scalar"
    return f"internalField nonuniform List<{kind}>\n{len(a)}\n(\n" + "\n".join(rows) + "\n);"


def prepare(case, selected_time=None):
    spec = json.loads((case/"restart.json").read_text())
    source = (case/spec["source"]).resolve()
    time = selected_time or spec.get("time")
    if time is None:
        times = [p.name for p in source.iterdir() if p.is_dir() and re.fullmatch(r"[0-9.eE+-]+", p.name) and float(p.name) > 0]
        if not times:
            raise ValueError(f"Run the predecessor first: {source}")
        time = max(times, key=float)
    target = case/"0"
    target.mkdir(exist_ok=True)
    method = spec["method"]
    if method in ("nearest3d", "idw2d"):
        import numpy as np
        from scipy.spatial import cKDTree
        old, new = centres(source,time), centres(case,"0")
        if method == "idw2d":
            distance, indices = cKDTree(old[:,:2]).query(new[:,:2],k=8)
            weights = 1 / np.maximum(distance, 1e-13)**2
            weights /= weights.sum(axis=1,keepdims=True)
        else:
            _, indices = cKDTree(old).query(new)
    ready = {}
    for name in spec["fields"]:
        text = read_field(source/time,name)
        if method == "copy":
            ready[name] = re.sub(r'(\blocation\s+)"[^"]*"', r'\1"0"', text)
            continue
        if method == "internal":
            block = text[text.index("internalField"):text.index("boundaryField")].strip()
        else:
            a = values(text)
            if method == "scale":
                a = a * {"U":5,"p":25}[name]
            elif method == "idw2d":
                a = np.einsum("ij,ijk->ik", weights, a[indices])
            elif method == "nearest3d":
                a = a[indices]
            else:
                raise ValueError(method)
            block = internal_block(a,11 if method=="nearest3d" else 12)
        template = (case/"initial-templates"/name).read_text()
        if template.count("@INTERNAL_FIELD@") != 1:
            raise ValueError(name)
        ready[name] = template.replace("@INTERNAL_FIELD@", block)
    for name, text in ready.items():
        (target/name).write_text(text,encoding="ascii")
    print(f"{source} at t={time} -> {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("--time",help="Select a reconstructed source time directory")
    args = parser.parse_args()
    prepare(args.case.resolve(),args.time)
