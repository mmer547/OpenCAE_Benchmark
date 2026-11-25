"""Recreate the benchmark input decks from the published Altair 2024 examples."""
from pathlib import Path, PurePosixPath
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO
import argparse
import hashlib
import json
import re

BASE = Path(__file__).resolve().parent
RELEASE_URL = "https://2024.help.altair.com/2024/simulation/demos/hwsolvers/radioss/"
def fields(deck, card):
    block = re.search(r"^" + re.escape(card) + r"\s*\n(.*?)(?=^/|\Z)", deck, re.M | re.S)[1]
    return [line for line in block.splitlines() if line.strip() and not line.startswith("#")]


def instrument(deck, case):
    if case.startswith("10_Bending"):
        nodes = [line.split() for line in fields(deck, "/NODE")]
        # Follow one longitudinal edge, including its fixed and rotating ends.
        xmax = max(float(row[1]) for row in nodes)
        ids = [int(row[0]) for row in nodes if float(row[1]) == xmax]
        title, variables = "Edge", "DY        DZ"
    elif case.startswith("0110"):
        ids = [int(n) for row in fields(deck, "/GRNOD/NODE/31")[1:] for n in row.split()]
        title, variables = "DrivenFace", "REACY"
    else:
        return deck
    group = ["/TH/NODE/110", title, variables] + [f"{n:10d}{0:10d}{title}{n}" for n in ids]
    return deck.replace("/END", "\n".join(group) + "\n/END", 1)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["impact", "friction", "pure-bending"])
    parser.add_argument("--archive",type=Path,help="Use an already downloaded official archive")
    parser.add_argument("--fine-impact",action="store_true")
    args = parser.parse_args()
    if args.fine_impact and args.case != "impact":
        parser.error("--fine-impact requires impact")
    case = BASE / args.case
    spec = json.loads((case/"input-source.json").read_text())
    if args.archive:
        data = args.archive.read_bytes()
    else:
        with urlopen(RELEASE_URL + spec["archive"], timeout=60) as response:
            data = response.read()
    if hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise ValueError("Unexpected archive checksum")
    output = case/("inputs-fine" if args.fine_impact else "inputs")
    with ZipFile(BytesIO(data)) as archive:
        for member in archive.namelist():
            path = PurePosixPath(member)
            if path.suffix.lower() != ".rad" or path.parts[0] != spec["directory"]:
                continue
            if ".." in path.parts or path.is_absolute():
                raise ValueError(member)
            target = output.joinpath(*path.parts[1:])
            target.parent.mkdir(parents=True,exist_ok=True)
            deck = archive.read(member).decode("latin1").replace("\r\n","\n")
            if path.name.endswith("_0000.rad"):
                deck = instrument(deck,spec["directory"])
            elif args.fine_impact:
                deck,a = re.subn(r"(?m)^\s*0\.02\s*$","0.0022",deck)
                deck,b = re.subn(r"(?m)^0\.02e-3\s*$","0.5e-6",deck)
                if (a,b)!=(1,1): raise ValueError("Unexpected Engine control deck")
            elif args.case == "friction":
                deck += "\n/TH/TITLE\n"
            target.write_bytes(deck.encode("latin1"))
    count = len(list(output.rglob("*_0000.rad")))
    if count != {"impact":3,"friction":1,"pure-bending":14}[args.case]:
        raise ValueError("Incomplete input set")
    print(count, "Starter / Engine pairs in", output)
