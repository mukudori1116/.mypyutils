import os
from pathlib import Path
import subprocess
from shutil import copy, copytree

root = Path(os.path.dirname(os.path.abspath(__file__)))
modify = root / "modify"
modify.mkdir(exist_ok=True)

raw_ligand_files = Path("./ligands/").glob("*")
for i, ligand in enumerate(raw_ligand_files):
    os.chdir(str(root))
    plig = str(ligand)
    name = ligand.stem

    # antechamber
    output = modify / (name + ".mol2")
    if not output.exists():
        subprocess.check_call([
            "antechamber",
            "-i", plig, "-fi", "pdb",
            "-o", str(output), "-fo", "mol2",
            "-at", "gaff2",
            "-c", "bcc",
            "-pf", "y",
            "-rn", "LIG"
        ])
    else:
        print(f"{output.name} has already exist")

    # parmchk2
    frcmod = modify/(name + ".frcmod")
    if not frcmod.exists():
        subprocess.check_call([
            "parmchk2",
            "-i", output, "-f", "mol2",
            "-o", str(frcmod),
        ])
    else:
        print(f"{frcmod.name} has already exist")

    # tleap
    topology = root / "top"
    topology.mkdir(exist_ok=True)
    leapin = topology/f"leap{i+1}.in"
    parm7 = topology/f"com{i+1}.parm7"
    rst7 = topology/f"com{i+1}.rst7"
    pdb = topology/f"com{i+1}.pdb"
    with open("./src/leap.in") as f:
        leapin_template = f.read()
    leapin_txt = leapin_template.replace("{{frcmod}}", str(frcmod)) \
        .replace("{{mol2}}", str(output)) \
        .replace("{{parm7}}", str(parm7)) \
        .replace("{{rst7}}", str(rst7)) \
        .replace("{{pdb}}", str(pdb))
    with open(str(leapin), "w") as f:
        f.write(leapin_txt)

    if not pdb.exists():
        subprocess.check_call(['tleap', '-f', leapin])

    # # move topology
    gro_dir = root / f"gromacs{i+1}"
    if not gro_dir.exists():
        copytree("./src/gromacs", str(gro_dir))
    grotop = gro_dir/"top"
    minimize = gro_dir/"minimize"
    heat = gro_dir/"heat"
    pr = gro_dir/"pr"

    if not (grotop/"com.rst7").exists():
        copy(str(parm7), grotop/"com.parm7")
        copy(str(rst7), grotop/"com.rst7")

    # convert topology
    if not (grotop/"com.gro").exists():
        os.chdir(str(grotop))
        subprocess.check_call(
            ["bash", "groconvert.sh", "-i", "com", "-o", "com", "-r"])
    if not (grotop/"index.ndx").exists():
        os.chdir(str(grotop))
        subprocess.check_call(
            f"gmx_mpi make_ndx -f com.gro < " + str(root/"index.in"), shell=True
        )

    os.chdir(str(minimize))
    subprocess.check_call("bash run.sh", shell=True)
    os.chdir(str(heat))
    subprocess.check_call("bash run.sh", shell=True)
    os.chdir(str(heat))
    subprocess.check_call("bash run.sh", shell=True)
