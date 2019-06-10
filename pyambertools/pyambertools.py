import os
from pathlib import Path
import subprocess
import shlex
import shutil
from typing import Optional, List
import re

from amber_utils import cutoff_mol2_charges


class AmberToolsWrapper:

    with open("src/leap.in") as f:
        LEAPIN = f.read()

    def __init__(self, path):
        root = Path(path)
        if not root.is_dir():
            raise ValueError("Path must be directory")
        else:
            self.root = Path(path).resolve()
            root.mkdir(exist_ok=True)

        self.src = self.root/"src"
        self.src.mkdir(exist_ok=True)
        self.ligads = self.root/"ligads"
        self.ligads.mkdir(exist_ok=True)
        self.rec = self.root/"rec"
        self.rec.mkdir(exist_ok=True)
        self.modify = self.root / "modify"
        self.modify.mkdir(exist_ok=True)
        self.top = self.root / "top"
        self.top.mkdir(exist_ok=True)
        self.traj = self.root / "traj"
        self.traj.mkdir(exist_ok=True)
        self.gromacs = self.root/"gromacs"
        if not self.gromacs.exists():
            shutil.copytree("./src/gmxconf", str(self.gromacs))

    def antechamber(self, ligand_files: Optional[List[Path]]=None) -> List[Path]:
        outputs = list()
        if ligand_files is None:
            ligand_files = self.ligads.glob("*")
        else:
            for ligand in ligand_files:
                input_file = str(ligand)
                input_format = ligand.suffix
                output_file = str(self.modify/(ligand.stem + ".mol2"))
                cmd = f"antechamber -i {input_file} -fi {input_format} " \
                    + f"-o {output_file} -fo mol2 " \
                    + f"-at gaff2 -c bcc -pf y -rn LIG"
                proc = subprocess.Popen(
                    shlex.split(cmd),
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                # fix charge
                cutoff_mol2_charges(output_file, output_file)
                outputs.append(self.modify/(ligand.stem + ".mol2"))
        return outputs

    def parmchk2(self, fixed_ligands: Optional[List[Path]]=None) -> List[Path]:
        if fixed_ligands is None:
            fixed_ligands = self.modify.glob("*.mol2")
        else:
            for flig in fixed_ligands:
                ligand_parameter = flig.stem + ".frcmod"
                cmd = f"parmchk2 -i {str(flig)} -f mol2 -o {ligand_parameter}"
                proc = subprocess.Popen(
                    shlex.split(cmd),
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
        return list(self.modify.glog("*.frcmod"))

    def generate_leapfile(
        self, frcmod: Path, ligand: Path, output_name: str, output_file_path: Path,
    ) -> Path:
        lig_parm7 = self.top/ligand.stem + ".parm7"
        lig_rst7 = self.top/ligand.stem + ".rst7"
        rec_parm7 = self.top/"rec.parm7"
        rec_rst7 = self.top/"rec.rst7"
        com_nonsolv_parm7 = self.top/(output_name + "_nonsolv.parm7")
        com_nonsolv_rst7 = self.top/(output_name + "_nonsolv.rst7")
        leapin = self.LEAPIN \
            .replace("{{frcmod}}", str(frcmod)) \
            .replace("{{mol2}}", str(ligand)) \
            .replace("{{rec}}", str(self.rec/"rec.pdb")) \
            .replace("{{lig_parm7}}", str(lig_parm7)) \
            .replace("{{lig_rst7}}", str(lig_rst7)) \
            .replace("{{rec_parm7}}", str(rec_parm7)) \
            .replace("{{rec_rst7}}", str(rec_rst7)) \
            .replace("{{parm7}}", output_name + ".parm7") \
            .replace("{{rst7}}", output_name + ".rst7") \
            .replace("{{pdb}}", output_name + ".pdb") \
            .replace("{{com_nonsolv_parm7}}", str(com_nonsolv_parm7)) \
            .replace("{{com_nonsolv_rst7}}", str(com_nonsolv_rst7))
        with output_file_path.open("w") as f:
            f.write(leapin)
        return output_file_path

    def generate_leapfiles(self) -> List[Path]:
        ligands = self.modify.glob("*.mol2")
        for ligand in ligands:
            frcmod = self.modify/(ligand.stem + ".frcmod")
            number = re.match(f"\d+", ligand.stem)
            if number:
                output_name = "com" + number.string
                output_file_path = self.top/f"leap{number.string}.in"
            else:
                raise ValueError("Ligand file name must contain number.")
            self.generate_leapfile(
                frcmod=frcmod, ligand=ligand,
                output_name=output_name, output_file_path=output_file_path)
        return list(self.top.glob("*.in"))

    def tleap(self, leapfiles: List[Path]):
        for leap in leapfiles:
            cmd = f"tleap -f {str(leap)}"
            proc = subprocess.Popen(
                shlex.split(cmd),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))


if __name__ == "__main__":

    # Antechamber
    atw = AmberToolsWrapper(os.getcwd())
    fixed_ligands = atw.antechamber()

    # parmchk2
    ligand_params = atw.parmchk2()

    # tleap
    leap_files = atw.generate_leapfiles()
    tleap = atw.tleap(leap_files)
