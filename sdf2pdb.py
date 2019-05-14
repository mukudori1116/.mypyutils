import sys
import pybel
import openbabel

mols = pybel.readfile("sdf", sys.argv[1])
convert_format = input("Convert and split into: ").strip()
conv = openbabel.OBConversion()
conv.SetInAndOutFormats("sdf", convert_format)

for i, mol in enumerate(mols):
    obmol = mol.OBMol
    fname = f"lig{i+1}.{convert_format}"
    conv.WriteFile(obmol, f"ligands/{fname}")
    print(f"Generate {fname}")
print("Done")
