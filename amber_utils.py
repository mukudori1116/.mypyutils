import os
import re
from glob import glob
import numpy as np
import pandas as pd
from myutils import grab


def split_decomp(file: str):
    dir = os.path.dirname(file)
    base = os.path.basename(file)
    name, ext = base.split(".")
    os.makedirs(f'{dir}/decomp', exist_ok=True)
    with open(file) as f:
        loop = 0
        line = f.readline()
        while line:
            prev_list_length = 0
            stop = False
            while line:
                if stop:
                    loop += 1
                    break
                text = ""
                for _ in range(1000000):
                    title = re.match(r"^.+\:", line)
                    if title:
                        text += title.string
                    elif len(line.split(",")) == 0:
                        pass
                    elif len(line.split(",")) >= prev_list_length:
                        text += line
                        prev_list_length = len(line.split(","))
                    else:
                        stop = True
                        break
                    line = f.readline()
                with open(f"{dir}/decomp/{name}.{loop}", "a") as g:
                    g.write(text)

    return glob(f"{dir}/decomp/*")


def residue_map(topology: str, Cterm=0):
    with open(topology) as f:
        lines = list(map(lambda x: x.strip(), f.readlines()))
        start = lines.index("%FLAG RESIDUE_LABEL") + 2
        end = lines.index("%FLAG RESIDUE_POINTER")
        residues = " ".join(lines[start:end]).split()
        resmap = pd.concat([
            pd.Series(range(Cterm, len(residues) + Cterm),
                      dtype=np.uint16).rename('#'),
            pd.Series(residues, dtype=np.dtype('U3')).rename('Residue')
        ], axis=1)
    return resmap


def get_interaction(reader, residue: int=1,):
    data = pd.DataFrame()
    for df in reader:
        splited1 = df['Resid 1'].str.split(expand=True)
        df['Residue1'] = splited1.iloc[:, 0].astype('U3')
        df['#1'] = splited1.iloc[:, 1].astype('uint16')

        splited2 = df['Resid 2'].str.split(expand=True)
        df['Residue2'] = splited2.iloc[:, 0].astype('U3')
        df['#2'] = splited2.iloc[:, 1].astype('uint16')
        df_fix = df.drop(['Resid 1', 'Resid 2', 'Internal'],
                         axis=1).iloc[:, [0, 6, 7, 8, 9, 1, 2, 3, 4, 5]]
        data = data.append(df_fix[df_fix['#1'] == residue])

    info = data[data['Frame #'] == 1].iloc[:, [1, 2, 3, 4]]
    nframe = data['Frame #'].max()
    ndata = int(len(data)/nframe)
    arr = data.values[:, 5:].reshape((nframe, ndata, 5)).astype('float64')
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    return (info, mean, std)


def decomp_summary(decomp):

    def fix_df(datalist):
        df = pd.DataFrame(datalist)
        resdf = pd.concat(
            [
                df.iloc[2:, 0].str.split(expand=True).rename(
                    columns={0: "res1", 1: "n1"}).astype({"res1": "U3", "n1": np.uint16}),
                df.iloc[2:, 1].str.split(expand=True).rename(
                    columns={0: "res2", 1: "n2"}).astype({"res2": "U3", "n2": np.uint16})
            ],
            axis=1)
        val_mean_df = df.iloc[2:, 2::3].astype(np.float64)
        val_std_df = df.iloc[2:, 3::3].astype(np.float64)
        val_sem_df = df.iloc[2:, 4::3].astype(np.float64)
        val_mean_df.columns = columns
        val_std_df.columns = columns
        val_sem_df.columns = columns
        mean = pd.concat([resdf, val_mean_df], axis=1).reset_index(drop=True)
        std = pd.concat([resdf, val_std_df], axis=1).reset_index(drop=True)
        sem = pd.concat([resdf, val_sem_df], axis=1).reset_index(drop=True)
        return (mean, std, sem)

    with open(decomp) as f:
        line = f.readline()
        tmp_list = []
        data = {}
        parent, child = None, None
        columns = ['Internal', 'van der Waals', 'Electrostatic',
                   'PolarSolvation', 'Non-PolarSolv', 'TOTAL']
        while line:
            if re.match(r"^.+\:", line):
                if tmp_list:
                    mean, std, sem = fix_df(tmp_list)
                    data[parent][child]['mean'] = mean
                    data[parent][child]['std'] = std
                    data[parent][child]['sem'] = sem
                    tmp_list = []
                line = "".join(line.split(","))
                keys1 = re.match(r"Complex|Receptor|Ligand|DELTA", line)
                keys2 = re.match(r"Total|Sidechain|Backbone", line)
                if keys1:
                    parent = keys1.group(0)
                    data[parent] = {}
                elif keys2:
                    child = keys2.group(0)
                    data[parent][child] = {}
            else:
                if len(line.split(",")) != 1:
                    tmp_list.append(line.split(","))
            line = f.readline()
        mean, std, sem = fix_df(tmp_list)
        data[parent][child]['mean'] = mean
        data[parent][child]['std'] = std
        data[parent][child]['sem'] = sem

    return data


def cutoff_mol2_charges(input_file, output_file):
    header, body, footer = grab(input_file, "@<TRIPOS>ATOM", "@<TRIPOS>BOND")
    atoms = map(lambda x: x.split(), body.split("\n"))
    df = pd.DataFrame(atoms).dropna()

    charges_origin = df.iloc[:, -1].values
    charges_fix = df.iloc[:, -1].astype(np.float128).round(5)
    delta = charges_fix.sum()
    charges_fix[charges_fix == charges_fix.max()] -= delta
    charges_fix_5 = map(lambda x: "{:.5f}".format(x), charges_fix.values)
    for origin, fix in zip(charges_origin, charges_fix_5):
        body = body.replace(origin, fix)
    with open(output_file, "w") as f:
        f.write(header+body+footer)
