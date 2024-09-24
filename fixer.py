"""This script is suppose to fix the broken chain in pdb files.

    Usage
    =====
        ./fixer.py -f <path_directory_pdb_files>
"""

__author__ = "Ravy LEON FOUN LIN"
__date__ = "27-08-2024"
__version__ = "1.0.0"

import numpy as np
from glob import glob
from tqdm import tqdm
import sys
import argparse
from Bio.PDB import PDBParser
from modeller import Environ, log
from modeller.automodel import *
from modeller import *
import pymolPy3

parser = argparse.ArgumentParser()
parser.add_argument("-f", help = "Path to the directory where are found the pdb files.")

args = parser.parse_args()

if args.f == None :
    parser.print_help()
    sys.exit()
else :
    PATH = args.f


acides_amines_3_lettres = [
    "ALA",  # Alanine
    "ARG",  # Arginine
    "ASN",  # Asparagine
    "ASP",  # Acide aspartique
    "CYS",  # Cystéine
    "GLU",  # Acide glutamique
    "GLN",  # Glutamine
    "GLY",  # Glycine
    "HIS",  # Histidine
    "ILE",  # Isoleucine
    "LEU",  # Leucine
    "LYS",  # Lysine
    "MET",  # Méthionine
    "PHE",  # Phénylalanine
    "PRO",  # Proline
    "SER",  # Sérine
    "THR",  # Thréonine
    "TRP",  # Tryptophane
    "TYR",  # Tyrosine
    "VAL"   # Valine
]

acides_amines_dict = {
    "ALA": "A",  # Alanine
    "ARG": "R",  # Arginine
    "ASN": "N",  # Asparagine
    "ASP": "D",  # Acide aspartique
    "CYS": "C",  # Cystéine
    "GLU": "E",  # Acide glutamique
    "GLN": "Q",  # Glutamine
    "GLY": "G",  # Glycine
    "HIS": "H",  # Histidine
    "ILE": "I",  # Isoleucine
    "LEU": "L",  # Leucine
    "LYS": "K",  # Lysine
    "MET": "M",  # Méthionine
    "PHE": "F",  # Phénylalanine
    "PRO": "P",  # Proline
    "SER": "S",  # Sérine
    "THR": "T",  # Thréonine
    "TRP": "W",  # Tryptophane
    "TYR": "Y",  # Tyrosine
    "VAL": "V"   # Valine
    }

def wrap(seq, w):
    """
    Wrap a text sequence into multiple lines of size w
    Input: Single string
    Output: List of string of size w
    """

    return [seq[i:i + w] for i in range(0, len(seq), w)]


def parse_the_seq(file) :
    """Use to retrieve the full seq of PDB.

        Parameters
        ----------
            file : str
                Path to the PDB file.
        
        Returns
        -------
            dict_chain_seq : dict
            Dict which contains full and missing sequence of each chains.

    """

    #  Parser object.
    parser = PDBParser()

    #  PDB ID.
    pdb_id = file.split('/')[-1].split('.')[0]

    #  Structure object from Bio.PDB
    structure = parser.get_structure(pdb_id, file)

    #  Dictionnary to return.
    dict_chain_seq = dict()


    #  Iterate on model.
    for model in structure :
        
        #  Work only on the first model.
        if model.id != 0 :
            continue
        
        #  Iterate on each chains.
        for chain in model :
        
            #  Create intern dictionnary.
            dict_chain_seq[chain.id] = {'complete':"",'missing':""}

            #  Iterate on residue of query chain.
            for residue in chain.get_residues() :
            
                #  Check if it is amino acid or others.
                if residue.resname in acides_amines_3_lettres :
                    
                    #  Insert amino acid one letter code.
                    dict_chain_seq[chain.id]['complete'] +=  str(acides_amines_dict[residue.resname])

    #  Return the dictionnary.
    return dict_chain_seq






def parse_missin_seq(file,dict_chain_seq) :
    """Complete the dictionnary with the missing residues sequence.

        Parameters
        ----------
            file : str
                Path to the PDB file.
            
            dict_chain_seq : dictionnary
                Contain sequence full and missing for each chains.
        
        Returns
        -------
            dict_chain_seq : dict
                But with the missing sequence part.

    """

    #  Retrieve lines of the query PDB file.
    with open(file,'r') as f :
        lines = [i.strip() for i in f.readlines()]

    #  Dict with position of missin residue.
    dict_chain_to_missin_id = dict()

    #  Retrive for each chain the positions of missing residues
    for line in lines :
        #  Work on REMARK 465 zone.
        if line.startswith("REMARK 465     ") :
            split_remark_465 = line.strip().split()
            
            #  Retrieve chain.
            chain = split_remark_465[3]
            
            #  Retrieve residue ID.
            ID = split_remark_465[-1]

            #  Fullfill the dictionnary.
            if chain not in dict_chain_to_missin_id.keys() :
                dict_chain_to_missin_id[chain] = list()
                dict_chain_to_missin_id[chain].append(int(ID))
            else :
                dict_chain_to_missin_id[chain].append(int(ID))

    #  Complete the intial dictionnary.
    for chain in dict_chain_seq.keys() :
        
        #  Complete the sequence in the query dictionnary.
        for index in range(len(dict_chain_seq[chain]['complete'])) :
            if index+1 in dict_chain_to_missin_id[chain] :
                dict_chain_seq[chain]['missing'] += '-'
            else :
                dict_chain_seq[chain]['missing'] += dict_chain_seq[chain]['complete'][index]
    
    return dict_chain_seq


def write_ali_file(file : str, dict_chain_seq : dict) :
    """Write the ali file needed by modeller. 
    
        Parameters
        ----------
            dict_chain_seq : dict
                Contains sequence informations for each chains.
        
        Returns
        -------
            str
            Write a file
    """
    pdb_id = file.split('/')[-1].split('.')[0]

    for chain in dict_chain_seq.keys() :
        
        with open(f"{pdb_id}_chain_{chain}.ali", "w") as f :
            f.write(f">P1;{pdb_id}_{chain}_complete\n")
            f.write(f"structureX:{pdb_id}_{chain}_complete:1:.:.:.::::\n")
            
            for line in wrap(dict_chain_seq[chain]['missing']+'*',76):
                f.write(f"{line}\n")
            
            f.write(f">P1;{pdb_id}_{chain}_fill\n")
            f.write(f"sequence:::::::::\n")
            for line in wrap(dict_chain_seq[chain]['complete']+'*',76) :
                f.write(f"{line}\n")
            
    return 



def split_pdb_file_bychains(file : str, dict_chain_seq : dict) :
    """Create the separate PDB file for each chain.

        Parameters
        ----------
            file : str
                Path to the original PDB file.
            
            dict_chain_seq : dict
                Contains sequence informations for each chains.

        Returns
        -------
            str
            Create new PDB file.
    """

    pdb_id = file.split('/')[-1].split('.')[0]

    for chain in dict_chain_seq.keys() :
        pm = pymolPy3.pymolPy3(0)
        pm(f"load {file}")
        pm("remove hetatm")
        pm(f"remove not chain {chain}")
        pm(f"save {pdb_id}_{chain}_complete.pdb")
        pm("quit")
    
    return 0

def complete_modeller(dict_chain_seq : dict, file : str, chain : str) :
    """Complete the chain.

        Parameters
        ----------
            dict_chain_seq : dict
                Contain the chain and sequence informations.
            
            file : str
                Path to the original PDB file.

            chain : str
                Name of the file.
        
        Returns
        -------
            str
            New PDB file of complete chains.

    """

    #  Set parameter verbose.
    log.verbose()

    #  Set environment.
    env = environ()

    #  Directories with inputs.
    env.io.atoms_files_directory = [PATH]

    a = LoopModel(env, alnfile = file, knowns = f'4o5l_{chain}_complete', sequence = f'4o5l_{chain}_fill')

    a.starting_model = 1
    a.ending_model = 1

    a.loop.starting_model = 1
    a.loop.ending_model = 2
    a.loop.md_level = refine.fast
    
    a.make()

    print(dict_chain_seq)
    print(file)

    return 0

if __name__ == "__main__" :

    #  Retrieve full sequence and chains.
    dict_chain_seq = parse_the_seq("4o5l.pdb")

    #  Retrieve missing residues.
    dict_chain_seq = parse_missin_seq("4o5l.pdb",dict_chain_seq)

    #  Write the ali file.
    write_ali_file("4o5l.pdb",dict_chain_seq)

    #  Create pdb file for each chain
    split_pdb_file_bychains("4o5l.pdb",dict_chain_seq)

    for i in dict_chain_seq.keys() :
        if i != 'H' :
            continue
        else :
            complete_modeller(dict_chain_seq[i],f"4o5l_chain_{i}.ali",i)
