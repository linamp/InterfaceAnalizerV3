from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
import xmlrpc.client as xmlrpclib
cmd = xmlrpclib.ServerProxy('http://localhost:9123')
cmd.reinitialize()
cmd.do('cd /mnt/c/Users/linam/Desktop/PyRosetta/Pymol_InterfaceAnalyzer')
cmd.do('run Scripts/InterfaceResidues.py')


def get_multimer(pdb_name, cmd, reinit=True):
    if reinit:
        cmd.reinitialize()
    cmd.load(pdb_name, 'multimer')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('not polymer')

def get_dssp(pdb_name):
    dssp_tuple = dssp_dict_from_pdb_file(pdb_name)[0]
    # Transform into dictionary with the keys in the format of the PDB file and only the values of the secondary structure and amino acid
    dssp_dict = {}
    for key in list(dssp_tuple.keys()):
        (chain_number, res_number) = key
        chain_res = str(res_number[1]) + ' ' + str(chain_number) + ' ' 
        dssp_dict[chain_res] = dssp_tuple[key][0:2]
    return dssp_dict

def get_interface(pdb_name, ch1, ch2, cmd):
    cmd.do(("interfaceResidues multimer, cA=c. {}, cB=c. {}").format(ch1, ch2))
    if cmd.do('count_atoms interface') != 0:
        try:
            cmd.do("extract interf, interface")
            cmd.save('Inputs/temp/interf.pdb', 'interf')
            print('--- SUCCESS generating interface file for chain {} and {}'.format(ch1, ch2))
            return True
        except:
            print('--- ERROR generating interface file for chain {} and {}'.format(ch1, ch2))
            try:
                cmd.reinitialize()
                get_multimer(pdb_name, cmd)
                cmd.do("extract interf, interface")
                cmd.do(("interfaceResidues multimer, cA=c. {}, cB=c. {}").format(ch1, ch2))
                cmd.save('Inputs/temp/interf.pdb', 'interf')
                print('--- SUCCESS in second attempt for chain {} and {}'.format(ch1, ch2))
                return True
            except:
                print('--- ERROR in second attempt for chain {} and {}'.format(ch1, ch2))
            return False
    else:
        return False

def get_coords(file_name='../Inputs/temp/interf.pdb', structure_name='interf'):
    parser = PDBParser()
    '''
    This function gets extracts the coordinates of the CA atoms in a PDB file,
    as well as the names of the residues they belong to.
    
    Args:
        file_name (str): path to the pdb file
        structure_name (str): name of the structure in the pdb file
        
    Returns:
        xyz (list): list of lists with the coordinates of the CA atoms
        names (list): list of lists with the names of the residues
    '''

    struct = parser.get_structure(structure_name, file_name)
    for model in struct:
        xyz = []
        names = []
        for chain in model:
            xyz_chain = []
            names_chain = []
            for residue in chain:
                for atom in residue:
                    # Get only atoms with name CA
                    if atom.get_name() == 'CA':
                        xyz_chain.append(atom.get_coord())
                        names_chain.append(str (atom.parent.id[1]) +' '+ str(atom.parent.parent.id) +' ')
            xyz.append(xyz_chain)
            names.append(names_chain)
    return xyz, names      

def distance(v1, v2):
    #get distance between two vectors in 3D space
    return np.sqrt(np.sum((v1 - v2)**2))

def get_secstruc_aminoacid (dssp_dict, pdb_resn):
    '''
    Args:
        pdb_resn: string in the format of the PDB file (chain + ' ' + chain + ' ') such as '36 A '
        dssp_dict: dictionary with the keys in the format of the PDB file and only the values of the secondary structure and amino acid
    
    Returns:
        aa: amino acid
        ss: secondary structure
    '''
    aa = dssp_dict[pdb_resn][0]
    ss = dssp_dict[pdb_resn][1]
    return aa, ss

def fastPairwiseDistance (coords_chain1, names_chain1, coords_chain2, names_chain2, dssp_dict, pdb_name, cutoff):
    '''
    This function gets the pairwise distances between the CA atoms of two chains.
    It also gets the secondary structure and amino acid of the interacting residues.

    Args:
        coords_chain1 (list): output of get_coords() for chain 1
        names_chain1 (list): output of get_coords() for chain 1
        coords_chain2 (list): output of get_coords() for chain 2
        names_chain2 (list): output of get_coords() for chain 2
        dssp_dict (dict): output of get_dssp()
        pdb_name (str): name of the PDB file
        cutoff (int): cutoff distance for the interaction between residues (in Angstroms)
    
    Returns:
        results (list): list of lists with the pairwise distances, names of the residues, secondary structure and amino acid
    '''
    results = [0] * len(coords_chain1)
    for i in range(len(coords_chain1)):
        distances = []
        for j in range(len(coords_chain2)):
            dist = distance(coords_chain1[i], coords_chain2[j])
            distances.append(dist)
        # get index of minimum value distance to detect the interaction pair
        if min(distances) > cutoff:
            results[i] = None
            continue

        min_index = distances.index(min(distances))
        aa1, ss1 = get_secstruc_aminoacid(dssp_dict, names_chain1[i])
        aa2, ss2 = get_secstruc_aminoacid(dssp_dict, names_chain2[min_index])

        results[i] = [round(distances[min_index], 3),\
            names_chain1[i], aa1, ss1,\
                names_chain2[min_index], aa2, ss2, pdb_name]
    return results

def get_interface_distances(pdb_name, cmd):
    '''Get pairwise distances between chains in a multimer

    Args:
        pdb_name (str): name of the pdb file

    Returns:
        A matrix wherein each row contains:
        - the distance between the interacting residues
        - the name of each interacting residue in the PDB File
        - their secondary structure
        - their identity (amino acid)
        - the path to the PDB file from where the data was extracted
    '''
    get_multimer(pdb_name, cmd)
    ch_interfaces = []
    dssp_dict = get_dssp('../' + pdb_name)
    chains = cmd.get_chains('multimer')
    for i in range(len(chains)):
        for j in range(i+1, len(chains)):
            get_multimer(pdb_name, cmd)
            x = get_interface(pdb_name, chains[i], chains[j], cmd=cmd)
            if x:
                try:
                    xyz, names = get_coords(file_name='../Inputs/temp/interf.pdb', structure_name='interf')
                    ch_interfaces.append(fastPairwiseDistance(xyz[0], names[0], xyz[1], names[1], dssp_dict, pdb_name, cutoff=cutoff))
                    print("--- SUCCESS reading PDB interf file of chains {} and {}".format(chains[i], chains[j]))
                except:
                    print("--- No interface found between chains {} and {}".format(chains[i], chains[j]))
    return ch_interfaces

def get_interfaces_ss(ch_interfaces, ss1='H', ss2='H'):
    '''
    This function gets the minimum distance between two chains for a given pair of secondary structures.
    It also gets the secondary structure and amino acid of the interacting residues.

    Args:
        ch_interfaces (list): output of get_interface_distances()
        ss1 (str): first secondary structure
        ss2 (str): second secondary structure

    Returns:
        min_distance_ss (list): list with the minimum distances for each interface
        stats_ss (list): Similar to ch_interfaces but with only the interfaces with minimum distance in the given secondary structures
    '''
    min_distance_ss = []
    stats_ss = []

    for i in range(len(ch_interfaces)):
        results_ss = [x for x in ch_interfaces[i] if (x[3] == ss1 and x[6] == ss2) or (x[3] == ss2 and x[6] == ss1)]
        if len(results_ss) > 0:
            stats_ss.append(results_ss)
            min_distance = min([x[0] for x in results_ss])
            min_distance_ss.append(min_distance)
    if len(min_distance_ss) == 0:
        print('No interface was found in this file for {}-{} secondary structure'.format(ss1, ss2))
    # join all the entries of stats_ss into a single list
    stats_ss = [item for sublist in stats_ss for item in sublist]
    return min_distance_ss, stats_ss

def PairwiseDistance (coords_chain1, names_chain1, coords_chain2, names_chain2, dssp_dict, pdb_name, cutoff):
    '''
    Same as fastPairwiseDistance but iterates through all the residues in BOTH interfaces chains.
    '''
    results = [0] * (len(coords_chain1)+len(coords_chain2))
    seen = []

    # iterate through all the residues in chain 1
    for i in range(len(coords_chain1)):
        distances = []
        for j in range(len(coords_chain2)):
            dist = distance(coords_chain1[i], coords_chain2[j])
            distances.append(dist)
        # get index of minimum value distance to detect the interaction pair
        if min(distances) > cutoff:
            results[i] = None
            continue
        min_index = distances.index(min(distances))
        aa1, ss1 = get_secstruc_aminoacid(dssp_dict, names_chain1[i])
        aa2, ss2 = get_secstruc_aminoacid(dssp_dict, names_chain2[min_index])
        results[i] = [round(distances[min_index], 3),\
            names_chain1[i], aa1, ss1,\
                names_chain2[min_index], aa2, ss2, pdb_name]
        seen.append([names_chain1[i], names_chain2[min_index]])
    
    # iterate through all the residues in chain 2
    for i in range(len(coords_chain2)):
        distances = []
        for j in range(len(coords_chain1)):
            dist = distance(coords_chain2[i], coords_chain1[j])
            distances.append(dist)
        # get index of minimum value distance to detect the interaction pair
        if min(distances) > cutoff:
            results[i+len(coords_chain1)] = None
            continue
        min_index = distances.index(min(distances))
        if [names_chain1[min_index], names_chain2[i]] in seen:
            results[i+len(coords_chain1)] = None
            continue
        aa1, ss1 = get_secstruc_aminoacid(dssp_dict, names_chain1[min_index])
        aa2, ss2 = get_secstruc_aminoacid(dssp_dict, names_chain2[i])
        results[i+len(coords_chain1)] = [round(distances[min_index], 3),\
            names_chain1[min_index], aa1, ss1,\
                names_chain2[i], aa2, ss2, pdb_name]
    results = [x for x in results if x is not None]
    return results

def ChainInterfaceDistances(pdb_name, ch1, ch2, cmd, cutoff=8):
    '''Get pairwise distances between two specific chains in a multimer

    Args:
        pdb_name (str): name of the pdb file
        ch1 (str): chain 1
        ch2 (str): chain 2
        cutoff (int, optional): cutoff distance. Defaults to 8.

    Returns:
        A matrix wherein each row contains:
        - the distance between the interacting residues
        - the name of each interacting residue in the PDB File
        - their secondary structure
        - their identity (amino acid)
        - the path to the PDB file from where the data was extracted
    '''
    get_multimer(pdb_name, cmd)
    ch_interfaces = []
    dssp_dict = get_dssp('../' + pdb_name)
    x = get_interface(pdb_name, ch1, ch2, cmd=cmd)
    if x:
        try:
            xyz, names = get_coords(file_name='../Inputs/temp/interf.pdb', structure_name='interf')
            ch_interfaces.append(PairwiseDistance(xyz[0], names[0], xyz[1], names[1], dssp_dict, pdb_name, cutoff=cutoff))
            print("--- SUCCESS reading PDB interf file of chains {} and {}".format(ch1, ch2))
        except:
            print("--- No interface found between chains {} and {}".format(ch1, ch2))
    return ch_interfaces[0]

def colorResidueList (chain, residue_list, cmd):
    '''Color residues in a list of residues in a chain in red'''
    for residue in residue_list:
        cmd.do('color red, chain {} and resi {}'.format(chain, residue))
    return

def totalResiduesInterface (pdb_name, ch1, ch2, cmd, min_cutoff=5, max_cutoff=10):
    '''Get the total number of residues in the interface between two chains.
        In addition, check if there are any residues whose distance between C alpha
        carbon atoms is less than a minimum cutoff distance, meaning they could
        possibly clash.

    Args:
        pdb_name (str): name of the pdb file
        ch1 (str): chain 1
        ch2 (str): chain 2
        min_cutoff (int, optional): minimum cutoff distance. Residues might clash.
        max_cutoff (int, optional): maximum cutoff distance. Defaults to 10.

    Returns:
        total number of residues in the interface, number of residues that might clash
    '''
    
    interface_results = ChainInterfaceDistances(pdb_name, ch1, ch2, cmd, cutoff=max_cutoff)
    interface_results_df = pd.DataFrame(interface_results, columns=['Distance', 'Residue1', 'AA1', 'SS1', 'Residue2', 'AA2', 'SS2', 'PDB'])
    ch1_interface = interface_results_df.Residue1.unique()
    ch2_interface = interface_results_df.Residue2.unique()
    residues_interface = len(ch1_interface) + len(ch2_interface)
    
    clash_results = ChainInterfaceDistances('Inputs/temp/multimer.pdb', ch1, ch2, cmd, cutoff=min_cutoff)
    clash_results_df = pd.DataFrame(clash_results, columns=['Distance', 'Residue1', 'AA1', 'SS1', 'Residue2', 'AA2', 'SS2', 'PDB'])
    clash_ch1_interface = clash_results_df.Residue1.unique()
    clash_ch2_interface = clash_results_df.Residue2.unique()
    clash_residues_interface = len(clash_ch1_interface) + len(clash_ch2_interface)

    return residues_interface, clash_residues_interface
