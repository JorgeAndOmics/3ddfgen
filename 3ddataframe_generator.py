import WCN
import glob
import os
import numpy as np
import warnings
import statistics
import pandas as pd

def 3df_generator():
    # .pdb files are usually full of discontinuities in their protein chains, which raise warnings. As there is
    # currently no known method to block them specifically, warnings are filtered off the console.
    warnings.filterwarnings('ignore')

    def input_filepath():
        '''
        :return: str: Input filepath. Detect invalid inputs and ask for path again.
        '''
        filepath: str = input('Paste input directory') or '/home/lympha/RAS/isoraspdb/'
        try:
            filepath = str(filepath)
            return str(filepath)
        except:
            print('Wrong Input. Please, Try Again:')
            return input_filepath()

    def output_filepath():
        '''
        :return: str: Input filepath. Detect invalid inputs and ask for path again.
        '''
        filepath: str = input('Paste input directory') or '/home/lympha/RAS/csv'
        try:
            filepath = str(filepath)
            return str(filepath)
        except:
            print('Wrong Input. Please, Try Again:')
            return output_filepath()

    # List initialization
    pdb_names: list = []
    pdb_classification: list = []
    bfactor_wcn: list = []
    c_alpha_wcn: list = []
    all_atom_wcn: list = []
    average_pr_wcn: list = []

    # numpy can´t add values to an array directly. It generates a new array in memory each time an element is added
    # through the .append method. Two solutions exist for this: Preallocate memory before computation or generate a
    # list and subsequently transform in into an array

    # I get the files through a loop, parse them with the WCN object initializer and compute them through
    # the calculator created beforehand
    for pdb in glob.glob(os.path.join(input_filepath(), '*.pdb')): # Initialize the object and return method outputs
        loc: object = WCN.WCNObject(pdb)
        bf: np.ndarray = loc.getBfactors()
        ca: np.ndarray = loc.calculateCAlpha()
        aa: np.ndarray = loc.calculateAllAtom()
        ar: np.ndarray = loc.calculateAveragePerResidue()
        pdb_file_name: str = pdb.split('/')[-1] # Split the filepath name and return protein name string
        pdb_protein: str = pdb_file_name.split('.')[0]
        pdb_names.append(pdb_protein)
        if '_a' in pdb_protein:# Search protein name and return 'active' or 'inactive' string
            pdb_classification.append('Active')
        elif '_i' in pdb_protein:
            pdb_classification.append('Inactive')
        bfactor_wcn.append(bf)
        c_alpha_wcn.append(ca)
        all_atom_wcn.append(aa)
        average_pr_wcn.append(ar)

    # Calculations to dictionary
    ras_data = {'protein':pdb_names,
                'classification': pdb_classification,
                'b_factor': bfactor_wcn,
                'c_alpha': c_alpha_wcn,
                'all_atom': all_atom_wcn,
                'per_residue': average_pr_wcn,
                }

    # Generate dataframe and exports to file. ¡¡¡.csv saves everything as strings!!!
    ras_df: pd.DataFrame = pd.DataFrame(ras_data)
    ras_df.to_pickle(f'{output_filepath()}/raswholematrix.pkl')
    print(ras_df.info())


