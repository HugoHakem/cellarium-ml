import os
import sys
from pathlib import Path
from zipfile import ZipFile

import requests
from tqdm import tqdm

# Code Credits: https://github.com/snap-stanford/GEARS/blob/master/gears/utils.py

def print_sys(s: str):
    """
    system pring in error output
xR
    Args:
        s (str): the string to print
    """
    print(s, flush=True, file=sys.stderr)

def dataverse_download(url: str, save_path: Path):
    """
    Dataverse download helper with progress bar

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
    """

    if os.path.exists(save_path):
        print_sys('Found local copy...')
    else:
        print_sys("Downloading...")
        response = requests.get(url, stream=True)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
        with open(save_path, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()

def zip_data_download_wrapper(url: str, save_path: Path, data_path: Path):
    """
    Wrapper for zip file download

    Args:
        url (str): the url of the dataset
        save_path (str): the path where the file is donwloaded
        data_path (str): the path to save the extracted dataset
    """

    if os.path.exists(data_path.joinpath(save_path.stem)):
        print_sys('Found local copy...')
    else:
        dataverse_download(url, save_path.with_suffix('.zip'))
        print_sys('Extracting zip file...')
        with ZipFile(save_path.with_suffix('.zip'), 'r') as zip:
            zip.extractall(path = data_path)
        print_sys("Done!")


# def load(self, data_name = None, data_path = None):
#     """
#     Load existing dataloader
#     Use data_name for loading 'norman', 'adamson', 'dixit' datasets
#     For other datasets use data_path

#     Parameters
#     ----------
#     data_name: str
#         Name of dataset
#     data_path: str
#         Path to dataset

#     Returns
#     -------
#     None

#     """

#     if data_name in ['norman', 'adamson', 'dixit',
#                         'replogle_k562_essential',
#                         'replogle_rpe1_essential']:
#         ## load from harvard dataverse
#         if data_name == 'norman':
#             url = 'https://dataverse.harvard.edu/api/access/datafile/6154020'
#         elif data_name == 'adamson':
#             url = 'https://dataverse.harvard.edu/api/access/datafile/6154417'
#         elif data_name == 'dixit':
#             url = 'https://dataverse.harvard.edu/api/access/datafile/6154416'
#         elif data_name == 'replogle_k562_essential':
#             ## Note: This is not the complete dataset and has been filtered
#             url = 'https://dataverse.harvard.edu/api/access/datafile/7458695'
#         elif data_name == 'replogle_rpe1_essential':
#             ## Note: This is not the complete dataset and has been filtered
#             url = 'https://dataverse.harvard.edu/api/access/datafile/7458694'
#         data_path = os.path.join(self.data_path, data_name)
#         zip_data_download_wrapper(url, data_path, self.data_path)
#         self.dataset_name = data_path.split('/')[-1]
#         self.dataset_path = data_path
#         adata_path = os.path.join(data_path, 'perturb_processed.h5ad')
#         self.adata = sc.read_h5ad(adata_path)

#     elif os.path.exists(data_path):
#         adata_path = os.path.join(data_path, 'perturb_processed.h5ad')
#         self.adata = sc.read_h5ad(adata_path)
#         self.dataset_name = data_path.split('/')[-1]
#         self.dataset_path = data_path
#     else:
#         raise ValueError("data attribute is either norman, adamson, dixit "
#                             "replogle_k562 or replogle_rpe1 "
#                             "or a path to an h5ad file")

#     self.set_pert_genes()
#     print_sys('These perturbations are not in the GO graph and their '
#                 'perturbation can thus not be predicted')
#     not_in_go_pert = np.array(self.adata.obs[
#                                 self.adata.obs.condition.apply(
#                                 lambda x:not filter_pert_in_go(x,
#                                     self.pert_names))].condition.unique())
#     print_sys(not_in_go_pert)

#     filter_go = self.adata.obs[self.adata.obs.condition.apply(
#                             lambda x: filter_pert_in_go(x, self.pert_names))]
#     self.adata = self.adata[filter_go.index.values, :]
#     pyg_path = os.path.join(data_path, 'data_pyg')
#     if not os.path.exists(pyg_path):
#         os.mkdir(pyg_path)
#     dataset_fname = os.path.join(pyg_path, 'cell_graphs.pkl')

#     if os.path.isfile(dataset_fname):
#         print_sys("Local copy of pyg dataset is detected. Loading...")
#         self.dataset_processed = pickle.load(open(dataset_fname, "rb"))
#         print_sys("Done!")
#     else:
#         self.ctrl_adata = self.adata[self.adata.obs['condition'] == 'ctrl']
#         self.gene_names = self.adata.var.gene_name


#         print_sys("Creating pyg object for each cell in the data...")
#         self.create_dataset_file()
#         print_sys("Saving new dataset pyg object at " + dataset_fname)
#         pickle.dump(self.dataset_processed, open(dataset_fname, "wb"))
#         print_sys("Done!")
