import logging
import os
import pathlib
import pickle
import sys
from typing import Any, Optional, Union

import cobra 
import pandas as pd

class HiddenPrints:
    '''Supress package print messages.'''
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def read_metabolic_model(file_name: str) -> cobra.Model:
    """Read a metabolic model from sbml format with .xml extension

    Parameters
    ----------
    file_name : str
        'full/path/to/metabolic_model.xml'

    Returns
    -------
    m_model : cobra.Model
        cobrapy metabolic model
    """
    if os.path.splitext(file_name)[1] != '.xml':
        raise ValueError('Specified file must be an sbml model with extentsion ".xml"')

    logging.basicConfig()
    logger = logging.getLogger(cobra.__name__)
    logger.setLevel(logging.CRITICAL)

    with HiddenPrints():
        m_model = cobra.io.read_sbml_model(file_name)

    return m_model

def load_metabolic_model(model_file: Union[str, cobra.Model]) -> cobra.Model:
    """Reads in metabolic model, checking for appropriate type

    Parameters
    ----------
    model_file : Union[str, cobra.Model]
        cobrapy model or 'full/path/to/metabolic_model.xml'

    Returns
    -------
    m_model : cobra.Model
        cobrapy metabolic model
    """
    if isinstance(model_file, str):
        return read_metabolic_model(model_file)
    if isinstance(model_file, cobra.Model):
        return model_file
    else:
        raise TypeError('Model arg must either by a cobrapy model or specify a path to a sbml file of a cobrapy model')

def read_psim(psim_file: str, h5_key: Optional[str] = None) -> pd.DataFrame:
    """Read the PSIM from .csv or .h5 format

    Parameters
    ----------
    psim_file : str
        'full/path/to/psim.csv or psim.h5'
    h5_key : str, optional
        key associated with h5 file, by default None

    Returns
    -------
    psim : pd.DataFrame
        see PSIM_README.md
    """
    _, file_extension = os.path.splitext(psim_file)
    if file_extension == '.h5':
        if h5_key is None:
            psim = pd.read_hdf(psim_file)
        else:
            psim = pd.read_hdf(psim_file, key = h5_key) # key = 'corrected'
    elif file_extension == '.csv':
        psim = pd.read_csv(psim_file, index_col=0)
    else:
        raise TypeError('PSIM must be in .csv or .h5 format')

    if 'SP' in psim.columns:
        psim['SP'] = psim['SP'].apply(lambda x: bool(x))

    return psim


def load_psim(psim_file: Union[str, pd.DataFrame], **kwargs) -> pd.DataFrame: 
    """Reads in PSIM, checking for appropriate type

    Parameters
    ----------
    psim_file : Union[str, pd.DataFrame]
        the PSIM dataframe or 'full/path/to/psim.csv or psim.h5'

    Returns
    -------
    pd.DataFrame
        the PSIM dataframe
    """

    if isinstance(psim_file, str):
        return read_psim(psim_file, **kwargs)
    if isinstance(psim_file, pd.DataFrame):
        if 'SP' in psim_file.columns and psim_file.SP.dtype != bool:
            psim_file['SP'] = psim_file['SP'].apply(lambda x: bool(x))
        return psim_file
    else:
        raise TypeError('The specified psim_file must be a path to the dataframe of the pandas DataFrame')

def read_pickled_me_model(file_name: str):
    """Loads a pickled me_model. Saved from me_model.pickle

    Parameters
    ----------
    file_name : str
        'full/path/to/me_model.pickle'

    Returns
    -------
    ME_Model
        ME model object
    """

    with open(file_name, 'rb') as handle:
        me_model = pickle.load(handle)
    me_model.correct_object_tracking()  # lost in pickling/loadings
    return me_model


def write_metabolic_model(m_model: cobra.Model, file_name: str, **kwargs) -> None: 
    """Write a metabolic model to smbl format with .xml extension

    Parameters
    ----------
    m_model : cobra.Model
        cobrapy metabolic model
    file_name : str
        'full/path/to/metabolic_model.xml'
    **kwargs
        additional parameters to cobra.io.write_sbml_model
    """
    cobra.io.write_sbml_model(cobra_model=m_model, filename=file_name, **kwargs)


def write_pickled_object(object: Any, file_name: str) -> None:
    """Save an object as a pickled file

    Parameters
    ----------
    object : Any
        object to save
    file_name : str
        'full/path/to/file.pickle'
    """
    if '.' in file_name:
        p = pathlib.Path(file_name)
        extensions = "".join(p.suffixes)
        file_name = str(p).replace(extensions, '.pickle')
    else:
        file_name = file_name + '.pickle'

    with open(file_name, 'wb') as handle:
        pickle.dump(object, handle)
