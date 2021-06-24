#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This is a code snipped to match areas/subnetworks/color files with matrices which had some regions removed

"""

__author__ = "Eduarda"
__contact__ = "e.centeno@amsterdamumc.nl"
__date__ = "2021/6/10"
__status__ = "Finished"


####################
# Libraries        #
####################

# Third party imports 
import pandas as pd


regs_removed = [1, 2, 5, 9] # your removed ROIs as indexes

def corr_regs_remvd(idx_regs_rmvd=regs_removed, target='all'):
    """ Function that removes the indices(ROIs) that were removed from your original matrix.
    
    Parameters
    ----------
    idx_regs: list
        ROI indices to be discarded
        
    target: str
        The file in which regions should be removed. Default = 'all'
        If 'all', it will work on the variables available in the notebook (lineList, sublist, colorlist, colornumbs)
        
    
    Returns
    ------
    If target=='all':   lineList_rmvd: list
                        sublist_rmvd: list
                        colorlist_rmvd: list
                        colornumbs_rmvd: np.array
                        
    If target=='somefile':  file_rmvd: DataFrame
    
    
    """
    
    if target=='all':
        lineList_rmvd = list(pd.DataFrame(lineList).drop(regs_removed, axis=0)[0])
        sublist_rmvd = list(pd.DataFrame(sublist).drop(regs_removed, axis=0)[0])
        colorlist_rmvd = list(pd.DataFrame(colorlist).drop(regs_removed, axis=0)[0])
        colornumbs_rmvd = pd.DataFrame(colornumbs).drop(regs_removed, axis=0)[0].values
        
        return lineList_rmvd, sublist_rmvd, colorlist_rmvd, colornumbs_rmvd
    
    else:
        file_rmvd = pd.read_csv(target, header=None).drop(regs_removed, axis=0) # this might need to be adapted to your file
        return file_rmvd
