#!/usr/bin/python
"""
Enrico Ciraci 08/2022

ij_neighbors_ind - Returns the indices of neighbor cells in an n*m matrix

Input Parameters:
- n_array - numpy ndarray - input n*m matrix.
- method: neighbors index search method - def. all. - need to be improved.

Based on ixneighbors.m by Wolfgang Schwanghart
    For more info, see Flow Accumulation (upslope area) tool on MATLAB Central:
    https://www.mathworks.com/matlabcentral/fileexchange/
          14504-flow-accumulation-upslope-area?s_tid=FX_rc3_behav
--------

UPDATE HISTORY:
"""

from __future__ import print_function
import numpy as np


def ij_neighbors_ind(n_array: np.ndarray, method: str = 'getall') -> dict:
    """
    Returns the indices of neighbor cells in an n*m matrix
    :param n_array:n*m numpy array
    :param method: neighbors index search method - def. all
    :return: dict
    """
    # - Array Shape
    size = n_array.shape
    n_rows = size[0]        # - Number of Rows
    n_columns = size[1]       # - Number of Columns
    # - number of array elements
    n_rc = n_rows * n_columns

    # - matrix coordinates axes
    n_x, n_y = np.meshgrid(np.arange(n_columns), np.arange(n_rows))
    n_x_flat = n_x.flatten()
    n_y_flat = n_y.flatten()
    pt_coords = list(zip(n_y_flat, n_x_flat))

    # - Find input array NaN values
    nan_ind = np.isnan(n_array)

    # - define neighborhood size based on the selected search method
    if method == 'getall':
        nhood = 8
    else:
        # - Consider just Up/Down - Left/Right
        nhood = 4

    # - Generate index vector with the same number of elements
    # - of the input array - NOTE: initializing index values as float
    ind_range = np.arange(n_rc) + 1.
    # - Reshape the index vector to the shape of the input array
    index_arr = np.reshape(ind_range, size)
    # - Set to NaN index values of elements with NaN value in the input array.
    if index_arr[nan_ind].size > 0:
        index_arr[nan_ind] = np.nan

    # - Add Zero-Pad to index array
    index_arr_pad = np.full((size[0]+2, size[1]+2), 0, dtype=float)
    index_arr_pad[1:-1, 1:-1] = index_arr

    # - Neighbours index array
    i_nbr = np.zeros([len(ind_range.flatten()), nhood]).squeeze()
    # - Shift logical matrix index_arr_pad across the neighbors
    i_nbr[:, 0] = index_arr_pad[1:-1, 2:].flatten()    # shift to the right
    i_nbr[:, 1] = index_arr_pad[2:, 1:-1].flatten()    # shift down
    i_nbr[:, 2] = index_arr_pad[1:-1, :-2].flatten()   # shift to the left
    i_nbr[:, 3] = index_arr_pad[:-2, 1:-1].flatten()    # shift  up

    if nhood == 8:
        i_nbr[:, 4] = index_arr_pad[:-2, 2:].flatten()      # shift up and right
        i_nbr[:, 5] = index_arr_pad[:-2, :-2].flatten()     # shift up and left
        i_nbr[:, 6] = index_arr_pad[2:, 2:].flatten()    # shift down and right
        i_nbr[:, 7] = index_arr_pad[2:, :-2:].flatten()  # shift down and left

    # - Remove Pad Values
    i_nbr[i_nbr == 0] = np.nan
    # - Subtract one to convert the neighbours indexes into Python style index
    i_nbr -= 1
    index_arr -= 1

    return{'pt_coords': pt_coords, 'index_arr': index_arr,
           'index_arr_out': index_arr, 'index_nbr': i_nbr}


def main() -> None:
    # - Test Neighbors Search
    x_test = np.zeros([5, 3])
    x_test[3, 1] = np.nan
    # - compute ij_neighbours
    sample_test = 10
    nigh_ind = ij_neighbors_ind(x_test)

    # - Selected point coordinated Row, Column
    pt_index = nigh_ind['pt_coords'][sample_test]

    print(f'- Sample Point Coordinates: {pt_index}')
    print(f'- Sample Point Value: {x_test[pt_index]}')

    print(f'- Input Matrix Index Array: {pt_index}')
    print(nigh_ind['index_arr'])

    # -Search Results
    if np.isfinite(x_test[pt_index]):
        print(f'- Neighbors grid-points coordinates'
              f' (Available only for finite values):')
        print(nigh_ind['index_nbr'][sample_test])


if __name__ == '__main__':
    main()
