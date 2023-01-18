# -*- coding: utf-8 -*-
"""
Find a reduced bound for all 4^m 2^n designs with n<20,
using the DOP of the type-m MA design with n=20, and going
down to n < 20.
If there are several MA designs, then the set of all designs is considered to
find the worst DOP.

Created on Tue Mar  1 10:50:12 2022

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""
from itertools import chain, combinations
from typing import List, Tuple

# % Packages
import csw93
import numpy as np
import oapackage as oa
from mldoe.enumeration import selectIsomorphismClasses
from tqdm import tqdm

# % Declaration
m = 2
run_size = 64


def twlp(mixed_array: oa.array_link):
    """ Compute the aligned type-specific word length pattern"""
    dist_matrix = oa.distance_distribution_mixed(mixed_array, 0)
    runsize = mixed_array.n_rows
    array_class = oa.arraylink2arraydata(mixed_array)
    b_prime_array = oa.macwilliams_transform_mixed(
        dist_matrix,
        runsize,
        factor_levels_for_groups=array_class.factor_levels_column_groups(),
        verbose=0,
    )
    b_prime_mat = np.array(b_prime_array).astype(int)
    r, c = b_prime_mat.shape
    w_vec = np.zeros((r, r + c + -1), dtype=int)
    for index in range(r + c - 1):
        if index < r and index < c:
            index_combinations = [(i, index - i) for i in range(index + 1)]
        else:
            index_combinations = [
                (i, index - i)
                for i in range(index + 1)
                if ((i < r) and (index - i) < c)
            ]
        for coord in index_combinations:
            w_vec[coord[0], index] = b_prime_mat[coord[0], coord[1]]
    return w_vec


def make4lvl(mat, pairs, m_factor) -> np.array:
    """Generate m four-level factors from them pairs of two-level factors"""
    n_rows, n_cols = mat.shape
    # Mask for the two-level part
    mask = np.ones(n_cols, dtype=bool)
    for i in chain(*pairs):
        mask[i] = False
    # Isolate two-level part
    two_part = mat[:, mask]
    # Create the four-level columns
    four_part = np.zeros((n_rows, m_factor))
    for i in range(m_factor):
        c1, c2 = pairs[i]
        four_part[:, i] = 2 * mat[:, c1] + mat[:, c2]
    return np.concatenate((four_part, two_part), 1).astype(int)


def tfi_matrix(ar: oa.array_link, pseudo_factors: List[Tuple[int]]) -> np.array:
    """Create a two-factor interaction matrix of the array, with unfolded four-level factors"""
    # Get the matrix of the design
    matrix = ar.getarray()
    pf_list = list(chain(*pseudo_factors))
    pf_matrix = matrix[:, pf_list]
    r, c = ar.shape
    r, c = pf_matrix.shape
    size = int(c * (c - 1) * (1 / 2) + c * (c - 1) * (c - 2) * (1 / 6))
    interaction_matrix = np.zeros((r, size), dtype=int)
    # k-factor interactions
    i = 0
    fact_int = []
    for k in (2, 3):
        for x in combinations(range(c), k):
            interaction_matrix[:, i] = np.sum(matrix[:, x], axis=1)
            fact_int.append(",".join(map(str, x)))
            i += 1
    result_matrix = np.matmul(pf_matrix.T * 2 - 1, interaction_matrix * 2 - 1)
    return result_matrix


# All triplets of pseudo-factors pairs
# Create all pairs
pairs_ntuples = []
indices = list(range(2 * m))
if m == 1:
    pairs_ntuples = list(combinations(indices, 2))
else:
    for pair1 in list(combinations(indices, 2)):
        rem_indices1 = indices.copy()
        for x in pair1:
            rem_indices1.remove(x)
        if m == 2:
            pair2 = tuple(rem_indices1)
            pairs_ntuples.append((pair1, pair2))
        else:
            for pair2 in combinations(rem_indices1, 2):
                rem_indices2 = rem_indices1.copy()
                for y in pair2:
                    rem_indices2.remove(y)
                pair3 = tuple(rem_indices2)
                pairs_ntuples.append((pair1, pair2, pair3))
# Test for unicity
ntuples_keys = []
if m == 1:
    unique_ntuples = [pairs_ntuples]
else:
    for ntuple in pairs_ntuples:
        str_keys = ["".join(map(str, sorted(pair))) for pair in ntuple]
        ntuples_keys.append(";".join(sorted(str_keys)))
    _, unique_keys_id = np.unique(ntuples_keys, return_index=True)
    # Select unique triplets
    unique_ntuples = [pairs_ntuples[i] for i in unique_keys_id]


# Number of two-level designs to consider, based on CSW93 and bounded A3 value found
# for the moment
if m == 3:
    max_rank = 4
    target_A31 = 30
    max_n = 9
else:
    max_rank = 5
    if m == 1:
        target_A31 = 5
        max_n = 16
    elif m == 2:
        target_A31 = 16
        max_n = 12
    else:
        raise ValueError

n_factors = 20 + 2 * m
n_added = n_factors - int(np.log2(run_size))

# % Activation
if __name__ == "__main__":
    # Get all the two-level design matrices
    designs = [
        csw93.get_design(64, f"{n_factors}-{n_added}.{i + 1}") for i in range(max_rank)
    ]

    # All combinations of 2m factors (as pseudo-factor pairs)
    factor_combinations = []
    for comb in combinations(range(n_factors), 2 * m):
        for ntuple in unique_ntuples:
            factor_combinations.append(
                [(comb[pair[0]], comb[pair[1]]) for pair in ntuple]
            )

    # Only keep the target designs
    MA_designs = []
    for matrix in designs:
        for f in tqdm(factor_combinations, desc="Factor combinations"):
            # Create four-level factors out of pairs of two-level factors
            fat_mat = oa.array_link(make4lvl(matrix, f, m))
            # Reject if resolution is lower than 3 or if target WLP is already too high
            wlp = fat_mat.GWLP()
            if wlp[2] != 0 or wlp[3] > target_A31:
                continue
            else:
                # Compute type specific wlp to extract number of words of length 3
                # and type 1
                t_wlp = twlp(fat_mat)
                A3 = t_wlp[:, 3].tolist()
                if m == 1:
                    if A3[0] == 0 and A3[1] == target_A31:
                        MA_designs.append(fat_mat)
                elif m == 2:
                    if A3[0] == 0 and A3[1] == target_A31 and A3[2] == 0:
                        MA_designs.append(fat_mat)
                else:
                    if A3[-1] == 0:
                        MA_designs.append(fat_mat)

    # Remove isomorphic designs from the list
    MA_designs_classes, _ = selectIsomorphismClasses(MA_designs)
    _, MA_designs_idx = np.unique(MA_designs_classes, return_index=True)
    unique_MA_designs = [MA_designs[idx] for idx in MA_designs_idx]

    # Find the type-m MA design(s)
    # Matrix of all the TWLP
    twlp_list = [twlp(x)[:, 3].tolist() for x in unique_MA_designs]
    twlp_mat = np.zeros((len(twlp_list), m + 1), dtype=int)
    for i, x in enumerate(twlp_list):
        twlp_mat[i, :] = x
    # Lexicographic sorting of TWLP (from m to 0)
    twlp_sort_index = np.lexsort(twlp_mat.T)
    MA_twlp = twlp_list[twlp_sort_index[0]]
    # Extract designs that have MA twlp
    MA_designs = []
    for i, x in enumerate(unique_MA_designs):
        if twlp_list[i] == MA_twlp:
            MA_designs.append(x)

    # Loop until the bound is not needed anymore
    result = {"20": MA_twlp}

    for n in range(19, max_n - 1, -1):
        dop_designs = []
        # Create all DOP and sort them by WLP (A3 values only)
        for design in MA_designs:
            array = design.getarray()
            for i in range(array.shape[1] - m):
                # Remove only two-level factors
                dop_array = np.delete(array, i + m, axis=1)
                dop_designs.append(dop_array)
        # Find TWLP of DOP designs
        dop_twlp_list = [twlp(oa.array_link(x))[:, 3].tolist()
                         for x in dop_designs]
        dop_twlp_mat = np.array(dop_twlp_list)
        # Lexicographic sorting of TWLP (from m to 0)
        dop_twlp_sort_index = np.lexsort(dop_twlp_mat.T)
        MA_dop_twlp = dop_twlp_list[dop_twlp_sort_index[-1]]

        # Save DOP MA TWLP to dictionnary
        # Extract designs that have MA twlp
        dop_MA_designs = []
        for i, x in enumerate(dop_designs):
            if dop_twlp_list[i] == MA_dop_twlp:
                dop_MA_designs.append(oa.array_link(x))

        # Remove isomorphic designs from the list if more than one
        if len(dop_MA_designs) > 1:
            dop_MA_designs_classes, _ = selectIsomorphismClasses(
                dop_MA_designs)
            _, dop_MA_designs_idx = np.unique(
                dop_MA_designs_classes, return_index=True)
            unique_dop_MA_designs = [dop_MA_designs[idx]
                                     for idx in dop_MA_designs_idx]
        else:
            unique_dop_MA_designs = dop_MA_designs

        # Extract the designs
        result[f"{n}"] = MA_dop_twlp
        print(f"{n} : {MA_dop_twlp}")
        # Save them as unique_MA_designs
        MA_designs = unique_dop_MA_designs

    # Write the results to a file
    with open(f"A3_bounds/reduced_bound_m{m}.txt", "w") as f:
        for k, v in result.items():
            f.write(f"{k} : {v}\n")
