# -*- coding: utf-8 -*-
"""
Extend a file of regular 4^1 2^n designs into all non-isomorphic 4^1 2^(n+1) designs
Created on Thu Jun  3 11:17:15 2021

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""

# %% Packages
from typing import List

import argh
import oapackage as oa
from math import log2
from itertools import chain
from mldoe.design import gen_len
from mldoe.matrix import bmat
import numpy as np
from time import process_time, asctime
from mldoe.enumeration import selectIsomorphismClasses
from os import devnull
import sys
from tqdm import tqdm


# %% Helper functions
def get_pf(m: int) -> List[List[int]]:
    """
    Generate the pseudo-factor triplets for m four-level factors.

    :param m: Number of four-level factors
    :type m: int
    :raises ValueError: Not implemented for m > 2
    :return: List of the pseudo-factor triplets (as list of three columns of the form a, b and ab)
    :rtype: List[List[int]]

    """
    pf = [[1, 2, 3]]
    if m == 2:
        pf.append([4, 8, 12])
    elif m > 2:
        raise ValueError("Not implemented for m > 2")
    return pf


def get_cols(ar: oa.array_link, B: np.array, m: int, N: int, col_num: np.array):
    """ Retrieve the column numbers of an OA"""
    mat = (ar.getarray()[:, m:] * 2) - 1
    cols_mat = (B.T @ mat) // N
    return (col_num @ cols_mat)[0].tolist()


def find_added_cols(
    cols: List[int], pf: List[List[int]], N: int, res: int
) -> List[int]:
    """Find all the possible added columns, given a list of pre-existing columns"""
    return [
        i
        for i in range(1, N)
        if i not in cols and i not in list(chain(*pf)) and gen_len(i, pf) >= res - 1
    ]


def is_dop(
    ar: oa.array_link,
    ref_wlp: tuple,
    n: int,
    m: int,
    cols: List[int],
    dop_wlp_cache: dict,
) -> (bool, dict, int, int):
    """Check if a design's parent has ma over all of its DOP"""
    ma = True
    n_c = 0
    tot_dop = 0
    for i in range(m, n + m):
        tot_dop += 1
        # Columns of the DOP
        dop_cols = cols.copy()
        dop_cols.pop(i - m)
        # Lookup WLP in DOP cache
        dop_cache_key = ".".join(map(str, dop_cols))
        if dop_cache_key in dop_wlp_cache.keys():
            n_c += 1
            dop_wlp = dop_wlp_cache[dop_cache_key]
        else:
            dop = ar.deleteColumn(i)
            dop_wlp = dop.GWLP()[1:]
            dop_wlp_cache[dop_cache_key] = dop_wlp
        for idx, x in enumerate(dop_wlp):
            if ref_wlp[idx] < x:
                break
            elif ref_wlp[idx] > x:
                ma = False
                break
        if not ma:
            break
    return ma, dop_wlp_cache, n_c, tot_dop


def part_isoSel(part: List[oa.array_link]) -> List[oa.array_link]:
    """Isomorphism selection on a partition of candidate designs"""
    if len(part) > 1:
        index, _ = selectIsomorphismClasses(part, verbose=0)
        # Select one rep. per class
        _, zz = np.unique(index, return_index=True)
        zz.sort()
        return [part[idx] for idx in list(zz)]
    else:
        return part


def nauty_reduction(al: List[oa.array_link]) -> List[oa.array_link]:
    """Reduce a set of OA to only the non-isomorphic arrays, using NAUTY"""
    index, _ = selectIsomorphismClasses(al, verbose=0)
    # Select one rep. per class
    _, zz = np.unique(index, return_index=True)
    zz.sort()
    return [al[idx] for idx in list(zz)]


# %% Declaring
@argh.arg("-p", "--progress", choices={0, 1, 2})
def main(
    N: "Run size",
    m: "Number of four-level factors",
    n: "Number of two-level factors (after extension)",
    r: "Minimal resolution",
    silent: "No std output" = False,
    log: "Create a log file with the number of designs and the total time" = False,
    progress: "Displays a progress bar on the processing of designs (0: None, 1: + candidates, 2: + representatives" = 2,
):

    # Check that variables are integers
    try:
        N = int(N)
        m = int(m)
        n = int(n)
        r = int(r)
    except ValueError:
        print("ERROR: N, m and r must be integers")
        exit()

    # Verify values of N, m, n and r
    if N & N - 1 != 0 or N < 0:
        print("ERROR: Run size must be a positive power of two")
        exit()
    if n < 1 or n > N - 1:
        print("ERROR: Number of two-level columns must be between 1 and N")
        exit()
    if m < 1 or m > 3:
        print("ERROR: Number of four-level factors must be between 0 and 3")
        exit()
    if r < 3 or r > 5:
        print("ERROR: Resolution must be between 3 and 5")
        exit()

    # Display options used
    if log:
        log_filename = f"logs/log_N{N}_m{m}_r{r}_DOP.txt"
        print(f"Logging details into {log_filename}\n")

    # Initialization
    if not silent:
        header = " {} ".format("INITIALIZING")
        print("{:#^20}".format(header))

    # Retrieve filename
    start = process_time()
    start_time = asctime()
    in_file = f"arrays/DOP/{N}_{m}_{n - 1}_{r}.txt"

    # Open using oapackage
    print("Reading parent arrays...", end="\r")
    sel = oa.readarrayfile(in_file)

    # Print number of parents
    if not silent:
        print(f"{len(sel)} parent designs read from {in_file}\n")

    # Create b_matrix matrix
    k = int(log2(N))
    b_matrix = bmat(k, alt_coding=True)
    b_ref = bmat(k)

    # Start of extension
    if not silent:
        header = " {} ".format("EXTENSION")
        print("{:#^20}".format(header))
        print(f"Extending to {n} two-level factors")

    # Find all possible added columns and create candidate set
    col_num = np.arange(1, N)[None, :]
    pf = get_pf(m)
    overall_candidates = 0
    candidates = []
    wlp_cache = dict()  # Cache for the WLP values
    # Log of the number of cached DOP values
    n_cached = tot = 0
    if progress > 0:
        pbar = tqdm(
            sel,
            leave=True,
            desc="Parent designs processed",
            total=len(sel),
            ascii=True,
            unit="des.",
        )
    else:
        pbar = sel
    for ar in pbar:
        # Compute WLP of the parent
        ref_wlp = ar.GWLP()[1:]
        # Find the columns of the parent
        cols = get_cols(ar, b_matrix, m, N, col_num)
        # Add parent's WLP to cache
        parent_wlp_key = ".".join(map(str, cols))
        wlp_cache[parent_wlp_key] = ref_wlp
        # Enumerate possible added columns
        added_cols = find_added_cols(cols, pf, N, r)
        # pool = Pool()
        for col in added_cols:
            overall_candidates += 1
            # Create candidate design as np.array then oa.array_link
            des = np.concatenate((ar.getarray(), b_ref[:, col - 1][:, None]), axis=1)
            des_ar = oa.array_link(des)
            # Check for resolution of the candidate
            if any(des_ar.GWLP()[1:r]):
                continue
            # Compute DOP among MA and increase the DOP cache
            ma, wlp_cache, n_cached_temp, tot_temp = is_dop(
                des_ar, ref_wlp, n, m, cols + [col], wlp_cache
            )
            n_cached += n_cached_temp
            tot += tot_temp
            # if candidate has MA among all DOP, add to candidate list
            if ma:
                candidates.append(des_ar)

    # Print the extension details
    if not silent:
        print(f"{overall_candidates} potential candidates considered")
        print(f"{len(candidates)} selected using DOP selection and resolution check")
        if tot != 0:
            perc = round(n_cached / tot, 2) * 100
        else:
            perc = 0
        print(f"\t{n_cached} on {tot} ({perc}%) of WLP values of DOP cached\n")

    # Isomorphic reduction
    old_stdout = sys.stdout  # backup current stdout
    sys.stdout = open(devnull, "w")
    # Check for number of candidates
    if len(candidates) == 1:
        lst = candidates
        # n_part = n_single = 1
    else:
        lst = nauty_reduction(candidates)
    sys.stdout = old_stdout

    # Details on partitioning
    if not silent:
        header = " {} ".format("PARTITION")
        print("{:#^20}".format(header))
        print(f"{len(candidates)} candidates partitioned on WLP")
        # print(f'{n_part} sets: {n_single} singletons\n')

    # Write to file
    out_file = f"arrays/DOP/{N}_{m}_{n}_{r}.txt"
    print("Writing representatives arrays...", end="\r")
    oa.writearrayfile(out_file, oa.arraylist_t(lst), oa.ATEXT)
    end = process_time()
    end_time = asctime()
    # Print reduction results
    if not silent:
        header = " {} ".format("REDUCTION")
        print("{:#^20}".format(header))
        print(f"{len(lst)} representatives found")
        print(f"Printed to {out_file}\n")

    # Print time
    if not silent:
        header = " {} ".format("TIME")
        print("{:#^20}".format(header))
        print(f"Started at {start_time}")
        print(f"Ended at {end_time}")
        print(f"Total elapsed time: {format(end - start, '.2f')} seconds")

    # Create a log file
    if log:
        with open(log_filename, "a") as f:
            f.write(
                f"n : {n} - P : {len(sel)} - C : {len(candidates)} - R: {len(lst)} - "
            )
            f.write(f"Total time: {format(end - start, '.2f')} seconds\n")


# %% Dispatching
if __name__ == "__main__":
    argh.dispatch_command(main)
