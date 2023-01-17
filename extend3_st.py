# -*- coding: utf-8 -*-
"""
Extend a file of regular 4^1 2^n designs into all non-isomorphic 4^1 2^(n+1) designs
Created on Thu Jun  3 11:17:15 2021

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""

# %% Packages
import argh
import sys
import oapackage as oa
import numpy as np

from operator import itemgetter
from itertools import chain, combinations
from math import log2
from mldoe.matrix import bmat
from time import process_time, asctime
from mldoe.enumeration import selectIsomorphismClasses
from os import devnull
from typing import List, Tuple
from tqdm import tqdm

# %% Helper functions


def get_pf(m):
    """Create triplets of pseudo-factors for m four-level factors"""
    pf = [[1, 2, 3]]
    if m == 2:
        pf.append([4, 8, 12])
    elif m > 2:
        raise ValueError("Not implemented for m > 2")
    return pf


def get_cols(ar, B, N, m, col_num):
    """ Retrieve the column numbers of an OA"""
    mat = (ar.getarray()[:, m:] * 2) - 1
    cols_mat = (B.T @ mat) // N
    cols = (col_num @ cols_mat)[0].tolist()
    return cols


def mixed_generators(k: int, m: int) -> List[Tuple]:
    """
    Create the list of generators for k basic factors and m four-level factors,
    sorted by order and by type.
    """
    # Enumerate all generators
    letters = [chr(97 + i) for i in range(26)]
    bf = letters[:k]
    g = [["".join(t) for t in combinations(bf, r + 1)] for r in range(k)]
    gen = list(chain(*g))

    # Compute original order
    og_order = [len(i) for i in gen]
    order = og_order.copy()
    gen_type = [0] * len(gen)

    # Compute the type
    pf_set = [(letters[2 * i], letters[2 * i + 1]) for i in range(m)]
    for pf in pf_set:
        for i, x in enumerate(gen):
            if any((c in x for c in pf)):
                gen_type[i] += 1
                order[i] += 1 - sum(c in x for c in pf)

    # Sort by order then type
    z = list(zip(gen, order, gen_type, og_order))
    res = sorted(z, key=itemgetter(1, 2, 3))

    return res


def generators(gen_list: List[Tuple], cols: List[int] = None) -> List[str]:
    """
    Select the possible generators for extension, using the index of the last
    added columns.
    If no columns list is provided, all generators are considered.
    """
    # Turn generator list into a number list
    num_list = [sum([2 ** (ord(i) - 97) for i in c[0]]) for c in gen_list]
    # Index of the last
    if cols is None:
        lastcol_idx = 0
    else:
        lastcol_idx = num_list.index(cols[-1]) + 1
    # Output the generators with order > 1, located after the last added column
    return [i[0] for i in gen_list[lastcol_idx:] if i[1] > 1]


def as_nums(gen_list: List[str]) -> List[int]:
    """
    Convert a list of generators from characters to numbers
    :param gen_list: list of generators represented as characters
    :return: list of generators represented as numbers
    """
    return [sum([2 ** (ord(i) - 97) for i in c]) for c in gen_list]


def find_added_cols(cols, pf, N):
    """Find the possible added columns, given a list of pre-existing columns
    using the search-table method"""
    k = int(log2(N))
    m = len(pf)
    g_list = mixed_generators(k, m)
    possible_gens = generators(g_list, cols)
    return as_nums(possible_gens)


def nauty_reduction(al: List[oa.array_link]) -> List[oa.array_link]:
    """Reduce a set of OA to only the non-isomorphic arrays, using NAUTY"""
    index, _ = selectIsomorphismClasses(al, verbose=0)
    # Select one rep. per class
    _, zz = np.unique(index, return_index=True)
    zz.sort()
    return [al[idx] for idx in list(zz)]


# %% Delcaring


@argh.arg("-p", "--progress", choices={0, 1, 2})
def main(
    N: "Run size",
    m: "Number of four-level factors",
    n: "Number of two-level factors (after extension)",
    r: "Resolution",  # noqa:F821
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
        sys.exit("ERROR: N, m, n and r must be integers")

    # Verify values of N, m, n and r
    if (N & N - 1) != 0 or N < 0:
        sys.exit("ERROR: Run size must be a positive power of two")

    if n < 1 or n > N - 1:
        sys.exit("ERROR: Number of two-level columns must be between 1 and N")

    if m < 1 or m > 3:
        sys.exit("ERROR: Number of four-level factors must be between 0 and 3")

    if r < 3 or r > 5:
        sys.exit("ERROR: Resolution must be between 3 and 5")

    # Display options used
    if log:
        log_filename = f"logs/log_N{N}_m{m}_r{r}_ST.txt"
        print(f"Logging details into {log_filename}\n")

    # Initialization
    if not silent:
        header = " {} ".format("INITIALIZING")
        print("{:#^20}".format(header))

    # Retrieve filename
    start = process_time()
    start_time = asctime()
    in_file = f"arrays/ST/{N}_{m}_{n-1}_{r}.txt"

    # Open using oapackage
    print("Reading parents arrays...", end="\r")
    sel = oa.readarrayfile(in_file)

    # Print number of parents
    if not silent:
        print(f"{len(sel)} parent designs read from {in_file}\n")

    # Create B matrix
    k = int(log2(N))
    B = bmat(k, alt_coding=True)
    B_ref = bmat(k)

    # Print start of initialization
    if not silent:
        header = " {} ".format("EXTENSION")
        print("{:#^20}".format(header))
        print(f"Extending to {n} two-level factors")

    # Find all possible added columns and create candidate set
    col_num = np.array(range(1, N))[None, :]
    pf = get_pf(m)
    candidates = []
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
        cols = get_cols(ar, B, N, m, col_num)
        added_cols = find_added_cols(cols, pf, N)
        for col in added_cols:
            des = np.concatenate((ar.getarray(), B_ref[:, col - 1][:, None]), axis=1)
            des_ar = oa.array_link(des)
            if any(des_ar.GWLP()[1:r]):
                continue
            candidates.append(des_ar)

    # Print the extension details
    if not silent:
        print(f"{len(candidates)} candidates selected with the search-table\n")

    # Start the reduction
    if not silent:
        header = " {} ".format("REDUCTION")
        print("{:#^20}".format(header))

    # Isomorphic reduction
    old_stdout = sys.stdout  # backup current stdout
    sys.stdout = open(devnull, "w")
    if len(candidates) == 1:
        lst = candidates
        # n_part = n_single = 1
    else:
        lst = nauty_reduction(candidates)
    sys.stdout = old_stdout

    # Number of representatives
    if not silent:
        print(f"{len(candidates)} candidates partitioned on WLP")
        # print(f'{n_part} sets: {n_single} singletons\n')
        # print(f'{len(lst)} representatives found')

    # Write to file
    out_file = f"arrays/ST/{N}_{m}_{n}_{r}.txt"
    print("Writing representatives arrays...", end="\r")
    oa.writearrayfile(out_file, oa.arraylist_t(lst), oa.ATEXT)
    end = process_time()
    end_time = asctime()

    # Print reduction results
    if not silent:
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
            f.write(f"Total time: {format(end-start,'.2f')} seconds\n")


# %% Dispatching
if __name__ == "__main__":
    argh.dispatch_command(main)
