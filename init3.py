# -*- coding: utf-8 -*-
"""
Initialize the first batch of MLD using the oapackage and the mldoe package
Created on Thu Jun  3 11:02:13 2021

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""

# %% Packages
from typing import List
from math import log2
from itertools import chain
from mldoe.design import MLD

import argh
import os
import oapackage as oa


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


# %% Declaring
def main(
    N: "Run size",
    m: "Number of four-level factors",
    r: "Minimal resolution",
    method: "Type of method to use (one of 'ST', 'DOP', or 'MCS')",
    silent: "No std output" = False,
    log: "Create a log file with the number of designs and the total time" = False,
):
    # Check that variables are integers
    try:
        N = int(N)
        m = int(m)
        r = int(r)
    except ValueError:
        print("N, m and r must be integers")
        exit()

    # Verify values of N, m and r
    if N & N - 1 != 0 or N < 0:
        print("Run size must be a positive power of two")
        exit()
    if m < 1 or m > 3:
        print("Number of four-level factors must be between 0 and 3")
        exit()
    if r < 3 or r > 5:
        print("Resolution must be between 3 and 5")
        exit()

    # Check that the method exists
    if method not in ["ST", "DOP", "MCS"]:
        raise ValueError("Unknown method. Must be one of: 'ST', 'DOP', 'MCS'.")

    # Create arrays/ and logs/ fodlers if it doesn't exist already
    folders = ["arrays", "logs", f"arrays/{method}"]
    for f in folders:
        if not os.path.exists(f):
            os.mkdir(f)

    # Display options used
    if log:
        log_filename = f"logs/log_N{N}_m{m}_r{r}_{method}.txt"
        print(f"---Logging details into {log_filename}---")

    # Create the basic-factor
    pf = get_pf(m)

    # Create the basic design
    k = int(log2(N))
    basic_cols = [2 ** i for i in range(k) if 2 ** i not in list(chain(*pf))]
    n = len(basic_cols)
    basic_design = MLD(N, pf, basic_cols)

    # Print to txt file as OA
    ar = oa.array_link(basic_design.array)
    filename = f"arrays/{method}/{N}_{m}_{n}_{r}.txt"
    oa.writearrayfile(filename, oa.arraylist_t([ar]), oa.ATEXT)

    # Print output
    if not silent:
        header = " {} ".format("INITIALIZING")
        print("{:#^20}".format(header))
        print(f"Root array: {N} x {m + n}")
        print(f"Starting design: 4^{m} 2^{n}")
        print(f"Stored in: {filename}\n")

    # Create a log file
    if log:
        with open(log_filename, "w") as f:
            f.write(
                f"n : {n} - P : {0} - C : {1} - R: {1} - Total time: 0.00 seconds\n"
            )


# %% Dispatching
if __name__ == "__main__":
    argh.dispatch_command(main)
