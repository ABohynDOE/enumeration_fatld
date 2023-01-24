# -*- coding: utf-8 -*-
"""
Extend a file of regular 4^1 2^n designs into all non-isomorphic 4^1 2^(n+1) designs
Created on Wed Jun  2 11:52:11 2021

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""

from math import log2
from time import asctime, process_time

# %% Packages
import argh
import numpy as np
import oapackage as oa
from mldoe.matrix import bmat
from tqdm import tqdm

# %% Declaring


@argh.arg("-p", "--progress", choices={0, 1, 2})
def main(
    N: "Run size", # noqa
    m: "Number of four-level factors", # noqa
    n: "Number of two-level factors (after extension)", # noqa
    r: "Resolution",  # noqa
    silent: "No std output" = False, # noqa
    log: "Create a log file with the number of designs and the total time" = False, # noqa
    progress: "Displays a progress bar on the processing of designs (0: None, 1: + candidates, 2: + representatives" = 2, # noqa
):

    # Check that variables are integers
    try:
        N = int(N)
        m = int(m)
        n = int(n)
        r = int(r)
    except ValueError:
        print("N, m, n and r must be integers")

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

    # Display logging information
    if log:
        log_filename = f"logs/log_N{N}_m{m}_r{r}_MCS.txt"
        print(f"Logging details into {log_filename}\n")

    # Initialization
    if not silent:
        header = " {} ".format("INITIALIZING")
        print("{:#^20}".format(header))

    # Retrieve filename
    start = process_time()
    start_time = asctime()
    in_file = f"arrays/MCS/{N}_{m}_{n-1}_{r}.txt"

    # Open using oapackage
    print("Reading parents arrays...", end="\r")
    sel = oa.readarrayfile(in_file)

    # Print number of parents
    if not silent:
        print(f"{len(sel)} parent designs read from {in_file}\n")

    # Start of extension
    if not silent:
        header = " {} ".format("EXTENSION")
        print("{:#^20}".format(header))
        print(f"Extending to {n} two-level factors")

    # Create arrayclass
    s_lst = [4] * m + [2] * n
    arrayclass = oa.arraydata_t(s_lst, N, r - 1, m + n)

    # Extend array list
    if not silent:
        print("Extending the set of candidates...", end="\r")
    lst = oa.extend_arraylist(sel, arrayclass)
    if not silent:
        print(
            f"{len(lst)} candidates created through extension with the MCS algorithm\n"
        )

    # Start of reduction
    if not silent:
        header = " {} ".format("REDUCTION")
        print("{:#^20}".format(header))

    # Build B matrix
    k = int(log2(N))
    B = ((bmat(k) * 2) - 1).T

    # Keep only regular arrays
    reg_ar = []
    if progress > 0:
        pbar = tqdm(
            lst,
            leave=True,
            desc="Candidate designs processed",
            total=len(lst),
            ascii=True,
            unit="des.",
        )
    else:
        pbar = lst
    for ar in pbar:
        # Regularity check
        mat = (ar.getarray()[:, m:] * 2) - 1
        cols = (B @ mat // N).sum(1)
        if np.array_equal(cols, cols.astype(bool)):
            reg_ar.append(ar)

    # Print filtering results
    if not silent:
        print(f"{len(reg_ar)} regular arrays (representatives) found")

    # Write to file
    out_file = f"arrays/MCS/{N}_{m}_{n}_{r}.txt"
    print("Writing regular arrays...", end="\r")
    oa.writearrayfile(out_file, oa.arraylist_t(reg_ar), oa.ATEXT)
    end = process_time()
    end_time = asctime()

    # Print reduction result
    if not silent:
        print(f"Printed to {out_file}\n")

    # Print time
    if not silent:
        header = " {} ".format("TIME")
        print("{:#^20}".format(header))
        print(f"Started at {start_time}")
        print(f"Ended at {end_time}")
        print(f"Total elapsed time: {format(end-start,'.2f')} seconds")

    # Create a log file
    if log:
        with open(log_filename, "a") as f:
            f.write(f"n : {n} - P : {len(sel)} - C : {len(lst)} - R: {len(reg_ar)} - ")
            f.write(f"Total time: {format(end-start,'.2f')} seconds\n")


# %% Dispatching
if __name__ == "__main__":
    argh.dispatch_command(main)
    argh.dispatch_command(main)
