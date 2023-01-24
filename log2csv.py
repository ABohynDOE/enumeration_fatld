# -*- coding: utf-8 -*-
"""
Retrieve all log files for each test case and gather the results.

Created on Tue Jun  15 14:38:00 2021

@author: Alexandre Bohyn - alexandre dot bohyn [at] kuleuven dot be
"""
import glob
import re
from collections import defaultdict

import pandas as pd


def get_log_info(log_filename: str) -> (int, int, int, str):
    """Function to extract info from log file name"""
    n_val = int(re.search(r"N(\d+)", log_filename).group(1))
    m_val = int(re.search(r"m(\d+)", log_filename).group(1))
    r_val = int(re.search(r"r(\d+)", log_filename).group(1))
    method_str = re.search(r"_([A-Z]+)\.txt", log_filename).group(1)
    return n_val, m_val, r_val, method_str


if __name__ == "__main__":
    # List the dirs of all the methods
    log_files = glob.glob("logs/log_*.txt")

    col_names = ("n", "parents", "candidates", "representatives", "time")
    df_dict = defaultdict(pd.DataFrame)
    for log in log_files:
        # Metadata
        N, m, r, method = get_log_info(log)

        # Create dataframe from the log file
        df = pd.DataFrame()
        with open(log, "r") as f:
            for line in f.readlines():
                temp_data = dict()
                match = re.finditer(r":\s+(\d+\.*\d*)", line)
                for idx, val in enumerate(match):
                    temp_data[col_names[idx]] = [float(val.group(1))]
                    temp_data["method"] = method
                    temp_data_dict = pd.DataFrame.from_dict(temp_data)

                df = pd.concat([df, temp_data_dict], ignore_index=True)

        # Append df to dictionnary of df
        key = f"N{N}_m{m}_r{r}"
        df_dict[key] = pd.concat([df_dict[key], df], ignore_index=True)

    # Make each entry of the dict of df wide
    with pd.ExcelWriter(
        "results/global_results.xlsx"
    ) as writer:  # pylint: disable=abstract-class-instantiated
        for key in df_dict.keys():
            or_df = df_dict[key]
            or_df.to_excel(f"results/{key}.xlsx", index=False)
            wide_df = or_df.pivot_table(
                index="n",
                values=["parents", "candidates", "representatives", "time"],
                columns="method",
            )
            # Write it to the excel file
            wide_df.to_excel(writer, sheet_name=key)
            wide_df.to_excel(writer, sheet_name=key)
