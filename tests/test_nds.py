#!/usr/bin/env python3

import os
import pandas as pd
import re
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting

file_name = os.path.abspath('../resources/csv_files/population_test.csv')

df = pd.read_csv(file_name, header=0)

cnames = df.columns
constrs = [f for f in cnames if re.match('^constr_', f)]

for constr in constrs:
    df = df[df[constr] <= 0]

F = df.drop(constrs, axis=1).values

nds = NonDominatedSorting()
front = nds.do(F, only_non_dominated_front=True)
pf = F[front, :]
