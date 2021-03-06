#!/usr/bin/env python3

import swat_utilities.prepare_ops as prep

inputDIR = '/home/sebwaymaster/Desktop/Saginaw_Project/Saginaw_Model_WQ'
outputDIR = '/home/sebwaymaster/Desktop/Saginaw_Project/Saginaw_Model_WQ'
start_date = {'MONTH': 1, 'DAY': 1, 'IYEAR': 1998}

prep.create_ops(inputDIR, outputDIR, start_date)
prep.link_sdr(inputDIR, outputDIR)