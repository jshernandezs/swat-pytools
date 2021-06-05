#!/usr/bin/env python3

from swat_utilities.operations import import_mgt_excel
    
file_path = '../resources/rotations/NSF_BEACON.xlsx'
discard = ['Conversions', 'Reference Sheet']
names = None

rotations = import_mgt_excel(file_path, sheets=names, remove=discard)
