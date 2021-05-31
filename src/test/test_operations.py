#!/usr/bin/env python3
import sys
import os

main_dir = os.path.abspath('../main/')
sys.path.insert(1, main_dir)

from swat_utilities.operations import Operation, Schedule, Rotation

# creating operations

op = 1
param = {'MONTH': 5,
         'DAY': 5,
         'PLANT_ID': 19}
name = 'plant'

plant = Operation(op, param, name)

plant.set_mgt_param({'HEAT_UNITS': 1217})  # editing an existing operation

tillage = Operation(mgt_op=6,
                    mgt_params={'MONTH': 5,
                                'DAY': 1,
                                'TILL_ID': 6},
                    name='tillage')

# corn schedule

corn = Schedule([plant, tillage], 'corn_0')

fertilizer_N = Operation(3, {'MONTH': 5, 'DAY': 4, 'FERT_ID': 4, 'FRT_KG': 194}, 'urea')
fertilizer_P = Operation(3, {'MONTH': 5, 'DAY': 5, 'FERT_ID': 2, 'FRT_KG': 25}, 'phosphorus')
kill = Operation(5, {'MONTH': 11, 'DAY': 1}, 'kill')
chisel = Operation(6, {'MONTH': 11, 'DAY': 15, 'TILL_ID': 60}, 'chisel')

corn.add([fertilizer_N, fertilizer_P, kill, chisel])

corn.edit({'phosphorus': {'FRT_KG': 59}})

# soybean schedule

plant = Operation(1, {'MONTH': 5, 'DAY': 15, 'PLANT_ID': 56}, 'plant')
tillage = Operation(6, {'MONTH': 5, 'DAY': 14, 'TILL_ID': 6}, 'tillage')
fertilizer_P = Operation(3, {'MONTH': 5, 'DAY': 14, 'FERT_ID': 2, 'FRT_KG': 45}, 'phosphorus')
kill = Operation(5, {'MONTH': 10, 'DAY': 1}, 'kill')
chisel = Operation(6, {'MONTH': 10, 'DAY': 30, 'TILL_ID': 60}, 'chisel')

soybn = Schedule([plant, tillage, chisel, fertilizer_P, kill], 'soybn_0')

# Corn-soybean rotation

rotation = Rotation([corn, corn, corn, corn, corn, soybn], 'Corn')
rotation.print()
