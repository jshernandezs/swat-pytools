#!/usr/bin/env python3
"""
Create crop rotations using a set of options
"""

import datetime as dt
import copy
from swat_utilities.operations import Operation, Schedule, Rotation
from swat_utilities.ops_operations import OpsOperation, OpsSchedule


def set_rot(opt_rot, opt_sys, opt_rate, opt_cc, opt_td, opt_fs, params_fs, base_rates=None, rates_multipliers=None,
            cc_id=None, cc_hu=None, kill_cc_date=None):
    # Base annual management operations (schedules)
    till1 = Operation(6, {'MONTH': 5, 'DAY': 1, 'TILL_ID': 6}, 'Tillage1')
    till2 = Operation(6, {'MONTH': 11, 'DAY': 14, 'TILL_ID': 60}, 'Tillage2')
    fertilizer_n = Operation(3, {'MONTH': 5, 'DAY': 5, 'FERT_ID': 1, 'FRT_KG': 0}, 'Nitrogen')
    fertilizer_p = Operation(3, {'MONTH': 5, 'DAY': 5, 'FERT_ID': 2, 'FRT_KG': 0}, 'Phosphorus')
    plant = Operation(1, {'MONTH': 5, 'DAY': 5, 'PLANT_ID': 19, 'HEAT_UNITS': 1217}, 'Plant')
    kill = Operation(5, {'MONTH': 11, 'DAY': 1}, 'Kill')

    corn = set_system(Schedule([till1, fertilizer_p, fertilizer_n, plant, kill, till2], 'corn_0'), opt_sys,
                      fert_name='Phosphorus', plant_name='Plant', kill_name='Kill', till_name1='Tillage1',
                      till_name2='Tillage2', delta=[0, 1])
    corn = set_cover_crop(corn, opt_cc, cc_id, cc_hu, kill_cc_date, 14, 'Kill')

    till1 = Operation(6, {'MONTH': 5, 'DAY': 14, 'TILL_ID': 6}, 'Tillage1')
    fertilizer_p = Operation(3, {'MONTH': 5, 'DAY': 14, 'FERT_ID': 2, 'FRT_KG': 0}, 'Phosphorus')
    plant = Operation(1, {'MONTH': 5, 'DAY': 15, 'PLANT_ID': 56}, 'Plant')
    kill = Operation(5, {'MONTH': 10, 'DAY': 1}, 'Kill')    

    soybn = set_system(Schedule([till1, fertilizer_p, plant, kill], 'soybn_0'), opt_sys,
                       fert_name='Phosphorus', plant_name='Plant', kill_name='Kill', till_name1='Tillage1',
                       delta=[0, 1])
    soybn = set_cover_crop(soybn, opt_cc, cc_id, cc_hu, kill_cc_date, 14, 'Kill')

    if opt_rot == 1:  # Corn-soybean rotation
        name_rot = 'corn-soybean'

        corn = set_fert_rates(corn, 'Phosphorus', 'phosphorus', opt_rate, 'corn', name_rot,
                              base_rates, rates_multipliers)
        corn = set_fert_rates(corn, 'Nitrogen', 'nitrogen', opt_rate, 'corn', name_rot,
                              base_rates, rates_multipliers)

        soybn = set_fert_rates(soybn, 'Phosphorus', 'phosphorus', opt_rate, 'soybean', name_rot,
                               base_rates, rates_multipliers)

        rot = Rotation([corn, soybn], name_rot)
    elif opt_rot == 2:  # Corn-soybean-wheat rotation
        name_rot = 'corn-soybean-wheat'
        
        corn.remove(['Tillage2'])
        till2 = Operation(6, {'MONTH': 11, 'DAY': 14, 'TILL_ID': 56}, 'Tillage2')
        corn.add([till2])
        corn = set_system(corn, opt_sys, till_name2='Tillage2')

        till2 = Operation(6, {'MONTH': 10, 'DAY': 2, 'TILL_ID': 6}, 'Tillage2')
        fertilizer_p = Operation(3, {'MONTH': 10, 'DAY': 7, 'FERT_ID': 2, 'FRT_KG': 0}, 'Phosphorus2')
        plant = Operation(1, {'MONTH': 10, 'DAY': 7, 'PLANT_ID': 28, 'HEAT_UNITS': 1549}, 'Plant2')
        soybn.add([till2, fertilizer_p, plant])
        soybn.remove(['plant cover crop'])
        soybn = set_system(soybn, opt_sys, fert_name='Phosphorus2', 
                           plant_name='Plant2', till_name2='Tillage2', delta=[0, 1])

        fertilizer_n = Operation(3, {'MONTH': 4, 'DAY': 15, 'FERT_ID': 1, 'FRT_KG': 0}, 'Nitrogen')
        kill = Operation(5, {'MONTH': 7, 'DAY': 15}, 'Kill')
        till1 = Operation(6, {'MONTH': 10, 'DAY': 14, 'TILL_ID': 60}, 'Tillage1')
        wheat = Schedule([fertilizer_n, kill, till1, fertilizer_p], 'wheat_0')

        wheat = set_system(wheat, opt_sys, fert_name='Phosphorus2', 
                           kill_name='Kill', till_name1='Tillage1', delta=[0, 1])
        wheat = set_cover_crop(wheat, opt_cc, cc_id, cc_hu, kill_cc_date, plant_timing=92, source_date='Kill')
        wheat.remove(['kill cover crop'])

        check = wheat.get({'Phosphorus2': ['MONTH', 'DAY']})
        if (check['Phosphorus2']['MONTH'] == fertilizer_p.get_mgt_param(['MONTH'])['MONTH']) & \
                (check['Phosphorus2']['DAY'] == fertilizer_p.get_mgt_param(['DAY'])['DAY']):
            wheat.remove(['Phosphorus2'])
            soybn = set_fert_rates(soybn, 'Phosphorus2', 'phosphorus', opt_rate, 'wheat', name_rot,
                                   base_rates, rates_multipliers)
        else:
            soybn.remove(['Phosphorus2'])
            wheat = set_fert_rates(wheat, 'Phosphorus2', 'phosphorus', opt_rate, 'wheat', name_rot,
                                   base_rates, rates_multipliers)

        corn = set_fert_rates(corn, 'Phosphorus', 'phosphorus', opt_rate, 'corn', name_rot,
                              base_rates, rates_multipliers)
        corn = set_fert_rates(corn, 'Nitrogen', 'nitrogen', opt_rate, 'corn', name_rot,
                              base_rates, rates_multipliers)

        soybn = set_fert_rates(soybn, 'Phosphorus', 'phosphorus', opt_rate, 'soybean', name_rot,
                               base_rates, rates_multipliers)

        wheat = set_fert_rates(wheat, 'Nitrogen', 'nitrogen', opt_rate, 'wheat', name_rot,
                               base_rates, rates_multipliers)

        rot = Rotation([corn, soybn, wheat], name_rot)
    else:
        print('Rotation option not available')
        return

    td_params = set_td(opt_td)
    td_params.update({'NROT': [rot.nrot, 'replace', 'mgt']})

    if opt_fs == 1:
        mgt_op = 10
        params_fs = {k: params_fs[k] for k in ('MONTH', 'DAY', 'IYEAR')}
        name = 'no filter strip'
    elif opt_fs == 2:
        mgt_op = 4
        name = 'filter strip'
    else:
        print('Filter strip option not available')
        return

    fs = OpsOperation(mgt_op, params_fs, 'filter strip')
    filter_strip = OpsSchedule([fs], name)

    return rot, td_params, filter_strip


def set_system(sched, opt_sys, fert_name=None, plant_name=None, kill_name=None, till_name1=None, till_name2=None,
               delta=None):
    if delta is None:
        delta = [0, 0]

    sched_new = copy.deepcopy(sched)

    if opt_sys == 1:
        if fert_name is not None:
            if plant_name is not None:
                sched_new = set_date(sched_new, plant_name, fert_name, 'before', delta[0])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0}})
    elif opt_sys == 2:
        if fert_name is not None:
            if kill_name is not None:
                sched_new = set_date(sched_new, kill_name, fert_name, 'after', delta[1])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0}})
    elif opt_sys == 3:
        if fert_name is not None:
            if plant_name is not None:
                sched_new = set_date(sched_new, plant_name, fert_name, 'before', delta[0])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0}})
        if till_name1 is not None:
            sched_new.edit({till_name1: {'TILL_ID': 17}})
        if till_name2 is not None:
            if sched.name == 'corn_0':
                sched_new.edit({till_name2: {'TILL_ID': 60}})
            elif sched.name == 'soybn_0':
                sched_new.edit({till_name2: {'TILL_ID': 6}})
            else:           
                sched_new.edit({till_name2: {'TILL_ID': 3}})            
    elif opt_sys == 4:
        if fert_name is not None:
            if kill_name is not None:
                sched_new = set_date(sched_new, kill_name, fert_name, 'after', delta[1])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0}})
        if till_name1 is not None:
            sched_new.edit({till_name1: {'TILL_ID': 17}})
        if till_name2 is not None:
            if sched.name == 'corn_0':
                sched_new.edit({till_name2: {'TILL_ID': 60}})
            elif sched.name == 'soybn_0':
                sched_new.edit({till_name2: {'TILL_ID': 6}})
            else:           
                sched_new.edit({till_name2: {'TILL_ID': 3}})
    elif opt_sys == 5:
        if fert_name is not None:
            if plant_name is not None:
                sched_new = set_date(sched_new, plant_name, fert_name, 'before', delta[0])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0.01}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 3}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    elif opt_sys == 6:
        if fert_name is not None:
            if kill_name is not None:
                sched_new = set_date(sched_new, kill_name, fert_name, 'after', delta[1])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0.01}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 3}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    elif opt_sys == 7:
        if fert_name is not None:
            if plant_name is not None:
                sched_new = set_date(sched_new, plant_name, fert_name, 'before', delta[0])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0.01}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 4}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    elif opt_sys == 8:
        if fert_name is not None:
            if kill_name is not None:
                sched_new = set_date(sched_new, kill_name, fert_name, 'after', delta[1])
            sched_new.edit({fert_name: {'FRT_SURFACE': 0.01}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 4}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    elif opt_sys == 9:
        if fert_name is not None:
            if plant_name is not None:
                sched_new = set_date(sched_new, plant_name, fert_name, 'before', delta[0])
            sched_new.edit({fert_name: {'FRT_SURFACE': 1}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 4}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    elif opt_sys == 10:
        if fert_name is not None:
            if kill_name is not None:
                sched_new = set_date(sched_new, kill_name, fert_name, 'after', delta[1])
            sched_new.edit({fert_name: {'FRT_SURFACE': 1}})
        if till_name1 is not None:
            if sched.name == 'wheat_0':
                sched_new.remove(['Tillage1'])
            else:
                sched_new.edit({till_name1: {'TILL_ID': 4}})
        if till_name2 is not None:
            sched_new.remove(['Tillage2'])
    else:
        print('Option is not available yet')
        return

    return sched_new


def set_date(schedule, source, target, timing, delta=0):
    schedule = copy.deepcopy(schedule)

    date = schedule.get({source: ['MONTH', 'DAY']})
    date_old = dt.date(2001, date[source]['MONTH'], date[source]['DAY'])

    if timing == 'before':
        date_new = date_old - dt.timedelta(delta)
    elif timing == 'after':
        date_new = date_old + dt.timedelta(delta)
    else:
        print('timing should be "after" or "before"')
        return

    schedule.edit({target: {'MONTH': date_new.month, 'DAY': date_new.day}})
    return schedule


def set_fert_rates(schedule, target, fert_name, opt_rate, crop, rotation_name, rates=None, multipliers=None):
    if rates is None:
        rates = {'corn-soybean': {'corn': [129, 24.64],
                                  'soybean': [0, 14.21]},
                 'corn-soybean-wheat': {'corn': [162.5, 24.64],
                                        'soybean': [0, 14.21],
                                        'wheat': [60.5, 22.70]}}

    if multipliers is None:
        multipliers = {'nitrogen': [1, 1, 1, 1, 1],
                       'phosphorus': [0, 0.5, 1, 1.5, 2]}

    sched_new = copy.deepcopy(schedule)

    if (opt_rate > len(multipliers[fert_name.lower()])) or (opt_rate < 1):
        print('opt_rate out of range')
        return

    if fert_name.lower() == 'nitrogen':
        ind = 0
    elif fert_name.lower() == 'phosphorus':
        ind = 1
    else:
        print('fert_name must be either "nitrogen" or "phosphorus')
        return

    list_rotations = list(rates.keys())
    if rotation_name.lower() not in list_rotations:
        print('Check that the rotation name you are entering is in the following list:', *list_rotations, sep='\n')
        return

    list_crops = list(rates[rotation_name].keys())
    if crop.lower() not in list_crops:
        print('Check that the crop you are entering is in the following list:', *list_crops, sep='\n')
        return

    rate = multipliers[fert_name.lower()][opt_rate - 1] * rates[rotation_name.lower()][crop.lower()][ind]

    sched_new.edit({target: {'FRT_KG': rate}})
    return sched_new


def set_cover_crop(schedule, opt_cc, cc_id=None, cc_hu=None, kill_cc_date=None, plant_timing=None, source_date=None):
    if cc_id is None:
        cc_id = [30, 28, 178, 32]

    if cc_hu is None:
        cc_hu = [1250, 1549, 750, 1628]

    if plant_timing is None:
        plant_timing = 3

    if kill_cc_date is None:
        kill_cc_date = {'MONTH': 4, 'DAY': 14}

    if source_date is None:
        source_date = 'Kill'
        
    if (opt_cc > (len(cc_id) + 1)) or (opt_cc < 1):
        print('opt_cc out of range')
        return
        
    if opt_cc == 1:
        return schedule
    else:
        sched_new = copy.deepcopy(schedule)
        kill = Operation(8, kill_cc_date, 'kill cover crop')
    
        plant = Operation(1, {'MONTH': 10, 'DAY': 14, 'PLANT_ID': cc_id[opt_cc - 2], 'HEAT_UNITS': cc_hu[opt_cc - 2]},
                          'plant cover crop')
    
        sched_new.add([plant, kill])
        sched_new = set_date(sched_new, source_date, 'plant cover crop', 'after', plant_timing)
    
        return sched_new


# noinspection PyUnboundLocalVariable
def set_td(opt_tile):

    iflag = 0
    if opt_tile == 1:
        ddrain = 0  # (mm)
        iflag = 1
    elif opt_tile == 2:
        itdrn = 1
        iwtdn = 1
        sol_p_model = 0
        ismax = 0        
        re = 50  # (mm)
        sdrain = 15000  # (mm)
        drain_co = 12  # (mm day-1)
        pc = 0  # (mm hr-1)
        latksatf = 1
        sstmaxd = 20  # (mm)
    else:
        print('Option for tile drainage not available yet')
        return

    if iflag == 1:
        params = {'DDRAIN': [ddrain, 'replace', 'mgt']}
    else:
        params = {'ITDRN': [itdrn, 'replace', 'bsn'],
                  'IWTDN': [iwtdn, 'replace', 'bsn'],
                  'SOL_P_MODEL': [sol_p_model, 'replace', 'bsn'],
                  'ISMAX': [ismax, 'replace', 'bsn'],                  
                  'RE': [re, 'replace', 'sdr'],
                  'SDRAIN': [sdrain, 'replace', 'sdr'],
                  'DRAIN_CO': [drain_co, 'replace', 'sdr'],
                  'PC': [pc, 'replace', 'sdr'],
                  'LATKSATF': [latksatf, 'replace', 'sdr'],
                  'SSTMAXD': [sstmaxd, 'replace', 'sdr']}

    return params
