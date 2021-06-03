#!/usr/bin/env python3
"""
Classes for building schedules of management operations for SWAT
"""
import pandas as pd
import collections


class Operation:

    def __init__(self, mgt_op, mgt_params, name):

        if mgt_op == 1:  # Planting/Beginning of growing season
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'PLANT_ID': None,
                              'CURYR_MAT': '',
                              'HEAT_UNITS': 1800,
                              'LAI_INIT': 0,
                              'HI_TARG': 0,
                              'BIO_INIT': 0,
                              'BIO_TARG': 0,
                              'CNOP': 0}
            line = " {:2d} {:2d} {:8.3f} {:2d} {:4d}     {:2d} {:12.5f} {:6.2f} {:11.5f} {:4.2f} {:6.2f} {:5.2f}"
        elif mgt_op == 3:  # Fertilizer application
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'FERT_ID': None,
                              'FRT_KG': None,
                              'FRT_SURFACE': 0}
            line = " {:2d} {:2d} {:8.3f} {:2d} {:4d}        {:12.5f} {:6.2f}"
        elif mgt_op == 4:  # Pesticide application
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'PEST_ID': None,
                              'PST_KG': None,
                              'PST_DEP': 0}
            line = " {:2d} {:2d} {:8.3f} {:2d} {:4d}        {:12.5f} {:6.2f}"
        elif mgt_op == 5:  # Harvest and kill operation
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'CNOP': 0,
                              'HI_OVR': '',
                              'FRAC_HARVK': ''}
            line = " {:2d} {:2d} {:8.3f} {:2d}             {:12.5f} {:6.2f} {:11.5f}"
        elif mgt_op == 6:  # Tillage operation
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'TILL_ID': None,
                              'CNOP': 0}
            line = " {:2d} {:2d} {:8.3f} {:2d} {:4d}        {:12.5f}"
        elif mgt_op == 7:  # Harvest operation
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op,
                              'IHV_GBM': 0,
                              'HARVEFF': '',
                              'HI_OVR': ''}
            line = " {:2d} {:2d} {:8.3f} {:2d}      {:3d}    {:12.5f} {:6.2f}"
        elif mgt_op == 8:  # Kill operation
            mgt_params_ref = {'MONTH': None,
                              'DAY': None,
                              'HUSC': '',
                              'MGT_OP': mgt_op}
            line = " {:2d} {:2d} {:8.3f} {:2d}"
        elif mgt_op == 17:  # Skip a year operation
            mgt_params_ref = {'MGT_OP': mgt_op}
            line = "                {:2d}"
        else:
            print('The management operation you specified is not available yet')
            return

        mgt_params = set_params(mgt_params, mgt_params_ref)

        mgt_values = list(mgt_params.values())
        if len(mgt_values) > 1:
            if (mgt_values[0] is None) or (mgt_values[1] is None) and (len(mgt_values[2]) == 0):
                print('You must set MONTH and DAY or HUSC (fraction of total base zero heat units)')
                return
            if None in mgt_values:
                fix = [k for k, v in mgt_params.items() if v is None]
                print('You must provide a value for:', *fix, sep='\n')
                return

        self.mgt_params = mgt_params
        self.line_ref = line

        line = new_line(line, mgt_values)  # Modify line formatting based on empty values
        self.line = line.format(*mgt_values)

        self.name = name

    def set_mgt_param(self, params):
        mgt_params_ref = self.mgt_params
        self.mgt_params = set_params(params, mgt_params_ref)

        # update line
        line = new_line(self.line_ref, list(self.mgt_params.values()))
        self.line = line.format(*list(self.mgt_params.values()))

        return self.mgt_params

    def get_mgt_param(self, params):
        values = {k: self.mgt_params[k] for k in params}
        return values

    def mgt_param_list(self):
        list1 = list(self.mgt_params.keys())
        print(*list1, sep='\n')


class Schedule:

    def __init__(self, operations, name):
        self.operations = sort_operations(operations)
        self.schedule = [operation.line for operation in self.operations]
        self.name = name

    def add(self, operations):
        operations = list(flatten([self.operations, operations]))
        self.operations = sort_operations(operations)
        self.schedule = [operation.line for operation in self.operations]
        return self

    def remove(self, names):
        for name in names:
            mask = [i for i, x in enumerate(self.operations) if x.name == name]
            if not mask:
                continue
            del self.operations[mask[0]]
        self.schedule = [operation.line for operation in self.operations]
        return self

    def edit(self, operations):
        names = operations.keys()
        names_ref = [operation.name for operation in self.operations]
        for name in names:
            mask = [i for i, x in enumerate(names_ref) if x == name][0]
            self.operations[mask].set_mgt_param(operations[name])

        self.operations = sort_operations(self.operations)
        self.schedule = [operation.line for operation in self.operations]
        return self

    def get(self, operations):
        names = operations.keys()
        names_ref = [operation.name for operation in self.operations]
        values = {}
        for name in names:
            mask = [i for i, x in enumerate(names_ref) if x == name][0]
            operation = self.operations[mask]
            values.update({name: operation.get_mgt_param(operations[name])})
        return values

    def print(self):
        print(*self.schedule, sep='\n')


class Rotation:

    def __init__(self, sequence, name=''):
        self.name = name
        self.nrot = len(sequence)
        self.sequence = sequence
        self.schedules = unique(sequence)

        op_sched = []
        for item in sequence:
            op_sched.append(item.schedule)
            op_sched.append(Operation(17, {}, '').line)

        self.op_sched = list(flatten(op_sched))

    def print(self):
        print(*self.op_sched, sep='\n')


def set_params(params, mgt_params_ref):
    ref_keys = list(mgt_params_ref.keys())
    set_keys = list(params.keys())

    for key in set_keys:
        if key in ref_keys:
            mgt_params_ref[key] = params[key]
        else:
            print('Management parameter {:s} is not defined for the selected operation'.format(key))

    return mgt_params_ref


def new_line(line, values):
    line_list = line.split(' ')
    index = list()
    c = 0
    for item in line_list:
        if item != '':
            index.append(c)
            c += 1
        else:
            index.append(None)

    for ind, item in enumerate(values):
        if item == '':
            mask = [i for i, x in enumerate(index) if x == ind][0]
            old_format = line_list[mask]
            new_format = '{:' + old_format.split('.')[0].split(':')[1].split('d')[0] + 's}'
            line_list[mask] = new_format

    line = ' '.join(line_list)
    return line


def sort_operations(operations):
    dates = [(i, operation.mgt_params['MONTH'], operation.mgt_params['DAY']) for i, operation in
             enumerate(operations)]
    dates = pd.DataFrame(dates, columns=['ID', 'month', 'days']).sort_values(['month', 'days'],
                                                                             ascending=[True, True])
    order = dates['ID'].tolist()

    return [operations[i] for i in order]


def flatten(lt):
    for el in lt:
        if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)

    return unique_list


def import_mgt_excel(path_to_file, sheets=None, remove=None):
    raw_data = pd.read_excel(path_to_file, sheet_name=sheets, header=0, skiprows=2)
    
    if (remove is not None) & (sheets is None):
        for i in remove:
            del raw_data[i]
    
    rot_dict = dict()            
    for name, rot in raw_data.items():
        mgt_ids = rot['mgt op'].tolist()
        ind = [-1] + [i for i, x in enumerate(mgt_ids) if x == 17]
        
        sched_list = []
        for i in range(0, len(ind)-1):
            crop = rot.iloc[(ind[i]+1):(ind[i+1]), 1:]
            sched_list.append(build_sched(crop, i))
        
        rot_dict[name] = Rotation(sched_list, name=name)        
     
    return rot_dict


def build_sched(df, ind):
    mgt_ids = df['mgt op']
    counter = collections.Counter(mgt_ids)
    
    operations_list = []
    for i in counter.keys():
        operations = df.loc[df['mgt op'] == i, :]
        c = 0
        for _, op in operations.iterrows():
            operations_list.append(build_operation(op, c))
            c += 1
    sched = Schedule(operations_list, 'crop_' + str(int(ind)))
    return sched


def build_operation(df, ind):
    
    mgt_op = int(df['mgt op'])
    
    if mgt_op == 1:  # Planting/Beginning of growing season
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'PLANT_ID': int(df['mgt1i']),
                      'CURYR_MAT': '' if pd.isna(df['mgt3i']) else int(df['mgt3i']),
                      'HEAT_UNITS': 1800 if pd.isna(df['mgt4']) else df['mgt4'],
                      'LAI_INIT': 0 if pd.isna(df['mgt5']) else df['mgt5'],
                      'HI_TARG': 0 if pd.isna(df['mgt6']) else df['mgt6'],
                      'BIO_INIT': 0 if pd.isna(df['mgt7']) else df['mgt7'],
                      'BIO_TARG': 0 if pd.isna(df['mgt8']) else df['mgt8'],
                      'CNOP': 0 if pd.isna(df['mgt9']) else df['mgt9']}
        name = 'Plant_' + str(int(ind))
    elif mgt_op == 3:  # Fertilizer application
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'FERT_ID': int(df['mgt1i']),
                      'FRT_KG': df['mgt4'],
                      'FRT_SURFACE': 0 if pd.isna(df['mgt5']) else df['mgt5']}
        name = 'Fertilizer_' + str(int(ind))
    elif mgt_op == 4:  # Pesticide application
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'PEST_ID': int(df['mgt1i']),
                      'PST_KG': df['mgt4'],
                      'PST_DEP': 0 if pd.isna(df['mgt5']) else df['mgt5']}
        name = 'Pesticide_' + str(int(ind))
    elif mgt_op == 5:  # Harvest and kill operation
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'CNOP': 0 if pd.isna(df['mgt4']) else df['mgt4'],
                      'HI_OVR': '' if pd.isna(df['mgt5']) else df['mgt5'],
                      'FRAC_HARVK': '' if pd.isna(df['mgt6']) else df['mgt6']}
        name = 'Harvest/Kill_' + str(int(ind))
    elif mgt_op == 6:  # Tillage operation
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'TILL_ID': int(df['mgt1i']),
                      'CNOP': 0 if pd.isna(df['mgt4']) else df['mgt4']}
        name = 'Tillage_' + str(int(ind))
    elif mgt_op == 7:  # Harvest operation
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op,
                      'IHV_GBM': 0 if pd.isna(df['mgt2i']) else int(df['mgt2i']),
                      'HARVEFF': '' if pd.isna(df['mgt4']) else df['mgt4'],
                      'HI_OVR': '' if pd.isna(df['mgt5']) else df['mgt5']}
        name = 'Harvest_' + str(int(ind))
    elif mgt_op == 8:  # Kill operation
        mgt_params = {'MONTH': '' if pd.isna(df['mon']) else int(df['mon']),
                      'DAY': '' if pd.isna(df['day'])  else int(df['day']),
                      'HUSC': '' if pd.isna(df['HU']) else int(df['HU']),
                      'MGT_OP': mgt_op}
        name = 'Kill_' + str(int(ind))
    else:
        print('The management operation you specified is not available yet')
        return
    
    operation = Operation(mgt_op, mgt_params, name)
    return operation