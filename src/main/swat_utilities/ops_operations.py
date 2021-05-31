#!/usr/bin/env python3
"""
Classes for building schedules of 'ops' operations for SWAT
"""
import pandas as pd
import collections.abc as collections


class OpsOperation:

    def __init__(self, mgt_op, ops_params, name):

        if mgt_op == 4:  # Filter strip
            ops_params_ref = {'MONTH': None,
                              'DAY': None,
                              'IYEAR': None,
                              'MGT_OP': mgt_op,
                              'FILTER_I': 1,
                              'FILTER_RATIO': 50,
                              'FILTER_CON': 0.5,
                              'FILTER_CH': 0.25}
            line = ' {:2d} {:2d}     {:4d} {:2d} {:4d}     {:6.2f} {:12.5f} {:6.2f}'
        elif mgt_op == 10:  # Generic operation
            ops_params_ref = {'MONTH': None,
                              'DAY': None,
                              'IYEAR': None,
                              'MGT_OP': mgt_op,
                              'RO_NMP_FLAG': 0,
                              'RO_NMP_SED': 0,
                              'RO_NMP_PP': 0,
                              'RO_NMP_SP': 0,
                              'RO_NMP_PN': 0,
                              'RO_NMP_SN': 0,
                              'RO_NMP_BAC': 0}
            line = ' {:2d} {:2d}     {:4d} {:2d} {:4d}     {:6.2f} {:12.5f} {:6.2f}             {:8.2f} {:6.2f} {:5.2f}'
        else:
            print('The ops management operation you specified is not available yet')
            return

        ops_params = set_params(ops_params, ops_params_ref)

        ops_values = list(ops_params.values())
        if len(ops_values) > 1:
            if (ops_values[0] is None) or (ops_values[1] is None) or (ops_values[2] is None):
                print('You must set MONTH, DAY and IYEAR')
                return
            if None in ops_values:
                fix = [k for k, v in ops_params.items() if v is None]
                print('You must provide a value for:', *fix, sep='\n')
                return

        self.ops_params = ops_params
        self.line_ref = line

        line = new_line(line, ops_values)  # Modify line formatting based on empty values
        self.line = line.format(*ops_values)

        self.name = name

    def set_ops_param(self, params):
        ops_params_ref = self.ops_params
        self.ops_params = set_params(params, ops_params_ref)

        # update line
        line = new_line(self.line_ref, list(self.ops_params.values()))
        self.line = line.format(*list(self.ops_params.values()))

        return self.ops_params

    def get_ops_param(self, params):
        values = {k: self.ops_params[k] for k in params}
        return values

    def ops_param_list(self):
        list1 = list(self.ops_params.keys())
        print(*list1, sep='\n')


class OpsSchedule:

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
            mask = [i for i, x in enumerate(self.operations) if x.name == name][0]
            del self.operations[mask]
        self.schedule = [operation.line for operation in self.operations]
        return self

    def edit(self, operations):
        names = operations.keys()
        names_ref = [operation.name for operation in self.operations]
        for name in names:
            mask = [i for i, x in enumerate(names_ref) if x == name][0]
            self.operations[mask].set_ops_param(operations[name])

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
            values.update({name: operation.get_ops_param(operations[name])})
        return values

    def print(self):
        print(*self.schedule, sep='\n')


def set_params(params, ops_params_ref):
    ref_keys = list(ops_params_ref.keys())
    set_keys = list(params.keys())

    for key in set_keys:
        if key in ref_keys:
            ops_params_ref[key] = params[key]
        else:
            print('Ops management parameter {:s} is not defined for the selected operation'.format(key))

    return ops_params_ref


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
    dates = [(i, operation.ops_params['MONTH'], operation.ops_params['DAY'], operation.ops_params['IYEAR']) for
             i, operation in enumerate(operations)]
    dates = pd.DataFrame(dates, columns=['ID', 'month', 'days', 'year']).sort_values(['year', 'month', 'days'],
                                                                                     ascending=[True, True, True])
    order = dates['ID'].tolist()

    return [operations[i] for i in order]


def flatten(lt):
    for el in lt:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)

    return unique_list
