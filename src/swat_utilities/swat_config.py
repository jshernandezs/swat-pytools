#!/usr/bin/env python3
"""
SWAT utilities for changing text input files and executing the model
"""
import os
import platform
import re
import subprocess as sp
import pandas as pd
import sys
import datetime as dt


class RunID:
    """ create object containing the identifiers for each new SWAT run
    """
    __ids = []

    @staticmethod
    def get_id():
        possible_id = 0
        while possible_id < len(RunID.__ids) and RunID.__ids[possible_id]:
            possible_id += 1
        if possible_id >= len(RunID.__ids):
            RunID.__ids.append(True)
        else:
            RunID.__ids[possible_id] = True
        return possible_id

    @staticmethod
    def remove_id(ind):
        if ind < len(RunID.__ids):
            RunID.__ids[ind] = False

    @staticmethod
    def clear_id():
        RunID.__ids = []


class ModelSetup:
    """ create object with SWAT settings for running a model
    """

    def __init__(self, model_file):
        self.output_dir = '/tmp/output_swat'
        self.temp_dir = '/tmp/swat_runs'
        self.model_file = os.path.abspath(model_file)
        self.agr_treshold = ''
        self.hid_treshold = ''
        self.swat_dir = os.path.abspath('../resources')
        self.swat_exec_name = 'SWAT_Rev670'
        self.new_model_name = 'New_SWAT'
        self.subbasins = []
        self.hrus = []
        self.param = []
        self.out_reach = []
        self.fig_file = ''
        self.conc_output = []
        self.sim_ts = {}
        self.run_output_dir = ''
        self.verbose = True
        self.verbose_swat = False
        self.executed = False
        self.saveconc_files = []
        self.drainmod = False
        self.ops_table = []
        self.mgt_table = []
        self.subbasins_ops = []
        self.hrus_ops = []
        self.subbasins_mgt = []
        self.hrus_mgt = []

    def get_subbasins(self):

        out_reaches = self.out_reach
        fig_file = self.fig_file

        if len(out_reaches) == 0:
            subbasins_list = []
        else:
            subbasin_list = []
            route_list = []
            add_list = []

            if len(fig_file) == 0:
                print(
                    'Warning: You must provide the full path to the .fig file to get the subbasins for the watershed. '
                    'Changes (if any) will be applied to all subbasins...')
                subbasins_list = []
            else:
                with open(fig_file) as f:
                    lines = f.readlines()
                    for line in lines:
                        temp = line.split()
                        if temp[0] == "subbasin":
                            subbasin_list.append(temp)
                        elif temp[0] == "route":
                            route_list.append(temp)
                        elif temp[0] == "add":
                            add_list.append(temp)

                    subbasin_list = pd.DataFrame(subbasin_list)
                    route_list = pd.DataFrame(route_list)
                    add_list = pd.DataFrame(add_list)

                subbasins_list = []

                for out_reach in out_reaches:
                    subbasins = []
                    # find inflow element to out_reach in route_list
                    inflow = int(route_list.loc[pd.to_numeric(route_list[3]) == out_reach].to_numpy()[:, -1:])
                    # get list of subbasins draining to out_reach
                    subbasins = subbasin_root(inflow, subbasin_list, route_list, add_list, subbasins)
                    # append to general list
                    subbasins_list.append(subbasins)

        self.subbasins = subbasins_list

    def prepare_swat(self):
        # get attributes from SWAT configuration object
        param = self.param
        subbasins = self.subbasins
        hrus = self.hrus
        output_dir = self.output_dir
        temp_dir = self.temp_dir
        model_file = self.model_file
        new_model_name = self.new_model_name
        swat_dir = self.swat_dir
        swat_exec_name = self.swat_exec_name
        drainmod = self.drainmod
        mgt_table = self.mgt_table
        subbasins_mgt = self.subbasins_mgt
        hrus_mgt = self.hrus_mgt
        ops_table = self.ops_table
        subbasins_ops = self.subbasins_ops
        hrus_ops = self.hrus_ops

        # create output directory if it does not exist

        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except OSError:
            print('Tried to create directory {:s} that already exists...'.format(output_dir))

        try:
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
        except OSError:
            print('Tried to create directory {:s} that already exists...'.format(temp_dir))

        # unzip SWAT model in temp directory
        input_dir_file = os.path.abspath(model_file)
        model_dir_file = os.path.abspath(temp_dir + '/' + new_model_name)

        if not os.path.exists(model_dir_file):
            os.makedirs(model_dir_file)

        sp.run(['unzip', input_dir_file, '-d', model_dir_file], stdout=sp.DEVNULL, stderr=sp.DEVNULL, check=True)

        if drainmod:
            update_txt(model_dir_file)

        # copy/update SWAT executable to directory containing text input files
        swat_dir_file = os.path.abspath(swat_dir + '/' + swat_exec_name)

        if platform.system() == 'Linux':
            sp.run(['cp', swat_dir_file, model_dir_file, '--update'], check=True)
        elif platform.system() == 'Darwin':
            sp.run(['cp', '-f', swat_dir_file, model_dir_file], check=True)

        # prepare SWAT input text files
        sp.run(['cp', '-r', model_dir_file, output_dir], check=True)
        output_dir = os.path.abspath(output_dir + '/' + new_model_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # write SWAT input text files
        write_new_files(param, subbasins, hrus, model_dir_file, output_dir)
        write_mgt_tables(mgt_table, 'mgt', subbasins_mgt, hrus_mgt, output_dir)
        write_mgt_tables(ops_table, 'ops', subbasins_ops, hrus_ops, output_dir)

        # remove SWAT model files from temporal directory
        sp.run(['rm', '-rf', model_dir_file], check=True)

        self.run_output_dir = output_dir

        return output_dir

    def run_swat(self):
        """ Function to execute SWAT in linux through python3
        input(s):
            model_dir = directory for the SWAT TxtInOut folder
            swat_exec_name = name of SWAT executable
        """
        # get attributes from SWAT configuration object

        swat_exec_name = self.swat_exec_name
        model_dir = self.run_output_dir

        if len(model_dir) == 0:
            print('You must execute method ''prepareSWAT'' before executing ''runSWAT'' method...')
        else:
            # run SWAT
            if self.verbose:
                print("Executing SWAT model {:s} in directory {:s}...".format(self.new_model_name, model_dir))
            if self.verbose_swat:
                sp.run('./' + swat_exec_name, cwd=model_dir, check=True)
            else:
                sp.run('./' + swat_exec_name, cwd=model_dir, stdout=sp.DEVNULL, stderr=sp.DEVNULL, check=True)
            self.executed = True

    def remove_swat(self):

        model_dir = self.run_output_dir

        if len(model_dir) == 0:
            print('Model folder was not found...')
        else:
            # remove SWAT model files from directory
            sp.run(['rm', '-rf', model_dir], check=True)

    def read_output(self):

        model_dir = self.run_output_dir

        if not self.executed:
            print('You must execute method ''runSWAT'' method first...')
        else:
            self.fig_file = model_dir + '/fig.fig'
            # get 'watout' files from configuration file
            with open(self.fig_file) as f:
                lines = f.readlines()
                saveconc_files = []
                for i, line in enumerate(lines):
                    temp = line.split()
                    if temp[0] == "saveconc":
                        saveconc_files.append(lines[i + 1].strip())
            self.saveconc_files = saveconc_files

            results = {}
            for file in saveconc_files:
                with open(os.path.abspath(model_dir + '/' + file.lower())) as f:
                    lines = f.readlines()
                    # build header
                    aux = lines[5].split()
                    aux2 = [e + '/L' for e in aux[-2].split('/L') if e]
                    header = aux[:17] + aux2 + aux[-1:]
                    # obtain spreadsheet-form results
                    temp = []
                    for line in lines[6:]:
                        temp.append(line.split())
                    temp = pd.DataFrame(temp, columns=header)
                    results[file] = temp

            self.conc_output = results

    def simulated(self, outfile, variables='FLOWm^3/s', plot=True):
        if len(self.conc_output) == 0:
            print('You must execute method ''readOutput'' method first...')
        else:
            try:
                wat_file = self.conc_output[outfile]
            except ValueError:
                sys.exit('Error: The file {:s} was not found...'.format(outfile))

            try:
                if type(variables) == str:
                    variables = [variables]

                for variable in variables:
                    ts_sim = [float(value) for value in wat_file[variable]]
                    years = [int(year) for year in wat_file.iloc[:, 0]]
                    days = [int(day) for day in wat_file.iloc[:, 1]]
                    dates = [dt.datetime(year, 1, 1) + dt.timedelta(day - 1) for year, day in zip(years, days)]

                    ts_sim = pd.DataFrame(ts_sim, index=dates, columns=['Values'])

                    self.sim_ts[variable] = ts_sim

                    if plot:
                        ts_sim.plot(color='k', legend=False)
            except ValueError:
                sys.exit('Please, verify the names within the variables you want to plot '
                         '(names should be as they are presented in the output file)...')


def is_float(s):
    """Determine whether a string can be converted to a float number.
    input(s):
        s = string
    output(s):
        True/False
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_int(s):
    """Determine whether a string can be converted to an integer number.
    input(s):
        s = string
    output(s):
        True/False
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_all_subs_hrus(dir_list):
    hru_codes = [x.split('.')[0] for x in dir_list if (x.endswith('.{}'.format('hru')) and not x.startswith('output'))]

    sub_aux = [x[0:5] for x in hru_codes]
    hru_aux = [x[-4:] for x in hru_codes]

    sub_list = sorted(list(set(sub_aux)))
    hru_list = []

    for sub_i in sub_list:
        ind = [i for i, x in enumerate(sub_aux) if x == sub_i]
        temp = [hru_aux[i] for i in ind]
        hru_list.append(sorted(temp))

    sub_list = [int(x) for x in sub_list]
    hru_list = [[int(x) for x in y] for y in hru_list]

    return sub_list, hru_list


def sort_subs_hrus(subs, hrus, params):
    ind = [i for i, _ in sorted(enumerate(subs), key=lambda x: len(x[1]), reverse=True)]
    subs2 = [subs[i] for i in ind]
    hrus2 = [hrus[i] for i in ind]
    params2 = [params[i] for i in ind]

    return subs2, hrus2, params2


def simplify_subs(subs, hrus):  # UNUSED DEFINITION
    """Remove duplicated subbasins and related hrus from a list of subbasins.
    inputs(s):
        subs = list of lists of subbasins
        hrus = list of lists of hrus related to each subbasin
    output(s):
        subs2 = list of lists of subbasins with duplicates removed
        hrus2 = corresponding list of hrus for subs2
    """
    subs2 = list(subs)
    hrus2 = list(hrus)

    for i in range(0, len(subs) - 1):
        for j in range(i + 1, len(subs)):
            a = subs2[i]
            b = subs2[j]
            if len(a) <= len(b):
                c = [k for k in a if k in b]
                if len(c) > 0:
                    subs2[j] = [k for k in b if k not in c]
            else:
                c = [k for k in b if k in a]
                if len(c) > 0:
                    subs2[i] = [k for k in a if k not in c]

    for i, sub_i in enumerate(subs2):
        sub_i0 = subs[i]
        if len(sub_i) < len(sub_i0):
            hru_i = []
            hru_i0 = hrus[i]
            for element in sub_i:
                ind = [w for w, x in enumerate(sub_i0) if x == element]
                for w in ind:
                    hru_i.append(hru_i0[w])
            hrus2[i] = hru_i

    return subs2, hrus2


def replace_line(line, value, method, ext, num_format):
    """Replace old values in a line by new values.
    input(s):
        line = string containing values to be changed. Generally, it consists of two parts separated by either ':' or '|'
        value = number or string used to determine the new values in the given line
        method = indicates how the new values are going to be determined. So far, there are four options:
            'replace', the new value is the input 'value';
            'multiply', the new value is obtained by changing the original value in the line by a fraction given by the input 'value';
            'factor', the new value is the old value multiplied by the input 'value'
            'add', the new value is the old value plus the input 'value'
        ext = SWAT input file extension
        num_format = format of the number to be replaced
    output(s):
        new_line = string containing the original line with values modified according to the inputs 'value' and 'method'
    """

    if (ext == 'sol') | (ext == 'chm'):  # especial case for soil properties
        parts = line.split(':')
        num = parts[1].strip()
        new_value = []

        if is_float(num) | is_int(num):  # change single numeric values
            if method == 'replace':
                new_value = value
            elif method == 'multiply':
                new_value = (1 + value) * float(num)
            elif method == 'factor':
                new_value = value * float(num)
            elif method == 'add':
                new_value = value + float(num)

            part1 = '{:{}}'.format(new_value, num_format)

        else:  # change array of values
            if not num_format == 's':
                # get the number of positions from num_format
                n = int(num_format.split('.')[0])
                # split string of numbers based on N and convert to float
                nums = [float(num[i:i + n]) for i in range(0, len(num), n)]

                if method == 'replace':
                    new_value = [value for _ in nums]
                elif method == 'multiply':
                    new_value = [(1.0 + value) * x for x in nums]
                elif method == 'factor':
                    new_value = [value * x for x in nums]
                elif method == 'add':
                    new_value = [(value + x) for x in nums]

                part1 = ''.join(['{:{}}'.format(x, num_format) for x in new_value])

            else:  # change strings
                part1 = ' {:13s}'.format(value)

        new_line = '{part1}:{part2}\n'.format(part1=parts[0], part2=part1)

    elif ext == 'rte':  # especial case for routing (some parameters can be an array of values)
        parts = line.split('|')
        num = parts[0].strip()
        new_value = []

        if is_float(num) | is_int(num):  # change single numeric values
            if method == 'replace':
                new_value = value
            elif method == 'multiply':
                new_value = (1 + value) * float(num)
            elif method == 'factor':
                new_value = value * float(num)
            elif method == 'add':
                new_value = value + float(num)

            if is_int(num): # determine number format since num_format is not provided
                part0 = '{:14d}'.format(int(new_value))
            elif is_float(num):
                nd = 6 #abs(decimal.Decimal(num).as_tuple().exponent)
                part0 = '{:14.{}f}'.format(new_value, nd)
            else:
                part0 = '{:14s}'.format(value)

            new_line = '{part1}    |{part2}'.format(part1=part0, part2=parts[1])

        else:  # change array of values (format is provided)
            n = int(num_format.split('.')[0])  # get the number of positions from num_format
            nums = [float(num[i:i + n]) for i in
                    range(0, len(num), n)]  # split string of numbers based on N and convert to float

            if method == 'replace':
                new_value = [value for _ in nums]
            elif method == 'multiply':
                new_value = [(1.0 + value) * x for x in nums]
            elif method == 'factor':
                new_value = [value * x for x in nums]
            elif method == 'add':
                new_value = [(value + x) for x in nums]

            part0 = ''.join(['{:{}}'.format(x, num_format) for x in new_value])

            if len(parts) < 2:
                parts.append('\n')

            new_line = '{part1}|{part2}'.format(part1=part0, part2=parts[1])

    else:  # generic case (only single values; array of values requires an especial case)
        parts = line.split('|')
        num = parts[0].strip()
        new_value = []

        if method == 'replace':
            new_value = value
        elif method == 'multiply':
            new_value = (1 + value) * float(num)
        elif method == 'factor':
            new_value = value * float(num)
        elif method == 'add':
            new_value = value + float(num)

        if is_int(num):
            part0 = '{:16d}'.format(int(new_value))
        elif is_float(num):
            nd = 6 #abs(decimal.Decimal(num).as_tuple().exponent)
            part0 = '{:16.{}f}'.format(new_value, nd)
        else:
            part0 = '{:13s}   '.format(value)

        new_line = '{part1}    |{part2}'.format(part1=part0, part2=parts[1])

    return new_line


def write_ext_files(param_df, dir_list, subbasins, hrus, exts, input_dir, output_dir):
    # build reference lists of subbasin and hru codes
    sub_ref = ['{:05d}0000'.format(x) for x in subbasins]
    hru_ref = ['{:05d}{:04d}'.format(x, y) for i, x in enumerate(subbasins) for y in hrus[i]]

    for ext in exts:
        # get files in input directory with extension '.ext'
        files_all = [x for x in dir_list if (x.endswith('.{}'.format(ext)) and not x.startswith('output'))]
        param = param_df.loc[(param_df.ext == ext)].to_dict(orient='index')
        n_line = []
        txtformat = []

        if ext == 'sol':
            var_list = ['SNAM', 'HYDGRP', 'SOL_ZMX', 'ANION_EXCL', 'SOL_CRK', 'TEXTURE',
                       'SOL_Z', 'SOL_BD', 'SOL_AWC', 'SOL_K', 'SOL_CBN', 'SOL_CLAY', 'SOL_SILT',
                       'SOL_SAND', 'SOL_ROCK', 'SOL_ALB', 'USLE_K', 'SOL_EC', 'SOL_CAL', 'SOL_PH']
            n_line = [x for x in range(2, len(var_list) + 2)]
            txtformat = ['s', 's', '12.2f', '6.3f', '6.3f', 's',
                         '12.2f', '12.2f', '12.2f', '12.2f', '12.2f', '12.2f', '12.2f',
                         '12.2f', '12.2f', '12.2f', '12.2f', '12.2f', '12.2f', '12.2f']

        elif ext == 'chm':
            var_list = ['SOL_NO3', 'SOL_ORGN', 'SOL_SOLP', 'SOL_ORGP', 'PPERCO_SUB']
            n_line = [x for x in range(4, len(var_list) + 4)]
            txtformat = ['12.2f', '12.2f', '12.2f', '12.2f', '12.2f']

        else:
            var_list = []

        # create list of files to change
        files_all.sort()

        if len(files_all) > 1:
            crit = int(files_all[0][8])
            if crit == 0:
                files = ['{part1}.{part2}'.format(part1=x, part2=ext) for x in sub_ref]
            else:
                files = ['{part1}.{part2}'.format(part1=x, part2=ext) for x in hru_ref]
        else:
            files = files_all  # this is the case of .bsn and .wwq

        # modify list of files
        for file in files:
            with open(os.path.abspath(input_dir + '/' + file), 'r', encoding='ISO-8859-1') as f:
                data = f.readlines()
                param_names = list(param.keys())
                for param_name in param_names:
                    c = 0
                    if var_list:
                        ind = [i for i, x in enumerate(var_list) if param_name in x][0]
                        c = n_line[ind] - 1
                        num_format = txtformat[ind]
                    else:

                        if (ext == 'rte') & (param_name == 'CH_ERODMO'):
                            c = 23
                            num_format = '6.2f'
                        elif (ext == 'rte') & (param_name == 'HRU_SALT'):
                            c = 28
                            num_format = '6.2f'
                        else:
                            for line in data:
                                if re.search(param_name, line):
                                    num_format = ''
                                    break
                                else:
                                    c = c + 1

                    new_line = replace_line(data[c],
                                           param[param_name]['value'],
                                           param[param_name]['method'],
                                           param[param_name]['ext'],
                                           num_format)
                    data[c] = new_line

            with open(os.path.abspath(output_dir + '/' + file), "w") as f:
                f.writelines(data)


def prepare_subs_hrus(dir_list, subs, hrus):
    # get all subbasins and hrus
    sub_list, hru_list = get_all_subs_hrus(dir_list)
    # complete subs and hrus lists if necessary
    if len(subs) == 0:
        subs = [list(sub_list)]
        hrus = [list(hru_list)]
    else:
        subs = [sub_i if len(sub_i) > 0 else list(sub_list) for sub_i in subs]
        if len(hrus) == 0:
            hrus = list(hrus)
            for sub_i in subs:
                hru_i = []
                for element in sub_i:
                    ind = [i for i, x in enumerate(sub_list) if x == element]
                    for i in ind:
                        hru_i.append(hru_list[i])
                hrus.append(hru_i)
        else:
            for i, sub_i in enumerate(subs):
                hru_i = hrus[i]
                if len(hru_i) == 0:
                    for element in sub_i:
                        ind = [w for w, x in enumerate(sub_list) if x == element]
                        for w in ind:
                            hru_i.append(hru_list[w])
                    hrus[i] = hru_i

                    # subs, hrus = simplify_subs(subs,hrus) # uncomment if you want to remove duplicated subbasins and respective hrus

    return subs, hrus


def write_new_files(param_all, subs, hrus, input_dir, output_dir):
    """Write new SWAT text input files based on a list of parameters to be changed.
    input(s):
        param_all = dictionary containing a set of n parameters to be modified. The format is as follows:
            {'key_name_1':[list_1], ..., 'key_name_n':[list_n]}
            where:
                'keyname_i' is the ith string with the SWAT name of the ith parameter to be modified
                [list_i] is the ith list of three elements defining the inputs 'value','method', and 'ext' of the function 'replaceLine'
        input_dir = directory for the SWAT TxtInOut folder to be modified
        output_dir = output directory where the new text files will be written
    """
    if type(param_all) == dict:
        param_all = [param_all]

    if (len(param_all) < len(subs)) and (len(param_all) != 0 and len(param_all) != 1):
        sys.exit('List of parameters must have the same length as subbasins list or be empty')

    dir_list = os.listdir(input_dir)

    # get list of subbasins and hrus ready to write files
    subs, hrus = prepare_subs_hrus(dir_list, subs, hrus)
    param_list = []

    if len(param_all) < len(subs):
        if len(param_all) == 1:
            param_list = [param_all[0] for _ in subs]
        elif len(param_all) == 0:
            param_list = [[] for _ in subs]
    elif (len(param_all) > 1) and (len(param_all) > len(subs)):
        param_list = [param_all[i] for i, _ in enumerate(subs)]
    else:
        param_list = list(param_all)

    # sort param, subs and hrus according to size of subs
    subs, hrus, param_list = sort_subs_hrus(subs, hrus, param_list)

    for i, sub in enumerate(subs):
        hru = hrus[i]
        param = param_list[i]
        param_df = pd.DataFrame.from_dict(param,
                                          orient='index',
                                          columns=['value', 'method', 'ext'])
        exts = param_df.ext.unique().tolist()
        write_ext_files(param_df, dir_list, sub, hru, exts, input_dir, output_dir)


def write_mgt_tables(mgt_tables, ext, subs, hrus, output_dir):

    if type(mgt_tables) is not list:
        mgt_tables = [mgt_tables]

    if (len(mgt_tables) < len(subs)) and (len(mgt_tables) != 0 and len(mgt_tables) != 1):
        sys.exit('List of management tables must have the same length as subbasins list or be empty')

    dir_list = os.listdir(output_dir)

    # get list of subbasins and hrus ready to write files
    subs, hrus = prepare_subs_hrus(dir_list, subs, hrus)

    if len(mgt_tables) < len(subs):
        if len(mgt_tables) == 1:
            mgt_list = [mgt_tables[0] for _ in subs]
        elif len(mgt_tables) == 0:
            mgt_list = [[] for _ in subs]
    elif (len(mgt_tables) > 1) and (len(mgt_tables) > len(subs)):
        mgt_list = [mgt_tables[i] for i, _ in enumerate(subs)]
    else:
        mgt_list = list(mgt_tables)

    # sort param, subs and hrus according to size of subs
    subs, hrus, mgt_list = sort_subs_hrus(subs, hrus, mgt_list)

    for i, sub in enumerate(subs):
        hru = hrus[i]
        mgt_table = mgt_list[i]
        if type(mgt_table) == list and len(mgt_table) == 0:
            break
        else:
            write_mgt_files(mgt_table, ext, dir_list, sub, hru, output_dir)


# set of definitions for tracking back draining subbasins using .fig files

def is_empty_list(inflow, ref):
    aux = ref.loc[pd.to_numeric(ref[2]) == inflow]
    if len(aux) > 0:
        iflag = bool(0)
    else:
        iflag = bool(1)

    return iflag


def type_element(inflow, subbasin_list, route_list, add_list):
    temp = subbasin_list
    if is_empty_list(inflow, temp):
        temp = route_list
        if is_empty_list(inflow, temp):
            temp = add_list
            if is_empty_list(inflow, temp):
                sys.exit('Verify your .fig file since an element was not found...')
            else:
                type_out = 'add'
        else:
            type_out = 'route'
    else:
        type_out = 'subbasin'

    return type_out


def subbasin_root(inflow, subbasin_list, route_list, add_list, subbasins):
    type_inflow = type_element(inflow, subbasin_list, route_list, add_list)
    if type_inflow == 'add':
        inflow_elements = add_list.loc[pd.to_numeric(add_list[2]) == inflow].to_numpy()[:, -2:].reshape(-1, 1)
        for element in inflow_elements:
            subbasins = subbasin_root(int(element), subbasin_list, route_list, add_list, subbasins)
    elif type_inflow == 'route':
        element = route_list.loc[pd.to_numeric(route_list[2]) == inflow].to_numpy()[:, -1:]
        subbasins = subbasin_root(int(element), subbasin_list, route_list, add_list, subbasins)
    elif type_inflow == 'subbasin':
        subbasins.append(inflow)

    subbasins.sort()
    return subbasins


def update_txt(input_dir):

    file = 'basins.bsn'
    lines = {'r2adj': [1, ' R2ADJ: Curve number retention parameter adjust factor\n', '{:16d}'],
             'sstmaxd_bsn': [0.0, ' SSTMAXD_BSN: Static maximum depressional storage (mm)\n', '{:16.2f}'],
             'ismax': [1, ' ISMAX: 1 = Dynamic sstmaxd as a function of random roughness 0 = read from .bsn or .sdr\n', '{:16d}'],
             'iroutunit': [0, ' IROUTUNIT: \n', '{:16d}']}

    new_lines = ['{part1}    |{part2}'.format(part1=x[2].format(x[0]), part2=x[1]) for _, x in lines.items()]

    with open(os.path.abspath(input_dir + '/' + file), 'r', encoding='ISO-8859-1') as f:
        data = f.readlines() + new_lines

    with open(os.path.abspath(input_dir + '/' + file), "w") as f:
        f.writelines(data)

    dir_list = os.listdir(input_dir)
    files = [x for x in dir_list if (x.endswith('.{}'.format('sdr')) and not x.startswith('output'))]
    lines = {'sstmaxd': [0.0, ' SSTMAXD: Static maximum depressional storage (mm)\n', '{:10.2f}']}
    new_lines = ['{part1}    |{part2}'.format(part1=x[2].format(x[0]), part2=x[1]) for _, x in lines.items()]

    for file in files:

        with open(os.path.abspath(input_dir + '/' + file), 'r', encoding='ISO-8859-1') as f:
            data = f.readlines() + new_lines

        with open(os.path.abspath(input_dir + '/' + file), "w") as f:
            f.writelines(data)


def write_mgt_files(mgt_table, ext, dir_list, subbasins, hrus, output_dir):

    # build reference lists of hru codes
    hru_ref = ['{:05d}{:04d}'.format(x, y) for i, x in enumerate(subbasins) for y in hrus[i]]

    # create list of files to change
    files = ['{part1}.{part2}'.format(part1=x, part2=ext) for x in hru_ref]

    if ext == 'mgt':
        table = mgt_table.op_sched
        ind = 30
    elif ext == 'ops':
        table = mgt_table.schedule
        ind = 1

    table = [x + '\n' for x in table]

    for file in files:
        with open(os.path.abspath(output_dir + '/' + file), 'r', encoding='ISO-8859-1') as f:
            data = f.readlines()
            del data[ind:]
            data = data + table

        with open(os.path.abspath(output_dir + '/' + file), "w") as f:
            f.writelines(data)


# set of functions for obtaining pre-defined model evaluations

def run_single_model(swat_model, out_file, out_var, lock, runid):
    # assigning ID to swat output folder
    if not isinstance(lock, int):
        lock.acquire()
    ind = runid.get_id()
    swat_model.new_model_name = 'RUN{:04d}'.format(ind)
    # prepare SWAT run
    swat_model.prepare_swat()
    if not isinstance(lock, int):
        lock.release()
    # execute model
    swat_model.run_swat()
    # get output time series of variable of interest
    swat_model.read_output()
    swat_model.simulated(out_file, plot=False)
    swat_model.remove_swat()
    # get simulated time series of interest
    simulated = swat_model.sim_ts[out_var]
    # remove Run ID
    # runid.remove_id(ind)
    return simulated
