#!/usr/bin/env python3
"""
Write .ops files from scratch and modify .sub files accordingly given an existing SWAT model
@author: J. Sebastian Hernandez-Suarez
"""

import os


def create_ops(input_dir, output_dir, start_date):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    dir_list = os.listdir(input_dir)

    sub_files = [file for file in dir_list if file.endswith('.sub')]  # get files in input directory with extension '.sub'
    ops_files = [file.split('.')[0] + '.ops' for file in dir_list if
                file.endswith('.hru')]  # change '.hru' extension to '.ops'

    # add .ops file into last lines of .sub files

    for subfile in sub_files:
        with open(input_dir + '/' + subfile, 'r', encoding='ISO-8859-1') as file:
            data = file.readlines()
            nlines = len(data)
            newlines = list()
            for i in range(61, nlines):
                old_line = data[i]
                ops_part = old_line[0:13].split('.')[0] + '.ops'
                aux = list(old_line)
                aux[65:78] = list(ops_part)
                newlines.append("".join(aux))
            data[61:nlines] = newlines

        with open(output_dir + '/' + subfile, "w") as file:
            file.writelines(data)

    # create .ops files (generic BMP)
    for opsFile in ops_files:
        with open(input_dir + '/' + opsFile.split('.')[0] + '.hru', 'r') as file:
            first_line = file.readline()

        with open(output_dir + '/' + opsFile, "w") as file:
            month = start_date['MONTH']
            day = start_date['DAY']
            iyear = start_date['IYEAR']
            mgt_op = 10
            ro_nmp_flag = 0
            ro_nmp_sed = 0
            ro_nmp_pp = 0
            ro_nmp_sp = 0
            ro_nmp_pn = 0
            ro_nmp_sn = 0
            ro_nmp_bac = 0

            line = ' %2d %2d     %4d %2d %4d     %6.2f %12.5f %6.2f             %8.2f %6.2f %5.2f' % \
                   (month, day, iyear, mgt_op, ro_nmp_flag, ro_nmp_sed, ro_nmp_pp, ro_nmp_sp, ro_nmp_pn, ro_nmp_sn,
                    ro_nmp_bac)

            data = list((first_line, line))
            file.writelines(data)


def link_sdr(input_dir, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    dir_list = os.listdir(input_dir)

    sub_files = [file for file in dir_list if file.endswith('.sub')]  # get files in input directory with extension '.sub'

    # add .sdr file into last line of .sub files

    for subfile in sub_files:
        with open(input_dir + '/' + subfile, 'r', encoding='ISO-8859-1') as file:
            data = file.readlines()
            nlines = len(data)
            newlines = list()
            for i in range(61, nlines):
                old_line = data[i]
                sdr_part = old_line[0:13].split('.')[0] + '.sdr'
                aux = list(old_line)
                aux[91:105] = list(sdr_part)
                newlines.append("".join(aux))
            data[61:nlines] = newlines

        with open(output_dir + '/' + subfile, "w") as file:
            file.writelines(data)
