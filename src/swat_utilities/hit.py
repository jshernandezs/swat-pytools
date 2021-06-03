#!/usr/bin/env python3
"""
Class for computing hydrologic indices
"""
import numpy as np
import pandas as pd
import scipy.special as sp
from scipy import linalg
import warnings
import collections.abc as collections
import bisect
import re
from swat_utilities.performance import read_obs_file


# noinspection DuplicatedCode
class HydrologicIndex:
    def __init__(self, series, drainage_area=1, wy_month=10, obs_file=None):

        if obs_file is not None:
            obs = read_obs_file(obs_file)
            series = series[obs.index[0]:obs.index[-1]]

        dates = series.index
        year = dates.year
        month = dates.month
        day = np.array(dates.dayofyear.tolist())
        daym = np.array(dates.day.tolist())

        self.drainage_area = drainage_area

        self.day = day
        self.daym = daym
        self.month = month
        self.w_year = year + (month >= wy_month) if wy_month > 1 else year
        self.w_year_list = self.w_year.unique()
        self.discharge = series.values
        self.wy_month = wy_month

        self.mean_discharge = np.array([])
        self.median_discharge = np.array([])
        self.median_log10_discharge = np.array([])

        self.mean_by_w_year = np.array([])
        self.median_by_w_year = np.array([])
        self.min_by_w_year = np.array([])
        self.max_by_w_year = np.array([])
        self.std_by_w_year = np.array([])
        self.cv_by_w_year = np.array([])

        self.std_month_year = np.array([])
        self.mean_month_year = np.array([])
        self.median_month_year = np.array([])
        self.max_month_year = np.array([])
        self.min_month_year = np.array([])
        self.mean_by_month = np.array([])

        self.rolling3day_mean = np.array([])
        self.rolling7day_mean = np.array([])
        self.rolling30day_mean = np.array([])
        self.rolling90day_mean = np.array([])

        self.discharge_percentiles_list = np.array([])
        self.discharge_percentiles = np.array([])
        self.discharge_log10_percentiles = np.array([])
        self.discharge_exceedance_percentiles_list = np.array([])
        self.discharge_exceedance_percentiles = np.array([])

        self.peak_thresh_prct = np.array([])
        self.peak_thresh = np.array([])
        self.base_flow_index = np.array([])
        self.colwell_bounds = np.array([])
        self.colwell_mat = np.array([])

        list_dh = ['dh{:d}'.format(x) for x in range(1, 25)]
        list_dl = ['dl{:d}'.format(x) for x in range(1, 21)]
        list_fh = ['fh{:d}'.format(x) for x in range(1, 12)]
        list_fl = ['fl{:d}'.format(x) for x in range(1, 4)]
        list_ma = ['ma{:d}'.format(x) for x in range(1, 46)]
        list_mh = ['mh{:d}'.format(x) for x in range(1, 28)]
        list_ml = ['ml{:d}'.format(x) for x in range(1, 23)]
        list_ta = ['ta{:d}'.format(x) for x in range(1, 4)]
        list_th = ['th{:d}'.format(x) for x in range(1, 4)]
        list_tl = ['tl{:d}'.format(x) for x in range(1, 5)]
        list_ra = ['ra{:d}'.format(x) for x in range(1, 10)]

        self.list_all = list(flatten([list_dh, list_dl, list_fh, list_fl, list_ma, list_mh,
                                      list_ml, list_ta, list_th, list_tl, list_ra]))
        
        self.list_mag_seven = ['mag{:d}'.format(x) for x in range(1, 8)]

        list_dh = ['dh{:d}'.format(x) for x in list(flatten([np.arange(1, 15), np.arange(17, 25)]))]
        list_dl = ['dl{:d}'.format(x) for x in list(flatten([np.arange(1, 16), [20]]))]
        list_fh = ['fh{:d}'.format(x) for x in list(flatten([np.arange(3, 12)]))]
        list_fl = ['fl{:d}'.format(x) for x in list(flatten([[3]]))]
        list_ma = ['ma{:d}'.format(x) for x in list(flatten([np.arange(1, 4)]))]
        list_mh = ['mh{:d}'.format(x) for x in list(flatten([np.arange(13, 15), np.arange(18, 28)]))]
        list_ml = ['ml{:d}'.format(x) for x in list(flatten([[13], np.arange(17, 23)]))]
        list_ta = ['ta{:d}'.format(x) for x in list(flatten([[3]]))]
        list_th = ['th{:d}'.format(x) for x in list(flatten([[3]]))]
        list_tl = ['tl{:d}'.format(x) for x in list(flatten([np.arange(3, 5)]))]

        self.list_single = list(flatten([list_dh, list_dl, list_fh, list_fl, list_ma, list_mh,
                                         list_ml, list_ta, list_th, list_tl]))

        self.list_multiple = ['dh15_16', 'dl16_17', 'dl18_19', 'fh1_2', 'fl1_2',
                              'ma4_11', 'ma12_23', 'ma24_35', 'ma36_40', 'ma41_45',
                              'mh1_12', 'mh15_17', 'ml1_12', 'ml14_16', 'ra1_9',
                              'ta1_2', 'th1_2', 'tl1_2']

    def get_index(self, name, opt=False):

        with np.errstate(invalid='ignore', divide='ignore'):

            name = name.lower()

            if name in self.list_single:
                indices = None
                names = None
                index = getattr(self, 'get_{:s}'.format(name))()
            else:
                if name[:3] == 'mag':
                    names = ['mag' + str(i) for i in range(1, 8)]
                    number = name[3:]
                    iflag = 0
                    if int(number) in range(1, 8):
                        ind = int(number) - 1
                        indices = self.get_mag_seven()
                        index = indices[ind]
                        iflag = 1
                else:
                    group = name[:2]
                    number = name[2:]
                    iflag = 0
                    for item in self.list_multiple:
                        item_group = item[:2]
                        if group == item_group:
                            r = item[2:].split('_')
                            ran = range(int(r[0]), int(r[1]) + 1)
                            names = [item[:2] + str(i) for i in ran]
                            if int(number) in ran:
                                ind = int(number) - ran[0]
                                indices = getattr(self, 'get_{:s}'.format(item))()
                                index = indices[ind]
                                iflag = 1
                                break
                if iflag == 0:
                    index = None
                    print('Check the index name you provided')

            if opt:
                # report all resulting indices for avoiding repeated evaluations when using get_all_indices()
                return index, indices, names
            else:
                return index

    def get_list_indices(self, list_hi):
        
        # list_hi.sort()
        indices = {}
        groups = ['mag', 'dh', 'dl', 'fl', 'fh', 'ma', 'mh', 'ml', 'ra', 'ta', 'th', 'tl']
        for group in groups:
            subset = [hi for hi in list_hi if " ".join(re.findall("[a-zA-Z]+", hi)) == group]
            aux = []
            names = []
            for index in subset:
                if index not in names:
                    ind_value = self.get_index(index, opt=True)
                    if ind_value[1] is None:
                        indices[index.upper()] = ind_value[0]
                    else:
                        aux = np.array(ind_value[1])
                        names = np.array(ind_value[2])
                        indices[index.upper()] = aux[names == index.lower()][0]
                else:
                    indices[index.upper()] = aux[np.array(names) == index.lower()][0]

        return indices

    def get_all_hit_indices(self):

        hit = {}
        groups = ['dh', 'dl', 'fl', 'fh', 'ma', 'mh', 'ml', 'ra', 'ta', 'th', 'tl']
        root = np.array([x[:2] for x in self.list_all])
        list_all = np.array(self.list_all)
        for group in groups:
            subset = list_all[root == group]
            temp = []
            n = 0
            for index in subset:
                if n == 0:
                    ind_value = self.get_index(index, opt=True)
                    if ind_value[1] is None:
                        temp.append(ind_value[0])
                    else:
                        temp.append(ind_value[1])
                        n = len(ind_value[1]) - 1
                else:
                    n = n - 1
                    continue

            hit[group.upper()] = list(flatten(temp))
            
        return hit

    def get_mag_seven(self):

        discharge = self.discharge

        # compute L-moment ratios for discharge time-series
        complmom = lmom_ratios(discharge, nmom=4)
        complmom[1] = complmom[1] / complmom[0]
        # compute AR(1) correlation coef.
        ar1_cor = self.get_ar1()
        # compute amplitude and phase
        seasonality = self.get_seasonality()

        mag_seven = list(flatten([complmom, ar1_cor, seasonality]))

        return mag_seven

    def get_dh1(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        dh1 = np.nanmean(self.max_by_w_year)

        return dh1

    def get_dh2(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling3day_mean.size == 0:
            self.rolling3day_mean = np.array(np.convolve(discharge[:, 0], np.ones((3,)) / 3, 'same'))
            self.rolling3day_mean[:1] = np.nan
            self.rolling3day_mean[-1:] = np.nan

        max_rolling3day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling3day_mean[w_year == wy]) for wy in w_year_list])
        dh2 = np.nanmean(max_rolling3day_mean_by_w_year)

        return dh2

    def get_dh3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        max_rolling7day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dh3 = np.nanmean(max_rolling7day_mean_by_w_year)

        return dh3

    def get_dh4(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        max_rolling30day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dh4 = np.nanmean(max_rolling30day_mean_by_w_year)

        return dh4

    def get_dh5(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling90day_mean.size == 0:
            self.rolling90day_mean = np.array(np.convolve(discharge[:, 0], np.ones((90,)) / 90, 'same'))
            self.rolling90day_mean[:44] = np.nan
            self.rolling90day_mean[-44:] = np.nan

        max_rolling90day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling90day_mean[w_year == wy]) for wy in w_year_list])
        dh5 = np.nanmean(max_rolling90day_mean_by_w_year)

        return dh5

    def get_dh6(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        dh6 = np.nanstd(self.max_by_w_year, ddof=1) / np.nanmean(self.max_by_w_year) * 100

        return dh6

    def get_dh7(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling3day_mean.size == 0:
            self.rolling3day_mean = np.array(np.convolve(discharge[:, 0], np.ones((3,)) / 3, 'same'))
            self.rolling3day_mean[:1] = np.nan
            self.rolling3day_mean[-1:] = np.nan

        max_rolling3day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling3day_mean[w_year == wy]) for wy in w_year_list])
        dh7 = np.nanstd(max_rolling3day_mean_by_w_year, ddof=1) / np.nanmean(max_rolling3day_mean_by_w_year) * 100

        return dh7

    def get_dh8(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        max_rolling7day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dh8 = np.nanstd(max_rolling7day_mean_by_w_year, ddof=1) / np.nanmean(max_rolling7day_mean_by_w_year) * 100

        return dh8

    def get_dh9(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        max_rolling30day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dh9 = np.nanstd(max_rolling30day_mean_by_w_year, ddof=1) / np.nanmean(max_rolling30day_mean_by_w_year) * 100

        return dh9

    def get_dh10(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling90day_mean.size == 0:
            self.rolling90day_mean = np.array(np.convolve(discharge[:, 0], np.ones((90,)) / 90, 'same'))
            self.rolling90day_mean[:44] = np.nan
            self.rolling90day_mean[-44:] = np.nan

        max_rolling90day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling90day_mean[w_year == wy]) for wy in w_year_list])
        dh10 = np.nanstd(max_rolling90day_mean_by_w_year, ddof=1) / np.nanmean(max_rolling90day_mean_by_w_year) * 100

        return dh10

    def get_dh11(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        dh11 = np.nanmean(self.max_by_w_year) / self.median_discharge

        return dh11

    def get_dh12(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        max_rolling7day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dh12 = np.nanmean(max_rolling7day_mean_by_w_year) / self.median_discharge

        return dh12

    def get_dh13(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        max_rolling30day_mean_by_w_year = np.array(
            [np.nanmax(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dh13 = np.nanmean(max_rolling30day_mean_by_w_year) / self.median_discharge

        return dh13

    def get_dh14(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_month_year.size == 0:
            self.mean_month_year = get_fun_by_month_year('mean', discharge, w_year, w_year_list, month)

        dh14 = get_percentile(self.mean_month_year, 0.95) / np.nanmean(self.mean_month_year)

        return dh14

    def get_dh15_16(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 75]

        event_mask_per_year = [discharge[w_year == wy] > thresh for wy in w_year_list]
        total_events_duration_per_year = np.array([np.sum(x) for x in event_mask_per_year])
        n_events_per_year = [
            np.sum(np.diff(np.append(0, event_mask_per_year[i])) > 0) for i, _ in enumerate(w_year_list)]

        avg_event_duration_per_year = total_events_duration_per_year / np.array(n_events_per_year)
        avg_event_duration_per_year[np.array(n_events_per_year) == 0] = 0

        dh15 = np.nanmedian(avg_event_duration_per_year)
        dh16 = np.nanstd(avg_event_duration_per_year, ddof=1) / np.nanmean(avg_event_duration_per_year) * 100

        return dh15, dh16

    def get_dh17(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = self.median_discharge

        event_mask = discharge > thresh
        total_events_duration = np.sum(event_mask)
        n_events = np.sum(np.diff(np.append(0, event_mask)) > 0)

        dh17 = total_events_duration / n_events
        
        if n_events == 0:
           dh17 = 0

        return dh17

    def get_dh18(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 3 * self.median_discharge

        event_mask = discharge > thresh
        total_events_duration = np.sum(event_mask)
        n_events = np.sum(np.diff(np.append(0, event_mask)) > 0)

        dh18 = total_events_duration / n_events
        
        if n_events == 0:
           dh18 = 0

        return dh18

    def get_dh19(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 7 * self.median_discharge

        event_mask = discharge > thresh
        total_events_duration = np.sum(event_mask)
        n_events = np.sum(np.diff(np.append(0, event_mask)) > 0)

        dh19 = total_events_duration / n_events
        
        if n_events == 0:
            dh19 = 0

        return dh19

    def get_dh20(self):

        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 75]

        event_mask = discharge > thresh
        total_events_duration = np.sum(event_mask)
        n_events = np.sum(np.diff(np.append(0, event_mask)) > 0)

        dh20 = total_events_duration / n_events
        
        if n_events == 0:
            dh20 = 0

        return dh20

    def get_dh21(self):

        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 25]

        event_mask = discharge > thresh
        total_events_duration = np.sum(event_mask)
        n_events = np.sum(np.diff(np.append(0, event_mask)) > 0)

        dh21 = total_events_duration / n_events
        
        if n_events == 0:
            dh21 = 0

        return dh21

    def get_dh22(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        event_mask_per_year = [discharge[w_year == wy] < thresh for wy in w_year_list]

        median_event_duration_by_year = np.zeros((len(w_year_list), 1))
        for i in range(0, len(w_year_list)):
            event_start = np.array(np.nonzero(np.diff(np.append(0, event_mask_per_year[i])) > 0))
            event_end = np.array(np.nonzero(np.diff(np.append(event_mask_per_year[i], 0)) < 0))
            events_duration = event_end - event_start + 1
            median_event_duration_by_year[i] = np.nanmedian(events_duration)

        dh22 = np.nanmean(median_event_duration_by_year)

        return dh22

    def get_dh23(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        total_events_duration_per_year = np.array([np.sum(discharge[w_year == wy] > thresh) for wy in w_year_list])

        dh23 = np.nanmean(total_events_duration_per_year)

        return dh23

    def get_dh24(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        event_mask_per_year = [discharge[w_year == wy] < thresh for wy in w_year_list]

        max_event_duration_by_year = np.zeros((len(w_year_list), 1))
        for i in range(0, len(w_year_list)):
            event_start = np.array(np.nonzero(np.diff(np.append(0, event_mask_per_year[i])) > 0))
            event_end = np.array(np.nonzero(np.diff(np.append(event_mask_per_year[i], 0)) < 0))
            events_duration = event_end - event_start + 1

            if events_duration.size == 0:
                max_event_duration_by_year[i] = 0
            else:
                max_event_duration_by_year[i] = np.nanmax(events_duration)

        dh24 = np.nanmean(max_event_duration_by_year)

        return dh24

    def get_dl1(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        dl1 = np.nanmean(self.min_by_w_year)

        return dl1

    def get_dl2(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling3day_mean.size == 0:
            self.rolling3day_mean = np.array(np.convolve(discharge[:, 0], np.ones((3,)) / 3, 'same'))
            self.rolling3day_mean[:1] = np.nan
            self.rolling3day_mean[-1:] = np.nan

        min_rolling3day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling3day_mean[w_year == wy]) for wy in w_year_list])
        dl2 = np.nanmean(min_rolling3day_mean_by_w_year)

        return dl2

    def get_dl3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        min_rolling7day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dl3 = np.nanmean(min_rolling7day_mean_by_w_year)

        return dl3

    def get_dl4(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        min_rolling30day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dl4 = np.nanmean(min_rolling30day_mean_by_w_year)

        return dl4

    def get_dl5(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling90day_mean.size == 0:
            self.rolling90day_mean = np.array(np.convolve(discharge[:, 0], np.ones((90,)) / 90, 'same'))
            self.rolling90day_mean[:44] = np.nan
            self.rolling90day_mean[-44:] = np.nan

        min_rolling90day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling90day_mean[w_year == wy]) for wy in w_year_list])
        dl5 = np.nanmean(min_rolling90day_mean_by_w_year)

        return dl5

    def get_dl6(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        dl6 = np.nanstd(self.min_by_w_year, ddof=1) / np.nanmean(self.min_by_w_year) * 100
        return dl6

    def get_dl7(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling3day_mean.size == 0:
            self.rolling3day_mean = np.array(np.convolve(discharge[:, 0], np.ones((3,)) / 3, 'same'))
            self.rolling3day_mean[:1] = np.nan
            self.rolling3day_mean[-1:] = np.nan

        min_rolling3day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling3day_mean[w_year == wy]) for wy in w_year_list])
        dl7 = np.nanstd(min_rolling3day_mean_by_w_year, ddof=1) / np.nanmean(min_rolling3day_mean_by_w_year) * 100

        return dl7

    def get_dl8(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        min_rolling7day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dl8 = np.nanstd(min_rolling7day_mean_by_w_year, ddof=1) / np.nanmean(min_rolling7day_mean_by_w_year) * 100

        return dl8

    def get_dl9(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        min_rolling30day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dl9 = np.nanstd(min_rolling30day_mean_by_w_year, ddof=1) / np.nanmean(min_rolling30day_mean_by_w_year) * 100

        return dl9

    def get_dl10(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling90day_mean.size == 0:
            self.rolling90day_mean = np.array(np.convolve(discharge[:, 0], np.ones((90,)) / 90, 'same'))
            self.rolling90day_mean[:44] = np.nan
            self.rolling90day_mean[-44:] = np.nan

        min_rolling90day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling90day_mean[w_year == wy]) for wy in w_year_list])
        dl10 = np.nanstd(min_rolling90day_mean_by_w_year, ddof=1) / np.nanmean(min_rolling90day_mean_by_w_year) * 100

        return dl10

    def get_dl11(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        dl11 = np.nanmean(self.min_by_w_year) / self.median_discharge

        return dl11

    def get_dl12(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        min_rolling7day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling7day_mean[w_year == wy]) for wy in w_year_list])
        dl12 = np.nanmean(min_rolling7day_mean_by_w_year) / self.median_discharge

        return dl12

    def get_dl13(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling30day_mean.size == 0:
            self.rolling30day_mean = np.array(np.convolve(discharge[:, 0], np.ones((30,)) / 30, 'same'))
            self.rolling30day_mean[:14] = np.nan
            self.rolling30day_mean[-14:] = np.nan

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        min_rolling30day_mean_by_w_year = np.array(
            [np.nanmin(self.rolling30day_mean[w_year == wy]) for wy in w_year_list])
        dl13 = np.nanmean(min_rolling30day_mean_by_w_year) / self.median_discharge

        return dl13

    def get_dl14(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        dl14 = self.discharge_exceedance_percentiles[
                   self.discharge_exceedance_percentiles_list == 75] / self.median_discharge

        return dl14[0]

    def get_dl15(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        dl15 = self.discharge_exceedance_percentiles[
                   self.discharge_exceedance_percentiles_list == 90] / self.median_discharge

        return dl15[0]

    def get_dl16_17(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 25]

        event_mask_per_year = [discharge[w_year == wy] < thresh for wy in w_year_list]
        total_events_duration_per_year = np.array([np.sum(x) for x in event_mask_per_year])
        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] < thresh, int))) > 0) for wy in w_year_list]

        avg_event_duration_per_year = total_events_duration_per_year / np.array(n_events_per_year)
        avg_event_duration_per_year[np.array(n_events_per_year) == 0] = 0

        dl16 = np.nanmedian(avg_event_duration_per_year)
        dl17 = np.nanstd(avg_event_duration_per_year, ddof=1) / np.nanmean(avg_event_duration_per_year) * 100

        return dl16, dl17

    def get_dl18_19(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        zero_flow_days_per_year = np.array([np.sum(discharge[w_year == wy] == 0) for wy in w_year_list])

        dl18 = np.nanmean(zero_flow_days_per_year)
        dl19 = np.nanstd(zero_flow_days_per_year, ddof=1) / np.nanmean(zero_flow_days_per_year) * 100

        if np.isnan(dl19):
            dl19 = 0

        return dl18, dl19

    def get_dl20(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_month_year.size == 0:
            self.mean_month_year = get_fun_by_month_year('mean', discharge, w_year, w_year_list, month)

        dl20 = np.sum(self.mean_month_year.reshape(-1, 1) == 0)

        return dl20

    def get_fh1_2(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 75]

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh1 = np.nanmean(n_events_per_year)
        fh2 = np.nanstd(n_events_per_year, ddof=1) / np.nanmean(n_events_per_year) * 100

        return fh1, fh2

    def get_fh3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 3 * self.median_discharge

        n_events_per_year = [np.sum(discharge[w_year == wy] > thresh) for wy in w_year_list]

        fh3 = np.nanmean(n_events_per_year)

        return fh3

    def get_fh4(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 7 * self.median_discharge

        n_events_per_year = [np.sum(discharge[w_year == wy] > thresh) for wy in w_year_list]

        fh4 = np.nanmean(n_events_per_year)

        return fh4

    def get_fh5(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = self.median_discharge

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh5 = np.nanmean(n_events_per_year)

        return fh5

    def get_fh6(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 3 * self.median_discharge

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh6 = np.nanmean(n_events_per_year)

        return fh6

    def get_fh7(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        thresh = 7 * self.median_discharge

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh7 = np.nanmean(n_events_per_year)

        return fh7

    def get_fh8(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        thresh = self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 25]

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh8 = np.nanmean(n_events_per_year)

        return fh8

    def get_fh9(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        thresh = self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 75]

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh9 = np.nanmean(n_events_per_year)

        return fh9

    def get_fh10(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        thresh = np.nanmedian(self.min_by_w_year)

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh10 = np.nanmean(n_events_per_year)

        return fh10

    def get_fh11(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] > thresh, int))) > 0) for wy in w_year_list]

        fh11 = np.nanmean(n_events_per_year)

        return fh11

    def get_fl1_2(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        thresh = self.discharge_percentiles[self.discharge_percentiles_list == 25]

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] < thresh, int))) > 0) for wy in w_year_list]

        fl1 = np.nanmean(n_events_per_year)
        fl2 = np.nanstd(n_events_per_year, ddof=1) / np.nanmean(n_events_per_year) * 100

        return fl1, fl2

    def get_fl3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_discharge.size == 0:
            self.mean_discharge = np.nanmean(discharge)

        thresh = 0.05 * self.mean_discharge

        n_events_per_year = [
            np.sum(np.diff(np.append(0, np.array(discharge[w_year == wy] < thresh, int))) > 0) for wy in w_year_list]

        fl3 = np.nanmean(n_events_per_year)

        return fl3

    def get_ma1(self):

        discharge = self.discharge

        if self.mean_discharge.size == 0:
            self.mean_discharge = np.nanmean(discharge)

        ma1 = self.mean_discharge

        return ma1

    def get_ma2(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        ma2 = self.median_discharge

        return ma2

    def get_ma3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_by_w_year.size == 0:
            self.mean_by_w_year = np.array([np.nanmean(discharge[w_year == wy]) for wy in w_year_list])

        if self.std_by_w_year.size == 0:
            self.std_by_w_year = np.array([np.nanstd(discharge[w_year == wy], ddof=1) for wy in w_year_list])

        if self.cv_by_w_year.size == 0:
            self.cv_by_w_year = self.std_by_w_year / self.mean_by_w_year

        ma3 = np.nanmean(self.cv_by_w_year) * 100

        return ma3

    def get_ma4_11(self):

        discharge = self.discharge

        if self.mean_discharge.size == 0:
            self.mean_discharge = np.nanmean(discharge)

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        if self.median_log10_discharge.size == 0:
            self.median_log10_discharge = np.nanmedian(np.log10(discharge))

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_log10_percentiles.size == 0:
            self.discharge_log10_percentiles = get_percentile(np.log10(discharge),
                                                              self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        mean_prct = np.nanmean(self.discharge_log10_percentiles)
        std_prct = np.nanstd(self.discharge_log10_percentiles, ddof=1)

        ma4 = std_prct / mean_prct * 100
        ma5 = self.mean_discharge / self.median_discharge
        ma6 = self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 10] / \
              self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 90]
        ma7 = self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 20] / \
              self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 80]
        ma8 = self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 25] / \
              self.discharge_exceedance_percentiles[self.discharge_exceedance_percentiles_list == 75]
        ma9 = (self.discharge_log10_percentiles[self.discharge_percentiles_list == 90] -
               self.discharge_log10_percentiles[self.discharge_percentiles_list == 10]) / self.median_log10_discharge
        ma10 = (self.discharge_log10_percentiles[self.discharge_percentiles_list == 80] -
                self.discharge_log10_percentiles[self.discharge_percentiles_list == 20]) / self.median_log10_discharge
        ma11 = (self.discharge_log10_percentiles[self.discharge_percentiles_list == 75] -
                self.discharge_log10_percentiles[self.discharge_percentiles_list == 25]) / self.median_log10_discharge

        return ma4, ma5, ma6[0], ma7[0], ma8[0], ma9[0], ma10[0], ma11[0]

    def get_ma12_23(self):

        month = self.month
        discharge = self.discharge

        if self.mean_by_month.size == 0:
            self.mean_by_month = np.array([np.nanmean(discharge[(month == mm)]) for mm in month.unique()])

        mh12_23 = self.mean_by_month

        return mh12_23.tolist()

    def get_ma24_35(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_month_year.size == 0:
            self.mean_month_year = get_fun_by_month_year('mean', discharge, w_year, w_year_list, month)

        if self.std_month_year.size == 0:
            self.std_month_year = get_fun_by_month_year('std', discharge, w_year, w_year_list, month)

        cv_month_year = self.std_month_year / self.mean_month_year
        ma24_35 = np.nanmean(cv_month_year, axis=0) * 100

        return ma24_35.tolist()

    def get_ma36_40(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_month_year.size == 0:
            self.mean_month_year = get_fun_by_month_year('mean', discharge, w_year, w_year_list, month)

        mean_by_month = self.mean_month_year.reshape(-1, 1)

        median_mean_by_month = np.nanmedian(mean_by_month)
        mean_mean_by_month = np.nanmean(mean_by_month)

        p = np.array([10, 25, 75, 90]) / 100
        percentiles = get_percentile(mean_by_month, p)

        ma36 = (np.nanmax(mean_by_month) - np.nanmin(mean_by_month)) / median_mean_by_month
        ma37 = (percentiles[2] - percentiles[1]) / median_mean_by_month
        ma38 = (percentiles[3] - percentiles[0]) / median_mean_by_month
        ma39 = np.nanstd(mean_by_month, ddof=1) / mean_mean_by_month * 100
        ma40 = (mean_mean_by_month - median_mean_by_month) / median_mean_by_month

        return ma36, ma37, ma38, ma39, ma40

    def get_ma41_45(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge
        drainage_area = self.drainage_area

        if self.mean_by_w_year.size == 0:
            self.mean_by_w_year = np.array([np.nanmean(discharge[w_year == wy]) for wy in w_year_list])

        p = np.array([10, 25, 75, 90]) / 100
        percentiles = get_percentile(self.mean_by_w_year, p)

        mean_mean_by_w_year = np.nanmean(self.mean_by_w_year)
        median_mean_by_w_year = np.nanmedian(self.mean_by_w_year)

        ma41 = mean_mean_by_w_year / drainage_area
        ma42 = (np.nanmax(self.mean_by_w_year) - np.nanmin(self.mean_by_w_year)) / median_mean_by_w_year
        ma43 = (percentiles[2] - percentiles[1]) / median_mean_by_w_year
        ma44 = (percentiles[3] - percentiles[0]) / median_mean_by_w_year
        ma45 = (mean_mean_by_w_year - median_mean_by_w_year) / median_mean_by_w_year

        return ma41, ma42, ma43, ma44, ma45

    def get_mh1_12(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_month_year.size == 0:
            self.max_month_year = get_fun_by_month_year('max', discharge, w_year, w_year_list, month)

        mh1_12 = np.nanmean(self.max_month_year, axis=0)

        return mh1_12.tolist()

    def get_mh13(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_month_year.size == 0:
            self.max_month_year = get_fun_by_month_year('max', discharge, w_year, w_year_list, month)

        mh13 = np.nanstd(self.max_month_year.reshape(-1, 1), ddof=1) / np.nanmean(
            self.max_month_year.reshape(-1, 1)) * 100

        return mh13

    def get_mh14(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.median_by_w_year.size == 0:
            self.median_by_w_year = np.array([np.nanmedian(discharge[w_year == wy]) for wy in w_year_list])

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        mh14 = np.nanmedian(self.max_by_w_year / self.median_by_w_year)

        return mh14

    def get_mh15_17(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        if self.discharge_exceedance_percentiles_list.size == 0:
            self.discharge_exceedance_percentiles_list = np.array([100 - x for x in self.discharge_percentiles_list])

        if self.discharge_exceedance_percentiles.size == 0:
            self.discharge_exceedance_percentiles = self.discharge_percentiles

        discharge_exceedance_prct = get_percentile(discharge, 1 - 1 / 100)

        mh15 = discharge_exceedance_prct / self.median_discharge
        mh16 = self.discharge_exceedance_percentiles[
                   self.discharge_exceedance_percentiles_list == 10] / self.median_discharge
        mh17 = self.discharge_exceedance_percentiles[
                   self.discharge_exceedance_percentiles_list == 25] / self.median_discharge

        return mh15, mh16[0], mh17[0]

    def get_mh18(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        log10_max_by_w_year = np.log10(self.max_by_w_year)
        mh18 = np.nanstd(log10_max_by_w_year, ddof=1) / np.nanmean(log10_max_by_w_year) * 100

        return mh18

    def get_mh19(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        log10_max_by_w_year = np.log10(self.max_by_w_year)
        n_w_years = w_year_list.size

        np.sumLog10MaxDischargeByWYear = np.sum(log10_max_by_w_year)
        np.sumLog10MaxDischargeByWYearP2 = np.sum(log10_max_by_w_year ** 2)
        np.sumLog10MaxDischargeByWYearP3 = np.sum(log10_max_by_w_year ** 3)

        nom1 = n_w_years ** 2 * np.sumLog10MaxDischargeByWYearP3
        nom2 = 3 * n_w_years * np.sumLog10MaxDischargeByWYear * np.sumLog10MaxDischargeByWYearP2
        nom3 = 2 * np.sumLog10MaxDischargeByWYear ** 3

        denom = n_w_years * (n_w_years - 1) * (n_w_years - 2) * np.nanstd(log10_max_by_w_year, ddof=1) ** 3

        mh19 = (nom1 - nom2 + nom3) / denom

        return mh19

    def get_mh20(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge
        drainage_area = self.drainage_area

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        mh20 = np.nanmean(self.max_by_w_year) / drainage_area

        return mh20

    def get_mh21(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > self.median_discharge
        total_flow = np.sum(np.array([x - self.median_discharge for x in discharge[mask]]))
        n_events = np.sum(np.diff(np.append(0, mask)) > 0)

        mh21 = total_flow / n_events / self.median_discharge

        if np.isnan(mh21):
            mh21 = 0

        return mh21

    def get_mh22(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > (3 * self.median_discharge)
        total_flow = np.sum(np.array([x - (3 * self.median_discharge) for x in discharge[mask]]))
        n_events = np.sum(np.diff(np.append(0, mask)) > 0)

        mh22 = total_flow / n_events / self.median_discharge
        
        if np.isnan(mh22):
            mh22 = 0

        return mh22

    def get_mh23(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > (7 * self.median_discharge)
        total_flow = np.sum(np.array([x - (7 * self.median_discharge) for x in discharge[mask]]))
        n_events = np.sum(np.diff(np.append(0, mask)) > 0)

        mh23 = total_flow / n_events / self.median_discharge

        if np.isnan(mh23):
            mh23 = 0

        return mh23

    def get_mh24(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > self.median_discharge

        event_start = np.array(np.nonzero(np.diff(np.append(0, mask)) > 0))[0, :]
        event_end = np.array(np.nonzero(np.diff(np.append(mask, 0)) < 0))[0, :]

        if event_start.size != event_end.size:
            print("get_mh24(): events are not properly recognized")

        peaks_vec = [max(discharge[start:(end + 1)]) for start, end in zip(event_start, event_end)]

        mh24 = np.nanmean(peaks_vec) / self.median_discharge

        return mh24

    def get_mh25(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > (3 * self.median_discharge)

        event_start = np.array(np.nonzero(np.diff(np.append(0, mask)) > 0))[0, :]
        event_end = np.array(np.nonzero(np.diff(np.append(mask, 0)) < 0))[0, :]

        if event_start.size != event_end.size:
            print("get_mh25(): events are not properly recognized")

        peaks_vec = [max(discharge[start:(end + 1)]) for start, end in zip(event_start, event_end)]

        mh25 = np.nanmean(peaks_vec) / self.median_discharge

        return mh25

    def get_mh26(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        mask = discharge > (7 * self.median_discharge)

        event_start = np.array(np.nonzero(np.diff(np.append(0, mask)) > 0))[0, :]
        event_end = np.array(np.nonzero(np.diff(np.append(mask, 0)) < 0))[0, :]

        if event_start.size != event_end.size:
            print("get_mh26(): events are not properly recognized")

        peaks_vec = [max(discharge[start:(end + 1)]) for start, end in zip(event_start, event_end)]

        mh26 = np.nanmean(peaks_vec) / self.median_discharge

        if np.isnan(mh26):
            mh26 = 0

        return mh26

    def get_mh27(self):

        discharge = self.discharge

        if self.median_discharge.size == 0:
            self.median_discharge = np.nanmedian(discharge)

        if self.discharge_percentiles_list.size == 0:
            self.discharge_percentiles_list = np.linspace(95, 5, 19).reshape(-1, 1)

        if self.discharge_percentiles.size == 0:
            self.discharge_percentiles = get_percentile(discharge, self.discharge_percentiles_list / 100)

        mask = discharge > self.discharge_percentiles[self.discharge_percentiles_list == 75]

        event_start = np.array(np.nonzero(np.diff(np.append(0, mask)) > 0))[0, :]
        event_end = np.array(np.nonzero(np.diff(np.append(mask, 0)) < 0))[0, :]

        if event_start.size != event_end.size:
            print("get_mh27(): events are not properly recognized")

        peaks_vec = [max(discharge[start:(end + 1)]) for start, end in zip(event_start, event_end)]

        mh27 = np.nanmean(peaks_vec) / self.median_discharge

        return mh27

    def get_ml1_12(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_month_year.size == 0:
            self.min_month_year = get_fun_by_month_year('min', discharge, w_year, w_year_list, month)

        ml1_12 = np.nanmean(self.min_month_year, axis=0)

        return ml1_12.tolist()

    def get_ml13(self):

        month = self.month
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_month_year.size == 0:
            self.min_month_year = get_fun_by_month_year('min', discharge, w_year, w_year_list, month)

        ml13 = np.nanstd(self.min_month_year.reshape(-1, 1), ddof=1) / np.nanmean(
            self.min_month_year.reshape(-1, 1)) * 100

        return ml13

    def get_ml14_16(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_by_w_year.size == 0:
            self.mean_by_w_year = np.array([np.nanmean(discharge[w_year == wy]) for wy in w_year_list])

        if self.median_by_w_year.size == 0:
            self.median_by_w_year = np.array([np.nanmedian(discharge[w_year == wy]) for wy in w_year_list])

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        ml14 = np.nanmean(self.min_by_w_year / self.median_by_w_year)
        ml15 = np.nanmean(self.min_by_w_year / self.mean_by_w_year)
        ml16 = np.nanmedian(self.min_by_w_year / self.median_by_w_year)

        return ml14, ml15, ml16

    def get_ml17(self):

        if self.base_flow_index.size == 0:
            self.get_base_flow_index()

        ml17 = np.nanmean(self.base_flow_index)

        return ml17

    def get_ml18(self):

        if self.base_flow_index.size == 0:
            self.get_base_flow_index()

        ml18 = np.nanstd(self.base_flow_index, ddof=1) / np.nanmean(self.base_flow_index) * 100

        return ml18

    def get_ml19(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.mean_by_w_year.size == 0:
            self.mean_by_w_year = np.array([np.nanmean(discharge[w_year == wy]) for wy in w_year_list])

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        ml19 = np.nanmean(self.min_by_w_year / self.mean_by_w_year) * 100

        return ml19

    def get_ml20(self):

        discharge = self.discharge

        n5_days_blocks = int(np.floor(discharge.size / 5))

        if n5_days_blocks == 0:
            print('get_ml20(): not enough data')

        discharge = discharge[0:5 * n5_days_blocks].reshape(5, -1, order='F')
        min_blocks = np.nanmin(discharge, axis=0)

        mask_left = 0.9 * min_blocks[1:-1] < min_blocks[0:-2]
        mask_right = 0.9 * min_blocks[1:-1] < min_blocks[2:]
        mask = np.append(True, np.append(mask_left & mask_right, True))

        x = np.arange(0, n5_days_blocks)
        min_blocks[~mask] = np.interp(x[~mask], x[mask], min_blocks[mask])

        total_flow = np.sum(np.sum(discharge))
        total_block_flow = np.sum(min_blocks * 5)

        ml20 = total_block_flow / total_flow

        return ml20

    def get_ml21(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        ml21 = np.nanstd(self.min_by_w_year, ddof=1) / np.nanmean(self.min_by_w_year) * 100

        return ml21

    def get_ml22(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge
        drainage_area = self.drainage_area

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        ml22 = np.nanmean(self.min_by_w_year) / drainage_area

        return ml22

    def get_ra1_9(self):
        
        warnings.filterwarnings("ignore")
        eps = np.finfo(float).eps

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        d_discharge = np.diff(discharge, axis=0)
        positive_change_mask = d_discharge > 0
        negative_change_mask = d_discharge < 0

        positive_changes = d_discharge[positive_change_mask]
        negative_changes = np.abs(d_discharge[negative_change_mask])

        ra1 = np.nanmean(positive_changes)
        ra2 = np.nanstd(positive_changes, ddof=1) / ra1 * 100

        ra3 = np.nanmean(negative_changes)
        ra4 = np.nanstd(negative_changes, ddof=1) / ra3 * 100

        n_flow_records = len(discharge)
        ra5 = len(positive_changes) / n_flow_records

        log_discharge = np.log(discharge + eps)
        d_log_discharge = np.diff(log_discharge, axis=0)
        ra6 = np.nanmedian(d_log_discharge[positive_change_mask])
        ra7 = np.nanmedian(np.abs(d_log_discharge[negative_change_mask]))

        if ra7.size == 0 or np.isnan(ra7):
            ra7 = np.array(9999)

        pcm = -np.ones((len(positive_change_mask), 1))
        pcm[positive_change_mask] = 1
        ncm = np.zeros((len(negative_change_mask), 1))
        ncm[negative_change_mask] = 1

        change_mask = pd.DataFrame(pcm + ncm)
        change_mask.loc[change_mask[0] == -1] = None
        change_mask = change_mask.fillna(method='ffill').fillna(method='bfill').to_numpy()
        change_mask = np.append(change_mask[0], change_mask)

        d_changes = np.diff(change_mask)
        d_changes[d_changes != 0] = 1
        reversals = np.cumsum(np.append(0, d_changes))

        n_reversals_by_wy = np.array([np.max(reversals[w_year == wy]-reversals[w_year == wy][0]) for wy in w_year_list])

        ra8 = np.nanmean(n_reversals_by_wy)
        ra9 = np.nanstd(n_reversals_by_wy, ddof=1) / ra8 * 100
        if np.isnan(ra9):
            ra9 = 0

        return ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9

    def get_ta1_2(self):

        colwell_mat = self.get_colwell_mat().values

        x = np.sum(colwell_mat, axis=0)
        y = np.sum(colwell_mat, axis=1)
        z = np.sum(x)

        hx = -np.sum((x[x > 0] / z) * (np.log10(x[x > 0] / z)))
        hy = -np.sum((y[y > 0] / z) * (np.log10(y[y > 0] / z)))
        hxy = -np.sum((colwell_mat[colwell_mat > 0] / z) * (np.log10(colwell_mat[colwell_mat > 0] / z)))

        ta1 = 1 - hy / np.log10(11)  # Constancy
        ta2 = (1 - (hxy - hx) / np.log10(11)) * 100  # Predictability

        return ta1, ta2

    def get_ta3(self):

        month = self.month
        discharge = self.discharge
        wy_month = self.wy_month

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        periods1 = np.arange(wy_month, 13)
        periods2 = np.arange(1, wy_month if wy_month > 1 else 2)
        periods = np.concatenate((periods1, periods2))
        periods = np.array(unique_order(periods)).reshape(6, -1)

        n_events_day_by_period = [np.sum(discharge[(month == period[0]) | (month == period[1])] > thresh)
                                  for period in periods]

        ta3 = np.nanmax(n_events_day_by_period) / np.sum(n_events_day_by_period)

        return ta3

    def get_th1_2(self):

        day = self.day
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        max_days = np.array(
            [min(day[w_year == wy].reshape(-1, 1)[discharge[w_year == wy] == self.max_by_w_year[i]]) for i, wy in
             enumerate(w_year_list)])

        rad2day = 2 * np.pi / 365.25
        x = np.cos(max_days * rad2day)
        y = np.sin(max_days * rad2day)

        theta = np.arctan2(np.nanmean(y), np.nanmean(x))

        if theta < 0:
            theta = 2 * np.pi + theta

        th1 = np.round(theta / rad2day, 0)
        th2 = np.sqrt(2 * (1 - np.linalg.norm([np.nanmean(x), np.nanmean(y)]))) / rad2day

        return th1, th2

    def get_th3(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 40]

        th3 = np.nanmax(np.array([np.sum(discharge[w_year == wy] < thresh) / 365 for wy in w_year_list]))

        return th3

    def get_tl1_2(self):

        day = self.day
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.min_by_w_year.size == 0:
            self.min_by_w_year = np.array([np.nanmin(discharge[w_year == wy]) for wy in w_year_list])

        min_days = np.array(
            [min(day[w_year == wy].reshape(-1, 1)[discharge[w_year == wy] == self.min_by_w_year[i]]) for i, wy in
             enumerate(w_year_list)])

        rad2day = 2 * np.pi / 365.25
        x = np.cos(min_days * rad2day)
        y = np.sin(min_days * rad2day)

        theta = np.arctan2(np.nanmean(y), np.nanmean(x))

        if theta < 0:
            theta = 2 * np.pi + theta

        tl1 = np.round(theta / rad2day, 0)
        tl2 = np.sqrt(2 * (1 - np.linalg.norm([np.nanmean(x), np.nanmean(y)]))) / rad2day

        return tl1, tl2

    def get_tl3(self):

        month = self.month
        discharge = self.discharge
        wy_month = self.wy_month

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 80]

        periods1 = np.arange(wy_month, 13)
        periods2 = np.arange(1, wy_month if wy_month > 1 else 2)
        periods = np.concatenate((periods1, periods2))
        periods = np.array(unique_order(periods)).reshape(6, -1)

        n_events_day_by_period = [np.sum(discharge[(month == period[0]) | (month == period[1])] <= thresh)
                                  for period in periods]

        tl3 = np.nanmax(n_events_day_by_period) / np.sum(n_events_day_by_period)

        return tl3

    def get_tl4(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        self.get_peak_thresh(prct=np.array([0.4, 0.8]))

        thresh = self.peak_thresh[self.peak_thresh_prct == 80]

        tl4 = np.nanmax(
            np.array([np.sum(discharge[w_year == wy] > thresh) / (365 + int((wy % 4) == 0)) for wy in w_year_list]))

        return tl4

    def get_base_flow_index(self):

        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.rolling7day_mean.size == 0:
            self.rolling7day_mean = np.array(np.convolve(discharge[:, 0], np.ones((7,)) / 7, 'same'))
            self.rolling7day_mean[:3] = np.nan
            self.rolling7day_mean[-3:] = np.nan

        rolling7day_mean_by_w_year = np.array([self.rolling7day_mean[w_year == wy] for wy in w_year_list],
                                              dtype=object)

        if self.mean_by_w_year.size == 0:
            self.mean_by_w_year = np.array([np.nanmean(discharge[w_year == wy]) for wy in w_year_list])

        self.base_flow_index = np.array(
            [np.nanmin(a) / b for a, b in zip(rolling7day_mean_by_w_year, self.mean_by_w_year)])

        return self.base_flow_index

    def get_peak_thresh(self, prct):

        day = self.day
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        if self.max_by_w_year.size == 0:
            self.max_by_w_year = np.array([np.nanmax(discharge[w_year == wy]) for wy in w_year_list])

        if self.peak_thresh.size == 0:

            mean_daily_flow = np.zeros((len(w_year_list), 1))
            for i, wy in enumerate(w_year_list):
                masked_data = discharge[w_year == wy][:, 0]
                masked_days = day[w_year == wy]
                day_of_max_flow = masked_days[masked_data == self.max_by_w_year[i]][0]
                mean_daily_flow[i] = np.nanmean(masked_data[masked_days == day_of_max_flow])

            x = np.log10(self.max_by_w_year)
            y = np.log10(mean_daily_flow)[:, 0]

            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)

            thresh = get_percentile(x, prct)

            self.peak_thresh_prct = prct * 100
            self.peak_thresh = 10 ** p(thresh)

        peak_thresh = self.peak_thresh

        return peak_thresh

    def get_colwell_mat(self):

        daym = self.daym
        month = self.month
        discharge = self.discharge

        if self.mean_discharge.size == 0:
            self.mean_discharge = np.nanmean(discharge)

        eps = np.finfo(float).eps

        aux = [0.10, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25]

        if self.mean_discharge < 1:
            colwell_bounds = np.array(
                [(-0.9 * x + 1) * np.log10(self.mean_discharge + eps) for x in aux]).reshape(-1, 1)
        else:
            colwell_bounds = np.array([x * np.log10(self.mean_discharge + eps) for x in aux]).reshape(-1, 1)

        self.colwell_bounds = colwell_bounds

        cat = ["{:02d}".format(int(x)) for x in range(0, colwell_bounds.size + 1)]
        
        # quantize discharge values into pre-defined intervals
        q_discharge = [bisect.bisect_left(colwell_bounds, x) for x in np.log10(discharge)]
        q_discharge = ["{:02d}".format(int(x)) for x in q_discharge]
        q_discharge = pd.Categorical(q_discharge, categories=cat)

        dates = ["{:02d}_{:02d}".format(m, d) for m, d in zip(month, daym)]
        dates = pd.Categorical(dates, categories=np.unique(np.array(dates)))

        colwell_mat = pd.crosstab(q_discharge, dates)

        colwell_mat = colwell_mat.drop('02_29', axis=1)

        self.colwell_mat = colwell_mat

        return colwell_mat

    def get_ar1(self):

        month = self.month
        discharge = self.discharge

        if self.mean_by_month.size == 0:
            self.mean_by_month = np.array([np.nanmean(discharge[(month == mm)]) for mm in month.unique()])

        # Deseasonalize time-series using long-term monthly means
        ds_discharge = np.array([float(x - self.mean_by_month[m - 1]) for x, m in zip(discharge, month)])
        ds_discharge = (ds_discharge - ds_discharge.mean()) / ds_discharge.std(ddof=1)

        # from http://mpastell.com/pweave/_downloads/AR_yw.html
        ar1_cor = fit_ar(ds_discharge)

        return ar1_cor[0]

    def get_seasonality(self):

        day = self.day
        w_year = self.w_year
        w_year_list = self.w_year_list
        discharge = self.discharge

        # get days in water year
        w_day = [np.arange(1, len(day[w_year == wy]) + 1) for wy in w_year_list]
        w_ndays = [len(day[w_year == wy]) for wy in w_year_list]

        if w_ndays[0] < 366:
            aux = w_ndays[1:4]
            if 366 in aux:
                w_day[0] = np.arange(365 - w_ndays[0] + 1, 366)
            else:
                w_day[0] = np.arange(366 - w_ndays[0] + 1, 367)

        w_day = np.array(list(flatten(w_day)))

        # amplitude and phase computation
        decimal_year = np.array(w_year + w_day / 365.25)
        n = discharge.size
        y = ((discharge - discharge.mean()) / discharge.std(ddof=1)).reshape(n, 1)
        x = np.concatenate((np.ones((n, 1)), np.sin(2 * np.pi * decimal_year.reshape(n, 1)),
                            np.cos(2 * np.pi * decimal_year.reshape(n, 1))), axis=1)
        beta = np.linalg.solve(x.T.dot(x), x.T.dot(y))

        amplitude = np.sqrt(beta[1] ** 2 + beta[2] ** 2)
        phase = 365.25 * (np.arctan2(beta[1], beta[2])) / (2 * np.pi)
        phase = (phase + 365.25) % 365.25

        return amplitude, np.round(phase, 0)


def autocorr(x, lag=1):
    c = np.correlate(x, x, 'full')
    mid = len(c) // 2
    acov = c[mid:mid + lag]
    acor = acov / acov[0]

    return acor


def fit_ar(x, p=1):
    ac = autocorr(x, p + 1)
    r_mat = linalg.toeplitz(ac[:p])
    r = ac[1:p + 1]
    phi = linalg.inv(r_mat).dot(r)

    return phi


def get_percentile(values, prct):
    # create numpy array for discharge values
    values = np.array(values)
    # remove nan
    values = values[~np.isnan(values)]
    # sort values
    values_sorted = np.sort(values)
    # probabilities
    n = len(values)
    p = np.linspace(1, n, n) / (1 + n)
    # interpolate
    q = np.interp(prct, p, values_sorted)
    return q


def get_fun_by_month_year(fun, discharge, year, year_list, month):
    warnings.filterwarnings("ignore")

    nmonths = len(month.unique())
    nyears = len(year_list)
    result = np.zeros((nyears, nmonths,))
    for i, yy in enumerate(year_list):
        for j, mm in enumerate(month.unique()):
            try:
                if fun == 'mean':
                    result[i, j] = np.nanmean(discharge[(year == yy) & (month == mm)])
                elif fun == 'median':
                    result[i, j] = np.nanmedian(discharge[(year == yy) & (month == mm)])
                elif fun == 'min':
                    result[i, j] = np.nanmin(discharge[(year == yy) & (month == mm)])
                elif fun == 'max':
                    result[i, j] = np.nanmax(discharge[(year == yy) & (month == mm)])
                elif fun == 'std':
                    result[i, j] = np.nanstd(discharge[(year == yy) & (month == mm)], ddof=1)
                else:
                    print('Function does not exist')
            except ValueError:
                result[i, j] = np.nan

    return result


def unique_order(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def flatten(lt):
    for el in lt:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


# code from package lmoments3 by Sam Gillespie, 2014
def lmom_ratios(data, nmom=5):
    """
    Estimate `nmom` number of L-moments from a sample `data`.

    :param data: Sequence of (sample) data
    :type data: list or array-like sequence
    :param nmom: number of L-moments to estimate
    :type nmom: int
    :return: L-moment ratios like this: l1, l2, t3, t4, t5, .. . As in: items 3 and higher are L-moment ratios.
    :rtype: list
    """

    if nmom <= 5:
        return _samlmusmall(data, nmom)
    else:
        return _samlmularge(data, nmom)


def _samlmularge(x, nmom=5):
    try:
        x = np.asarray(x, dtype=np.float64)
        n = len(x)
        x.sort()
    except ValueError:
        raise ValueError("Input data to estimate L-moments must be numeric.")

    if nmom <= 0:
        raise ValueError("Invalid number of sample L-moments")

    if n < nmom:
        raise ValueError("Insufficient length of data for specified nmoments")

    # Calculate first order
    lm = [np.sum(x) / sp.comb(n, 1, exact=True)]

    if nmom == 1:
        return lm[0]

    # Setup comb table, where comb[i][x] refers to comb(x,i)
    comb = []
    for i in range(1, nmom):
        comb.append([])
        for j in range(n):
            comb[-1].append(sp.comb(j, i, exact=True))

    for mom in range(2, nmom + 1):
        coefl = 1.0 / mom * 1.0 / sp.comb(n, mom, exact=True)
        xtrans = []
        for i in range(0, n):
            coeftemp = []
            for j in range(0, mom):
                coeftemp.append(1)

            for j in range(0, mom - 1):
                coeftemp[j] = coeftemp[j] * comb[mom - j - 2][i]

            for j in range(1, mom):
                coeftemp[j] = coeftemp[j] * comb[j - 1][n - i - 1]

            for j in range(0, mom):
                coeftemp[j] = coeftemp[j] * sp.comb(mom - 1, j, exact=True)

            for j in range(0, int(0.5 * mom)):
                coeftemp[j * 2 + 1] = -coeftemp[j * 2 + 1]
            coeftemp = sum(coeftemp)
            xtrans.append(x[i] * coeftemp)

        if mom > 2:
            lm.append(coefl * sum(xtrans) / lm[1])
        else:
            lm.append(coefl * sum(xtrans))
    return lm


def _samlmusmall(x, nmom=5):
    try:
        x = np.asarray(x, dtype=np.float64)
        n = len(x)
        x = sorted(x)
    except ValueError:
        raise ValueError("Input data to estimate L-moments must be numeric.")

    if nmom <= 0 or nmom > 5:
        raise ValueError("Invalid number of sample L-moments")

    if n < nmom:
        raise ValueError("Insufficient length of data for specified nmoments")

    # First L-moment

    l1 = np.sum(x) / sp.comb(n, 1, exact=True)

    if nmom == 1:
        return l1

    # Second L-moment

    comb1 = range(n)
    coefl2 = 0.5 / sp.comb(n, 2, exact=True)
    sum_xtrans = sum([(comb1[i] - comb1[n - i - 1]) * x[i] for i in range(n)])
    l2 = coefl2 * sum_xtrans

    if nmom == 2:
        return [l1, l2]

    # Third L-moment

    comb3 = [sp.comb(i, 2, exact=True) for i in range(n)]
    coefl3 = 1.0 / 3.0 / sp.comb(n, 3, exact=True)
    sum_xtrans = sum([(comb3[i] - 2 * comb1[i] * comb1[n - i - 1] + comb3[n - i - 1]) * x[i] for i in range(n)])
    l3 = coefl3 * sum_xtrans / l2

    if nmom == 3:
        return [l1, l2, l3]

    # Fourth L-moment

    comb5 = [sp.comb(i, 3, exact=True) for i in range(n)]
    coefl4 = 0.25 / sp.comb(n, 4, exact=True)
    sum_xtrans = sum(
        [(comb5[i] - 3 * comb3[i] * comb1[n - i - 1] + 3 * comb1[i] * comb3[n - i - 1] - comb5[n - i - 1]) * x[i]
         for i in range(n)])
    l4 = coefl4 * sum_xtrans / l2

    if nmom == 4:
        return [l1, l2, l3, l4]

    # Fifth L-moment

    comb7 = [sp.comb(i, 4, exact=True) for i in range(n)]
    coefl5 = 0.2 / sp.comb(n, 5, exact=True)
    sum_xtrans = sum(
        [(comb7[i] - 4 * comb5[i] * comb1[n - i - 1] + 6 * comb3[i] * comb3[n - i - 1] -
          4 * comb1[i] * comb5[n - i - 1] + comb7[n - i - 1]) * x[i]
         for i in range(n)])
    l5 = coefl5 * sum_xtrans / l2

    return [l1, l2, l3, l4, l5]
