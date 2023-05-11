import glob
import spacepy.pycdf
import numpy as np


class struct:
    pass


def juice_read_cdfs(date_str, label, ver_str="01", base_dir="/db/JUICE/juice/datasets/"):

    yr_str = date_str[0:4]
    mn_str = date_str[4:6]
    dy_str = date_str[6:8]
    search_path = base_dir+yr_str+'/'+mn_str+'/'+dy_str + \
        '/JUICE_LU_RPWI-PPTD-'+label+'_'+date_str+'T??????_V'+ver_str+'.cdf'

    fname = glob.glob(search_path)
    if len(fname) > 0:
        err = 0
        ret = spacepy.pycdf.concatCDF([
            spacepy.pycdf.CDF(f) for f in glob.glob(search_path)])
    else:
        err = 1
        ret = 0

    return ret, err


def juice_gethk_hf(data):

    hk = struct()
    hk.epoch = data['Epoch'][...]
    hk.temp_rwi_u = data['LWT03337_CALIBRATED'][...]
    hk.temp_rwi_w = data['LWT03339_CALIBRATED'][...]
    hk.temp_hf_fpga = data['LWT0333B_CALIBRATED'][...]

    return hk


def juice_gethk_dpu(data):

    hk = struct()
    hk.epoch = data['Epoch'][...]
    hk.dpu_temp = data['LWT03437_CALIBRATED'][...]
    hk.lvps_temp = data['LWT03438_CALIBRATED'][...]
    hk.lp_temp = data['LWT03439_CALIBRATED'][...]
    hk.lf_temp = data['LWT0343A_CALIBRATED'][...]
    hk.hf_temp = data['LWT0343B_CALIBRATED'][...]
    hk.scm_temp = data['LWT0343C_CALIBRATED'][...]

    return hk


def juice_gethk_lvps(data):

    hk = struct()
    hk.epoch = data['Epoch'][...]
    hk.vol_hf_33 = data['LWT03358_CALIBRATED'][...]
    hk.vol_hf_85 = data['LWT03359_CALIBRATED'][...]
    hk.cur_hf_33 = data['LWT03362_CALIBRATED'][...]
    hk.cur_hf_85 = data['LWT03363_CALIBRATED'][...]
    hk.hf_on_off = data['LWT03372'][...]

    return hk


def juice_getdata_hf_sid02(cdf):

    data = struct()

    data.epoch = cdf['Epoch'][...]
    data.scet = cdf['SCET'][...]
    data.Eu_i = cdf['Eu_i'][...]
    data.Eu_q = cdf['Eu_q'][...]
    data.Ev_i = cdf['Ev_i'][...]
    data.Ev_q = cdf['Ev_q'][...]
    data.Ew_i = cdf['Ew_i'][...]
    data.Ew_q = cdf['Ew_q'][...]
    data.overflow = cdf['overflow'][...]
    data.pps_count = cdf['pps_count'][...]
    data.reduction = cdf['reduction'][...]
    data.sweep_start = cdf['sweep_start'][...]
    data.freq_width = cdf['freq_width'][...]
    data.frequency = cdf['frequency'][...]

    return data


def juice_getspec_hf_sid02(data, AFSW_ver=1):

    if AFSW_ver == 1:

        spec = juice_getspec_hf_sid02_ver01_(data)

     # check sweep start position
#     index = np.where(data.sweep_start == 1)
#      n_blk = (index[0][1]-index[0][0])
#       if n_blk == 32*512:
#            correction = False
#        else:
#            correction = True
#
#        if correction == False:
#            spec = juice_getspec_hf_sid02_ver01(data)
#        else:
#            spec = juice_getspec_hf_sid02_ver01_correction(data)

    return spec


def juice_getspec_hf_sid02_ver01(data):

    spec = struct()

    n_step = 512
    n_samp = 32
    n_set = int(len(data.Eu_i)/n_samp/n_step)

    frequency_array = data.frequency.reshape(n_set, n_step, n_samp)
    epoch_array = data.epoch.reshape(n_set, n_step, n_samp)
    spec.frequency = frequency_array[0, :, 0]
    spec.epoch = epoch_array[:, 0, 0]

    Eu_i_array = data.Eu_i.reshape(n_set, n_step, n_samp)
    Eu_q_array = data.Eu_q.reshape(n_set, n_step, n_samp)
    Ev_i_array = data.Ev_i.reshape(n_set, n_step, n_samp)
    Ev_q_array = data.Ev_q.reshape(n_set, n_step, n_samp)
    Ew_i_array = data.Ew_i.reshape(n_set, n_step, n_samp)
    Ew_q_array = data.Ew_q.reshape(n_set, n_step, n_samp)

    Eu_power = Eu_i_array**2 + Eu_q_array**2
    Ev_power = Ev_i_array**2 + Ev_q_array**2
    Ew_power = Ew_i_array**2 + Ew_q_array**2

    spec.Eu_power = (np.mean(Eu_power, axis=2)).transpose()
    spec.Ev_power = (np.mean(Ev_power, axis=2)).transpose()
    spec.Ew_power = (np.mean(Ew_power, axis=2)).transpose()

    return spec


def juice_getspec_hf_sid02_ver01_correction(data):

    spec = struct()

    n_step = 512
    n_samp = 50
    n_samp_ = 32

    index = np.where(data.sweep_start == 1)
    n = len(index[0])

    Eu_i = []
    Eu_q = []
    Ev_i = []
    Ev_q = []
    Ew_i = []
    Ew_q = []
    frequency = []
    epoch = []

    for i in range(n-1):
        if (index[0][i+1]-index[0][i]) == n_step*n_samp:
            frequency.append(data.frequency[index[0][i]:index[0][i+1]-1])
            epoch.append(data.epoch[index[0][i]])

            Eu_i.append(data.Eu_i[index[0][i]:index[0][i+1]])
            Eu_q.append(data.Eu_q[index[0][i]:index[0][i+1]])
            Ev_i.append(data.Ev_i[index[0][i]:index[0][i+1]])
            Ev_q.append(data.Ev_q[index[0][i]:index[0][i+1]])
            Ew_i.append(data.Ew_i[index[0][i]:index[0][i+1]])
            Ew_q.append(data.Ew_q[index[0][i]:index[0][i+1]])

    Eu_i = np.float_(Eu_i)
    Eu_q = np.float_(Eu_q)
    Ev_i = np.float_(Ev_i)
    Ev_q = np.float_(Ev_q)
    Ew_i = np.float_(Ew_i)
    Ew_q = np.float_(Ew_q)

    n_set = int(Eu_i.size/n_samp/n_step)
    n_set_ = int(len(data.Eu_i)/n_samp_/n_step)

    frequency_array = data.frequency.reshape(n_set_, n_step, n_samp_)
    spec.frequency = frequency_array[0, :, 0]
    spec.epoch = epoch

    Eu_i_array = Eu_i.reshape(n_set, n_step, n_samp)
    Eu_q_array = Eu_q.reshape(n_set, n_step, n_samp)
    Ev_i_array = Ev_i.reshape(n_set, n_step, n_samp)
    Ev_q_array = Ev_q.reshape(n_set, n_step, n_samp)
    Ew_i_array = Ew_i.reshape(n_set, n_step, n_samp)
    Ew_q_array = Ew_q.reshape(n_set, n_step, n_samp)

    Eu_power = Eu_i_array**2 + Eu_q_array**2
    Ev_power = Ev_i_array**2 + Ev_q_array**2
    Ew_power = Ew_i_array**2 + Ew_q_array**2

    spec.Eu_power = (np.mean(Eu_power, axis=2)).transpose()
    spec.Ev_power = (np.mean(Ev_power, axis=2)).transpose()
    spec.Ew_power = (np.mean(Ew_power, axis=2)).transpose()

    return spec

def juice_getspec_hf_sid02_ver01_(data):

    spec = struct()

    n_step = 512
    n_samp1 = 32
    n_samp2 = 50
    len1 = n_step * n_samp1
    len2 = n_step * n_samp2

    index = np.where(data.sweep_start == 1)
    n = len(index[0])

    Eu_power = []
    Ev_power = []
    Ew_power = []
    frequency = []
    epoch = []

    for i in range(n-1):
        index_len = index[0][i+1]-index[0][i]
        if index_len == len1:
            n_samp = n_samp1
        elif index_len == len2:
            n_samp = n_samp2
        else:
            continue    

        epoch.append(data.epoch[index[0][i]])

        Eu_i = data.Eu_i[index[0][i]:index[0][i+1]]
        Eu_q = data.Eu_q[index[0][i]:index[0][i+1]]
        Ev_i = data.Ev_i[index[0][i]:index[0][i+1]]
        Ev_q = data.Ev_q[index[0][i]:index[0][i+1]]
        Ew_i = data.Ew_i[index[0][i]:index[0][i+1]]
        Ew_q = data.Ew_q[index[0][i]:index[0][i+1]]

        Eu_i_array = Eu_i.reshape(n_step, n_samp)
        Eu_q_array = Eu_q.reshape(n_step, n_samp)
        Ev_i_array = Ev_i.reshape(n_step, n_samp)
        Ev_q_array = Ev_q.reshape(n_step, n_samp)
        Ew_i_array = Ew_i.reshape(n_step, n_samp)
        Ew_q_array = Ew_q.reshape(n_step, n_samp)

        Eu_power.append(np.mean(Eu_i_array**2 + Eu_q_array**2, axis=1))
        Ev_power.append(np.mean(Ev_i_array**2 + Ev_q_array**2, axis=1))
        Ew_power.append(np.mean(Ew_i_array**2 + Ew_q_array**2, axis=1))

        frequency = data.frequency[index[0][i]:index[0][i]+len1]
        frequency = frequency.reshape(n_step, n_samp1)
        frequency = frequency[:,0]

    Eu_power = np.float_(Eu_power)
    Ev_power = np.float_(Ev_power)
    Ew_power = np.float_(Ew_power)

    n_set = int(Eu_power.size/n_step)

    spec.frequency = frequency
    spec.epoch = epoch
    spec.Eu_power = Eu_power.reshape(n_set, n_step).transpose()
    spec.Ev_power = Ev_power.reshape(n_set, n_step).transpose()
    spec.Ew_power = Ew_power.reshape(n_set, n_step).transpose()

    return spec
