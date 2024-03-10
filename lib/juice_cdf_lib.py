import glob
import spacepy.pycdf
import numpy as np

class struct:
    pass

#---------------------------------------------------------------------
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

#---------------------------------------------------------------------
def juice_gethk_hf(data):

    hk = struct()
    hk.epoch = data['Epoch'][...]
    
    hk.heater_ena = data['LWT03314']
    hk.calsig_ena = data['LWT0332C']

    hk.deploy_pri_x=data['LWT0332E']
    hk.deploy_red_x=data['LWT0332F']
    hk.deploy_pri_y=data['LWT03330']
    hk.deploy_red_y=data['LWT03331']
    hk.deploy_pri_z=data['LWT03332']
    hk.deploy_red_z=data['LWT03333']
    hk.deploy_lock_stat=data['LWT03334']

    hk.temp_rwi_u = data['LWT03337_CALIBRATED'][...]
    hk.temp_rwi_w = data['LWT03339_CALIBRATED'][...]
    hk.temp_hf_fpga = data['LWT0333B_CALIBRATED'][...]

    return hk

#---------------------------------------------------------------------
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

#---------------------------------------------------------------------
def juice_gethk_lvps(data):

    hk = struct()
    hk.epoch = data['Epoch'][...]
    hk.vol_hf_33 = data['LWT03358_CALIBRATED'][...]
    hk.vol_hf_85 = data['LWT03359_CALIBRATED'][...]
    hk.cur_hf_33 = data['LWT03362_CALIBRATED'][...]
    hk.cur_hf_85 = data['LWT03363_CALIBRATED'][...]
    hk.hf_on_off = data['LWT03372'][...]

    return hk

#---------------------------------------------------------------------
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

#---------------------------------------------------------------------
def juice_getspec_hf_sid02(data, AFSW_ver=1):

    if AFSW_ver == 1:
        spec = juice_getspec_hf_sid02_ver01(data)

    return spec

#---------------------------------------------------------------------
def juice_getspec_hf_sid02_ver01(data):

    spec = struct()

    dt = 1.0/148000.0   # 1/bandwidth [sec]    (fixed for ver.1 SW SID2)
    n_step = 512        # number of sweep step (fixed for ver.1 SW SID2)
    n_samp1 = 32
    n_samp2 = 50
    n_samp3 = 64
    len1 = n_step * n_samp1
    len2 = n_step * n_samp2
    len3 = n_step * n_samp3

    sweep_start_1d = [x for row in data.sweep_start for x in row]
    Eu_i_1d = [x for row in data.Eu_i for x in row]
    Eu_q_1d = [x for row in data.Eu_q for x in row]
    Ev_i_1d = [x for row in data.Ev_i for x in row]
    Ev_q_1d = [x for row in data.Ev_q for x in row]
    Ew_i_1d = [x for row in data.Ew_i for x in row]
    Ew_q_1d = [x for row in data.Ew_q for x in row]
    frequency_1d = [x for row in data.frequency for x in row]

    index_1d = np.where(sweep_start_1d)
    index = np.where(data.sweep_start == 1)
    n = len(index_1d[0])
    print("number of sweeps : ", n)

    epoch = []
    Eu_power = []
    Ev_power = []
    Ew_power = []
    frequency = []

    for i in range(n-1):
        index_len = index_1d[0][i+1]-index_1d[0][i]
        if index_len == len1:
            n_samp = n_samp1
        elif index_len == len2:
            n_samp = n_samp2
        elif index_len == len3:
            n_samp = n_samp3
        else:
            print("invalid length : ", index_len)
            print("skipped")
            continue

        print("data length : ", n_samp)
        epoch.append(data.epoch[index[0][i]])

        Eu_i = np.array(Eu_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Eu_q = np.array(Eu_q_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ev_i = np.array(Ev_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ev_q = np.array(Ev_q_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ew_i = np.array(Ew_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ew_q = np.array(Ew_q_1d[index_1d[0][i]:index_1d[0][i+1]])

        Eu_i_array = Eu_i.reshape(n_step, n_samp)
        Eu_q_array = Eu_q.reshape(n_step, n_samp)
        Ev_i_array = Ev_i.reshape(n_step, n_samp)
        Ev_q_array = Ev_q.reshape(n_step, n_samp)
        Ew_i_array = Ew_i.reshape(n_step, n_samp)
        Ew_q_array = Ew_q.reshape(n_step, n_samp)

        # low resolution power spectra
        Eu_power.append(np.mean(Eu_i_array**2 + Eu_q_array**2, axis=1))
        Ev_power.append(np.mean(Ev_i_array**2 + Ev_q_array**2, axis=1))
        Ew_power.append(np.mean(Ew_i_array**2 + Ew_q_array**2, axis=1))

        frequency = np.array(frequency_1d[index_1d[0][i]:index_1d[0][i]+len1])
        frequency = frequency.reshape(n_step, n_samp1)
        frequency = frequency[:, 0]

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

#---------------------------------------------------------------------
def juice_getspec_hf_sid02_highres_ver01(data):

    spec = struct()

    dt = 1.0/148000.0   # 1/bandwidth [sec]    (fixed for ver.1 SW SID2)
    n_step = 512        # number of sweep step (fixed for ver.1 SW SID2)
    n_samp1 = 32
    n_samp2 = 50
    n_samp3 = 512
    len1 = n_step * n_samp1
    len2 = n_step * n_samp2
    len3 = n_step * n_samp3

    sweep_start_1d = [x for row in data.sweep_start for x in row]
    Eu_i_1d = [x for row in data.Eu_i for x in row]
    Eu_q_1d = [x for row in data.Eu_q for x in row]
    Ev_i_1d = [x for row in data.Ev_i for x in row]
    Ev_q_1d = [x for row in data.Ev_q for x in row]
    Ew_i_1d = [x for row in data.Ew_i for x in row]
    Ew_q_1d = [x for row in data.Ew_q for x in row]
    frequency_1d = [x for row in data.frequency for x in row]

    index_1d = np.where(sweep_start_1d)
    index = np.where(data.sweep_start == 1)
    n = len(index_1d[0])
    print("number of sweeps : ",n)

    epoch = []

    Eu_power_high = []
    Ev_power_high = []
    Ew_power_high = []
    frequency_high = []

    for i in range(n-1):
        index_len = index_1d[0][i+1]-index_1d[0][i]
        if index_len == len1:
            n_samp = n_samp1
        elif index_len == len2:
            n_samp = n_samp2
        elif index_len == len3:
            n_samp = n_samp3
        else:
            print("invalid length : ", index_len)
            print("skipped")
            continue

        print("data length : ", n_samp)
        epoch.append(data.epoch[index[0][i]])

        Eu_i = np.array(Eu_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Eu_q = np.array(Eu_q_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ev_i = np.array(Ev_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ev_q = np.array(Ev_q_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ew_i = np.array(Ew_i_1d[index_1d[0][i]:index_1d[0][i+1]])
        Ew_q = np.array(Ew_q_1d[index_1d[0][i]:index_1d[0][i+1]])

        Eu_i_array = Eu_i.reshape(n_step, n_samp)
        Eu_q_array = Eu_q.reshape(n_step, n_samp)
        Ev_i_array = Ev_i.reshape(n_step, n_samp)
        Ev_q_array = Ev_q.reshape(n_step, n_samp)
        Ew_i_array = Ew_i.reshape(n_step, n_samp)
        Ew_q_array = Ew_q.reshape(n_step, n_samp)

        # high resolution power spectra
        for ii in range(n_step):
            y = Eu_i_array[ii][:] + Eu_q_array[ii][:]*1j
            s = np.fft.fft(y)
            power = np.power(np.abs(s)/(n_samp/2), 2.0)
            Eu_power_high.append(power)

            y = Ev_i_array[ii][:] + Ev_q_array[ii][:]*1j
            s = np.fft.fft(y)
            power = np.power(np.abs(s)/(n_samp/2), 2.0)
            Ev_power_high.append(power)

            y = Ew_i_array[ii][:] + Ew_q_array[ii][:]*1j
            s = np.fft.fft(y)
            power = np.power(np.abs(s)/(n_samp/2), 2.0)
            Ew_power_high.append(power)

            if (ii == 0):
                freq = np.fft.fftfreq(n_samp, d=dt) + frequency[ii]
                frequency_high.append(freq)


    Eu_power_high = np.array(Eu_power_high)
    Ev_power_high = np.array(Ev_power_high)
    Ew_power_high = np.array(Ew_power_high)

    n_set = int(Eu_power_high.size/n_step)

    spec.frequency_high = np.array(frequency_high)
    spec.Eu_power_high = Eu_power_high.reshape(n_set, n_step).transpose()
    spec.Ev_power_high = Ev_power_high.reshape(n_set, n_step).transpose()
    spec.Ew_power_high = Ew_power_high.reshape(n_set, n_step).transpose()

    return spec

#---------------------------------------------------------------------
def juice_getdata_hf_sid03(cdf):

    spec = struct()

    spec.epoch = cdf['Epoch'][...]
    spec.scet = cdf['SCET'][...]
    spec.frequency = cdf['frequency'][...]

    spec.EuEu = cdf['EuEu'][...]
    spec.EvEv = cdf['EvEv'][...]
    spec.EwEw = cdf['EwEw'][...]

    fill_value = -1e30 

    idx = np.where(spec.EuEu == fill_value)
    spec.EuEu[idx] = np.nan
    idx = np.where(spec.EvEv == fill_value)
    spec.EvEv[idx] = np.nan
    idx = np.where(spec.EwEw == fill_value)
    spec.EwEw[idx] = np.nan

    idx = np.where(spec.frequency == fill_value)
    idx = sorted(idx,reverse=True)

    for index in idx:
        del spec.frequency[index]
        for row in spec.EuEu:
            del row[index]
        for row in spec.EvEv:
            del row[index]
        for row in spec.EwEw:
            del row[index]

    return spec
