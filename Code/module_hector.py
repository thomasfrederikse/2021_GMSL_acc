# -------------------------------------------------------------------------
# Functions to run Hector from the main program
# Check whether 'hector_path' points to the path where the Hector binaries
# are installed. These can be obtained from http://segal.ubi.pt/hector/
#
# The research was carried out in part at the Jet Propulsion Laboratory,
# California Institute of Technology, under a contract with the National
# Aeronautics and Space Administration
#
# (c) All rights reserved
# -------------------------------------------------------------------------
import numpy as np
import os

def time2mjd(time):
    year = np.floor(time)
    day = 365.25*(time-year)
    mjd = (365.25*year+day) - 678942
    return mjd


def time2mjd_month(time):
    year = np.floor(time)
    day  = 365.25*(time-year)
    mjd = (365.25*year+day) - 678942
    mjd_month = np.zeros(len(mjd))
    mjd_mdiff = np.ones(len(mjd)-1)*np.mean(np.diff(mjd))
    mjd_month[0] = mjd[0]
    mjd_month[1:] = mjd_month[0] + np.cumsum(mjd_mdiff)
    return mjd_month


def tstep2mjd(tstep):
    year = np.floor(tstep)
    day = 365.25 * (tstep - year)
    totday = year * 365.25 + day
    mjd = totday - 678942
    return mjd


def est_trend(time,height,NoiseModel='None',est_acc=False,AR=0,MA=0,SA=False,SSA=False,jumps=False,jumplist = np.zeros(1),cust_sample=False,sample_period=0.0,monthly_data=False):
    result = {}
    time = time.flatten()
    height = height.flatten()

   # Prepare time and determine sample time
    if monthly_data: # Time step is 30.4375
        mjd = time2mjd_month(time)
    else: # Custom
        mjd = time2mjd(time)

    # remove nans if existing
    acc_idx = (np.isfinite(time) & np.isfinite(height))
    mjd   = mjd[acc_idx]
    height = height[acc_idx]

    # Prepare jumps
    if jumps:
        jump_mjd = tstep2mjd(jumplist)

    # Set up file names
    hector_path = os.getenv('HOME')+'/Code/Hector/hector_2/'
    fn_config = hector_path+'estimatetrend.ctl'
    fn_input  = hector_path+'input.mom'
    fn_output = hector_path+'output.txt'

    # Write input file
    if cust_sample:
        headerline = 'sampling period '+str(sample_period)
        mjd = np.round(mjd,0)
    else:
        sample_period = np.min(np.diff(mjd))
        headerline = 'sampling period '+str(sample_period)

    if jumps:
        for jump in range(len(jumplist)):
            headerline = headerline +' \noffset ' + str(jump_mjd[jump]) +' 0'
    np.savetxt(fn_input, np.transpose([mjd,height]), fmt=['%.4f','%.4f'], header=headerline)

    # Define config file
    config_list = []
    config_list.append('DataFile input.mom\n')
    config_list.append('DataDirectory ./\n')
    config_list.append('OutputFile trend.out\n')
    config_list.append('interpolate no\n')
    config_list.append('firstdifference no\n')
    config_list.append('PhysicalUnit m\n')
    if est_acc:
        config_list.append('DegreePolynomial 2\n')
    else:
        config_list.append('DegreePolynomial 1\n')
    if SA:
        config_list.append('seasonalsignal yes\n')
    else:
        config_list.append('seasonalsignal no\n')
    if SSA:
        config_list.append('halfseasonalsignal yes\n')
    else:
        config_list.append('halfseasonalsignal no\n')
    if jumps:
        config_list.append('estimateoffsets yes\n')
    else:
        config_list.append('estimateoffsets no\n')
    if NoiseModel=='Powerlaw':
        config_list.append('NoiseModels Powerlaw White\n')
    elif NoiseModel=='GGM':
        config_list.append('NoiseModels GGM White\n')
    elif NoiseModel=='ARMA':
        config_list.append('NoiseModels ARMA White\n')
        config_list.append('AR_p '+str(int(AR))+'\n')
        config_list.append('MA_q '+str(int(MA))+'\n')
    elif NoiseModel == 'None':
        print('No noise model, assuming white noise')
        config_list.append('NoiseModels White\n')
    else:
        raise ValueError('Unknown noise model')

    # Write config file
    open(fn_config,'w+').writelines(config_list)

    # Compute trend
    return_dir = os.getcwd()
    os.chdir(hector_path)
    os.system('./estimatetrend >'+fn_output)
    os.chdir(return_dir)

    # Read file and save data
    output_data = [line.rstrip('\n') for line in open(fn_output)]
    trend_str = next((s for s in output_data if 'trend:' in s), None)
    result['trend'] = np.zeros(2)
    result['trend'][0] = float(trend_str.split()[1])
    result['trend'][1] = float(trend_str.split()[3])
    if est_acc:
        acc_str = next((s for s in output_data if 'quadratic' in s), None)
        result['acc'] = np.zeros(2)
        result['acc'][0] = 2*float(acc_str.split()[3])
        result['acc'][1] = 2*float(acc_str.split()[5])
    if jumps:
        result['offset'] = np.zeros([len(jumplist), 3])
        str_found = (s for s in output_data if 'offset at' in s)
        for j in range(len(jumplist)):
            offset_str = next(str_found)
            result['offset'][j,0] = jumplist[j]
            result['offset'][j,1] = float(offset_str.split()[4])
            result['offset'][j,2] = float(offset_str.split()[6])
    return result

def generate_ggm_noise(time,height,num_ens):
    mjd = time2mjd_month(time)
    acc_idx = (np.isfinite(time) & np.isfinite(height))
    mjd = mjd[acc_idx]
    height_acc = height[acc_idx]

    # File names for Hector
    hector_path = os.getenv('HOME')+'/Code/Hector/hector_2/'

    fn_config = hector_path + 'estimatetrend.ctl'
    fn_input  = hector_path + 'input.mom'
    fn_output = hector_path + 'output.txt'

    sample_period = np.min(np.diff(mjd))
    headerline = 'sampling period ' + str(sample_period)
    np.savetxt(fn_input, np.transpose([mjd, height_acc]), fmt=['%.4f', '%.4f'], header=headerline)

    # Define config file
    config_list = []
    config_list.append('DataFile input.mom\n')
    config_list.append('DataDirectory ./\n')
    config_list.append('OutputFile trend.out\n')
    config_list.append('interpolate no\n')
    config_list.append('firstdifference no\n')
    config_list.append('PhysicalUnit m\n')
    config_list.append('DegreePolynomial 2\n')
    config_list.append('seasonalsignal no\n')
    config_list.append('halfseasonalsignal no\n')
    config_list.append('estimateoffsets no\n')
    config_list.append('NoiseModels GGM White\n')

    # Write config file
    open(fn_config, 'w+').writelines(config_list)

    # Compute trend
    return_dir = os.getcwd()
    os.chdir(hector_path)
    os.system('./estimatetrend >' + fn_output)
    os.chdir(return_dir)

    # Read trend data and noise parameters
    output_data = [line.rstrip('\n') for line in open(fn_output)]
    trend_str = next((s for s in output_data if 'trend:' in s), None)
    trend_orig = np.zeros(2)
    trend_orig[0] = float(trend_str.split()[1])
    trend_orig[1] = float(trend_str.split()[3])

    GGM_par = {}
    GGM_par['fraction'] = float((next((s for s in output_data if 'fraction  = ' in s), None)).split()[2])
    GGM_par['sigma']    = float((next((s for s in output_data if 'sigma     = ' in s), None)).split()[2])
    GGM_par['sigma_p']    = float((next((s for s in output_data if 'sigma     = ' in s), None)).split()[3][5:])

    GGM_par['d']        = float((next((s for s in output_data if 'd         = ' in s), None)).split()[2])
    GGM_par['kappa']        = float((next((s for s in output_data if 'kappa     = ' in s), None)).split()[2])
    GGM_par['1-phi']        = float((next((s for s in output_data if '1-phi     = ' in s), None)).split()[2])

    White_par = {}
    str_found = (s for s in output_data if 'fraction  =' in s)
    next(str_found)
    White_par['fraction'] = float((next(str_found)).split()[2])
    str_found = (s for s in output_data if 'sigma     =' in s)
    next(str_found)
    White_par['sigma'] = float((next(str_found)).split()[2])
    White_par['Driving']        = float((next((s for s in output_data if 'STD of the driving noise:' in s), None)).split()[5])

    fn_config = hector_path + 'simulatenoise.ctl'

    # Generate noise
    config_list = []
    config_list.append('SimulationDir .\n')
    config_list.append('SimulationLabel ggm\n')
    config_list.append('NumberOfSimulations '+str(1)+'\n')
    config_list.append('NumberOfPoints '+str(num_ens*len(time))+'\n')
    config_list.append('SamplingPeriod 30\n')
    config_list.append('TimeNoiseStart 300\n')
    config_list.append('NoiseModels GGM White\n')
    config_list.append('PhysicalUnit mm\n')
    open(fn_config,'w+').writelines(config_list)

    run_command = " printf '"+"{:.6f}".format(White_par['Driving'])+"\n"+"{:.6f}".format(GGM_par['fraction'])+"\n"+"{:.6f}".format(White_par['fraction'])+"\n"+"{:.6f}".format(GGM_par['d'])+"\n"+"{:.6f}".format(GGM_par['1-phi'])+"\n '| ./simulatenoise"
    return_dir = os.getcwd()
    os.chdir(hector_path)
    os.system(run_command)
    os.chdir(return_dir)

    height_perturbed = np.loadtxt(hector_path + 'ggm_0.mom',usecols=0).reshape([num_ens,len(time)])
    height_perturbed[:,~acc_idx] = np.nan # Noise at accepted time steps
    return height_perturbed

