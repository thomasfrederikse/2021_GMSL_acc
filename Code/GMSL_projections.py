# ------------------------------------------------------------------------
# Generate projections and hindcasts of sea level based on altimetry GMSL
#
# Part of "Extrapolating Empirical Models of Satellite-Observed Global Mean
# Sea Level to Estimate Future Sea Level Change
# R. S. Nerem, T. Frederikse, and B. D. Hamlington
#
# This code does the following:
# 1. Extract trend and acceleration from GSFC GMSL curve
# 2. Generate an ensemble using uncertainties on
#    - GIA
#    - Measurement uncertainties
#    - Serially-correlated noise
# 3. Extrapolate each ensemble member
#
# The research was carried out in part at the Jet Propulsion Laboratory,
# California Institute of Technology, under a contract with the National
# Aeronautics and Space Administration
#
# (c) All rights reserved
# -------------------------------------------------------------------------

import numpy as np
import os
from scipy.interpolate import interp1d
import scipy.io as scio
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import module_hector as hector


def main():
    settings = def_settings()
    gsl = read_gsl(settings)
    gia_ensemble = read_gia(settings)
    enso_corr_ens, pinatubo_corr_ens = read_enso_pinatubo(settings)
    ptb_gsl_alt, ptb_gsl_ggm = compute_gsl_ensemble(gsl, enso_corr_ens, pinatubo_corr_ens, settings)
    obs_ensemble, trend_ensemble, recon_ensemble, extrap_ensemble = compute_trend_acc_extrap_ensemble(gsl, ptb_gsl_alt, ptb_gsl_ggm, gia_ensemble, enso_corr_ens, pinatubo_corr_ens, settings)
    gmsl_f2020, gmsl_srocc = read_f2020_srocc(extrap_ensemble, settings)
    trend_stats, tseries_stats = compute_statistics(obs_ensemble, trend_ensemble, recon_ensemble, extrap_ensemble, settings)

    plot_figure_1_gmsl(tseries_stats, settings)
    plot_figure_3_histograms(trend_ensemble, extrap_ensemble, settings)
    plot_figure_4_extrap_vs_tg(tseries_stats, gmsl_f2020, settings)
    plot_figure_5_extrap_vs_srocc(tseries_stats, gmsl_srocc, settings)
    return


def def_settings():
    settings = {}
    # Directories
    settings["dir_project"] = os.getenv('HOME')+'/Projects/2021_GMSL_acc/'
    settings["dir_data"] = settings["dir_project"] + 'Data/'
    settings["dir_figures"] = settings["dir_project"] + 'Figures/'

    # Input files
    settings["fn_gmsl_gia"] = settings["dir_data"] + 'gmsl_gia.txt'
    settings['fn_gmsl'] = settings["dir_data"] + 'GMSL_TPJAOS_5.0_199209_202011.txt'
    settings['fn_enso'] = settings["dir_data"] + 'enso_correction.mat'
    settings['fn_pinatubo'] = settings["dir_data"] + 'gmsla_le_1991_fit.txt'
    settings['fn_SROCC_rcp26'] = settings["dir_data"] + 'SROCC_gmsl_26.txt'
    settings['fn_SROCC_rcp45'] = settings["dir_data"] + 'SROCC_gmsl_45.txt'
    settings['fn_SROCC_rcp85'] = settings["dir_data"] + 'SROCC_gmsl_85.txt'
    settings['fn_f2020'] = settings["dir_data"] + 'global_timeseries_measures.nc'

    # Define time steps
    settings['time_sat'] = np.arange(1993+1/24, 2020+23/24, 1/12)
    settings['time_extrap'] = np.arange(1960+1/24, 2101+1/24, 1/12)

    # Number of ensembles for Monte Carlo simulation
    settings['num_ens'] = 100000
    return settings


def read_gsl(settings):
    # Read global-mean geocentric sea level from GSFC:
    #
    # Beckley, B. D., Callahan, P. S., Hancock, D. W., Mitchum, G. T., & Ray, R. D. (2017).
    # On the “Cal-Mode” Correction to TOPEX Satellite Altimetry and Its Effect on the
    # Global Mean Sea Level Time Series.
    # Journal of Geophysical Research: Oceans, 122(11), 8371–8384. https://doi.org/10.1002/2017JC013090
    #
    # Data: https://doi.org/10.5067/GMSLM-TJ142
    time_read = np.loadtxt(settings['fn_gmsl'], comments='HDR')[:, 2]
    gsl_read = np.loadtxt(settings['fn_gmsl'], comments='HDR')[:, 5]
    gsl = np.zeros(len(settings['time_sat']))
    for idx, t in enumerate(settings['time_sat']):
        acc_idx = (time_read >= t-1/24) & (time_read < t+1/24)
        gsl[idx] = gsl_read[acc_idx].mean()
    gsl = remove_seasonal(settings['time_sat'], gsl)
    return gsl


def read_gia(settings):
    # Read pre-processed GIA altimetry correction estimates from:
    #
    # Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., & Blewitt, G. (2018).
    # GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science.
    # Geophysical Research Letters, 45(5), 2203–2212. https://doi.org/10.1002/2017GL076644
    gmsl_gia_raw = np.loadtxt(settings["fn_gmsl_gia"], delimiter=';')
    probability = gmsl_gia_raw[:, 0]
    gia_ensemble_trend = gmsl_gia_raw[:, 1]
    # Compute GIA random samples using inverse transform sampling
    sortidx = np.argsort(gia_ensemble_trend)
    gia_ensemble_trend = gia_ensemble_trend[sortidx]
    probability = probability[sortidx]
    interpolant = interp1d(np.cumsum(probability), gia_ensemble_trend, kind='linear', bounds_error=False, fill_value=(gia_ensemble_trend[0], gia_ensemble_trend[-1]))
    gia_ensemble = interpolant(np.random.uniform(0, 1, settings['num_ens']))
    return gia_ensemble


def read_enso_pinatubo(settings):
    # Read corrections for ENSO and Pinatubo from
    #
    # Nerem, R. S., Beckley, B. D., Fasullo, J. T., Hamlington, B. D., Masters, D., & Mitchum, G. T. (2018).
    # Climate-change–driven accelerated sea-level rise detected in the altimeter era.
    # Proceedings of the National Academy of Sciences, 115(9), 2022–2025. https://doi.org/10.1073/pnas.1717312115
    enso_file = scio.loadmat(settings['fn_enso'])
    enso_time = enso_file['time'].squeeze()
    enso_gmsl = enso_file['gmsl'].squeeze()
    enso_corr = np.interp(settings['time_sat'], enso_time, enso_gmsl)
    pinatubo_file = np.loadtxt(settings['fn_pinatubo'])
    pinatubo_corr = interp1d(pinatubo_file[:, 0], pinatubo_file[:, 1], kind='slinear', fill_value="extrapolate")(settings['time_sat'])
    # Generate ensemble by assuming an acceleration uncertainty of 0.01 mmyr^-2 in both corrections
    enso_corr_ens = (enso_corr + 0.5 * np.random.normal(0, 0.01, settings['num_ens'])[:, np.newaxis] * ((settings['time_sat'] - settings['time_sat'].mean())**2)[np.newaxis, :]).astype(np.float32)
    pinatubo_corr_ens = (pinatubo_corr + 0.5 * np.random.normal(0, 0.01, settings['num_ens'])[:, np.newaxis] * ((settings['time_sat'] - settings['time_sat'].mean())**2)[np.newaxis, :]).astype(np.float32)
    return enso_corr_ens, pinatubo_corr_ens


def read_f2020_srocc(extrap_ensemble, settings):
    # Read Frederikse et al. (2020) GMSL recontsruction and SROCC GMSL scenarios
    gmsl_f2020 = {}
    fh = Dataset(settings['fn_f2020'], 'r')
    fh.set_auto_mask(False)
    gmsl_f2020['time'] = fh['time'][:] / 365.25 + 1900.049
    gmsl_f2020['gmsl'] = np.vstack([fh['global_average_sea_level_change_lower'][:], fh['global_average_sea_level_change'][:], fh['global_average_sea_level_change_upper'][:]]).T

    gmsl_srocc = {}
    rcp26 = np.loadtxt(settings['fn_SROCC_rcp26'], skiprows=1, max_rows=95)
    rcp45 = np.loadtxt(settings['fn_SROCC_rcp45'], skiprows=1, max_rows=95)
    rcp85 = np.loadtxt(settings['fn_SROCC_rcp85'], skiprows=1, max_rows=95)

    bline = extrap_ensemble[:, (settings['time_extrap'] > 1985) & (settings['time_extrap'] < 2006)].mean()
    rcp26[:, 1:] *= 1000
    rcp45[:, 1:] *= 1000
    rcp85[:, 1:] *= 1000
    rcp26[:, 1:] += bline
    rcp45[:, 1:] += bline
    rcp85[:, 1:] += bline

    gmsl_srocc['time'] = rcp26[:, 0] + 0.5
    gmsl_srocc['gmsl_rcp26'] = rcp26[:, 1:]
    gmsl_srocc['gmsl_rcp45'] = rcp45[:, 1:]
    gmsl_srocc['gmsl_rcp85'] = rcp85[:, 1:]
    return gmsl_f2020, gmsl_srocc


def compute_gsl_ensemble(gsl,enso_corr_ens,pinatubo_corr_ens,settings):
    # Generate an ensemble of global-mean geocentric sea level curve
    # uncertainties around the mean, corrected for ENSO/Pinatubo

    # Term 1: Ablain et al: measurement uncertainties
    t_rand_highf  = gen_autocov(1.7,1.5, 1.2, 2/12, settings)
    t_rand_medf   = gen_autocov(1.3,1.2, 1, 1, settings)
    t_rand_wettr  = gen_autocov(1.1,1.1, 1.1, 5, settings)
    t_rand_lorbit = gen_autocov(1.12,0.5, 0.5, 10, settings)
    t_rand_intmis = gen_randjump(settings)
    t_rand_dorbit = gen_drift(0.1,settings)
    t_rand_tpx = gen_tpxdrift(settings)
    ptb_gsl_alt  = (t_rand_highf + t_rand_medf + t_rand_wettr + t_rand_lorbit + t_rand_intmis + t_rand_dorbit+t_rand_tpx).astype(np.float32)

    # Term 2:  Serially-correlated noise
    ptb_gsl_ggm = hector.generate_ggm_noise(settings['time_sat'], gsl-enso_corr_ens.mean(axis=0)-pinatubo_corr_ens.mean(axis=0), settings['num_ens'])
    return ptb_gsl_alt,ptb_gsl_ggm


def compute_trend_acc_extrap_ensemble(gsl,ptb_gsl_alt,ptb_gsl_ggm,gia_ensemble,enso_corr_ens,pinatubo_corr_ens,settings):
    # Compute the trend and acceleration and the resulting extrapolation

    # Create an ensemble of GMSL estimates: GMSL = GSL + GIA - Pinatubo - ENSO + pertubations
    obs_ensemble = gsl[np.newaxis,:] + ptb_gsl_alt + ptb_gsl_ggm - enso_corr_ens - pinatubo_corr_ens - gia_ensemble[:,np.newaxis] * (settings['time_sat']-settings['time_sat'].mean())

    # Remove 2000-2020 mean
    baseline_idx = (settings['time_sat']>2000)&(settings['time_sat']<2020)
    obs_ensemble -= obs_ensemble[:,baseline_idx].mean(axis=1)[:,np.newaxis]

    # Estimate the trend and acceration and
    trend_ensemble   = np.zeros([settings['num_ens'],3])
    recon_ensemble   = np.zeros([settings['num_ens'],len(settings['time_sat'])])
    extrap_ensemble  = np.zeros([settings['num_ens'],len(settings['time_extrap'])])
    amat = np.ones([len(settings['time_sat']),3])
    amat[:,1] = settings['time_sat'] - np.mean(settings['time_sat'])
    amat[:,2] = 0.5*(settings['time_sat'] - np.mean(settings['time_sat']))**2

    amat_extrap = np.ones([len(settings['time_extrap']),3])
    amat_extrap[:,1] = settings['time_extrap'] - np.mean(settings['time_sat'])
    amat_extrap[:,2] = 0.5*(settings['time_extrap'] - np.mean(settings['time_sat']))**2

    amat_T  = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))

    for i in range(settings['num_ens']):
        trend_ensemble[i,:] = np.dot(amat_sq, np.dot(amat_T, obs_ensemble[i,:]))
        recon_ensemble[i,:] = amat@trend_ensemble[i,:]
        extrap_ensemble[i, :] = amat_extrap @ trend_ensemble[i, :]
    return(obs_ensemble, trend_ensemble, recon_ensemble, extrap_ensemble)


def compute_statistics(obs_ensemble, trend_ensemble, recon_ensemble, extrap_ensemble, settings):
    # Compute statistics

    # Trends, accelerations, GMSL in 2050 and 2100
    acc_idx_2050 = (np.floor(settings['time_extrap']) == 2050)
    acc_idx_2100 = (np.floor(settings['time_extrap']) == 2100)

    trend_stats = {}
    trend_stats['percentiles'] = [5, 17, 50, 83, 95]
    trend_stats["trend"] = np.percentile(trend_ensemble[:, 1], trend_stats['percentiles'])
    trend_stats["accel"] = np.percentile(trend_ensemble[:, 2], trend_stats['percentiles'])
    trend_stats["gsml_2050"] = np.percentile(extrap_ensemble[:, acc_idx_2050].mean(axis=1) / 10, trend_stats['percentiles'])  # Units cm
    trend_stats["gsml_2100"] = np.percentile(extrap_ensemble[:, acc_idx_2100].mean(axis=1) / 10, trend_stats['percentiles'])  # Units cm

    # Time series
    tseries_stats = {}
    tseries_stats['percentiles'] = [5,50,95]
    tseries_stats['obs'] = np.percentile(obs_ensemble,tseries_stats['percentiles'],axis=0)
    tseries_stats['recon'] = np.percentile(recon_ensemble,tseries_stats['percentiles'],axis=0)
    tseries_stats['extrap'] = np.percentile(extrap_ensemble,tseries_stats['percentiles'],axis=0)
    return(trend_stats,tseries_stats)


def plot_figure_1_gmsl(tseries_stats,settings):
    plt.figure(figsize=(12,4))
    plt.subplot(1,2,1)
    plt.fill_between(settings['time_sat'],tseries_stats['obs'][0,:],tseries_stats['obs'][2,:],color='C0',alpha=0.4,linewidth=0)
    plt.fill_between(settings['time_sat'],tseries_stats['recon'][0,:],tseries_stats['recon'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.plot(settings['time_sat'],tseries_stats['recon'][1,:],'k',linewidth=2,label='Quadratic fit')
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL')
    plt.grid()
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.xlim([1993,2020.7])
    plt.ylim([-65,50])
    plt.ylabel('Sea level (mm)',fontsize=9)

    plt.subplot(1,2,2)
    plt.fill_between(settings['time_sat'],tseries_stats['obs'][0,:],tseries_stats['obs'][2,:],color='C0',alpha=0.4,linewidth=0)
    plt.fill_between(settings['time_extrap'],tseries_stats['extrap'][0,:],tseries_stats['extrap'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.plot(settings['time_extrap'],tseries_stats['extrap'][1,:],'k',linewidth=2,label='Extrapolation using quadratic fit')
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL')
    plt.grid()
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.legend(fontsize=9,loc=4)
    plt.xlim([1993,2051])
    plt.ylim([-70, 290])
    plt.ylabel('Sea level (mm)',fontsize=9)
    plt.tight_layout()
    plt.savefig(settings["dir_figures"]+'figure_1_gmsl.png')
    return

def plot_figure_3_histograms(trend_ensemble,extrap_ensemble,settings):
    # Plot histograms
    acc_idx_2050 = (np.floor(settings['time_extrap']) == 2050)

    plt.figure(figsize=(12,4))
    # Panel 1: trend histogram
    plt.subplot(1,3,1)
    plt.hist(trend_ensemble[:,1],np.arange(2.3,4,0.02),alpha=0.5)
    pct_bnd = [2.84, 3.52]
    acc_hst = (trend_ensemble[:, 1] > pct_bnd[0]) & (trend_ensemble[:, 1] < pct_bnd[1])
    plt.hist(trend_ensemble[acc_hst,1],np.arange(2.3,4,0.02),color='C0')
    plt.xlim([2.4, 3.9])
    plt.ylim([0, 4050])
    plt.xlabel('Trend (mm yr$^{-1}$)',fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(np.arange(0,4000,1000),[],fontsize=9)
    plt.grid()

    # Panel 2: acceleration histogram
    plt.subplot(1,3,2)
    plt.hist(trend_ensemble[:,2],np.arange(-0.1,0.3,0.004),alpha=0.5,color='C0')
    pct_bnd = [0.012, 0.148]
    acc_hst = (trend_ensemble[:, 2] > pct_bnd[0]) & (trend_ensemble[:, 2] < pct_bnd[1])
    plt.hist(trend_ensemble[acc_hst,2],np.arange(-0.1,0.3,0.004),color='C0')
    plt.xlabel('Acceleration (mm yr$^{-2}$)',fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(np.arange(0,4000,1000),[],fontsize=9)
    plt.xlim([-0.08, 0.23])
    plt.ylim([0, 4050])
    plt.grid()

    # Panel 3: GMSL in 2050 histogram
    plt.subplot(1,3,3)
    plt.hist(extrap_ensemble[:,acc_idx_2050].mean(axis=1)/10,np.arange(0,40,0.35),alpha=0.5,color='C0')
    pct_bnd = [14.7, 25.9]
    acc_hst = (extrap_ensemble[:, acc_idx_2050].mean(axis=1)/10 > pct_bnd[0]) & (extrap_ensemble[:, acc_idx_2050].mean(axis=1)/10 < pct_bnd[1])
    plt.hist((extrap_ensemble[:,acc_idx_2050].mean(axis=1)/10)[acc_hst],np.arange(0,40,0.35),color='C0')
    plt.xlabel('Sea level in 2050 (cm)',fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(np.arange(0,4000,1000),[],fontsize=9)
    plt.xlim([8, 33])
    plt.ylim([0, 4050])
    plt.grid()

    plt.tight_layout()
    plt.savefig(settings["dir_figures"]+'figure_3_histograms.png')
    return


def plot_figure_4_extrap_vs_tg(tseries_stats,gmsl_f2020,settings):
    plt.figure(figsize=(6,4))
    plt.fill_between(settings['time_extrap'],tseries_stats['extrap'][0,:],tseries_stats['extrap'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.fill_between(gmsl_f2020['time'],gmsl_f2020['gmsl'][:,0],gmsl_f2020['gmsl'][:,2],color='C1',alpha=0.3,linewidth=0)
    plt.plot(gmsl_f2020['time'],gmsl_f2020['gmsl'][:,1],'C1',linewidth=2,label='Observed GMSL (tide gauges)')
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL (altimetry)')
    plt.plot(settings['time_extrap'],tseries_stats['extrap'][1,:],'k',linewidth=2,label='Quadratic fit extrapolated')
    plt.ylabel('Sea level (mm)',fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.legend(fontsize=9,loc=2)
    plt.xlim([1940,2021])
    plt.ylim([-150,50])
    plt.grid()
    plt.tight_layout()
    plt.savefig(settings["dir_figures"]+'figure_4_extrap_vs_tg.png')

def plot_figure_5_extrap_vs_srocc(tseries_stats,gmsl_srocc,settings):
    plt.figure(figsize=(12,4))
    plt.subplot(1,3,1)
    plt.fill_between(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp26'][:,0],gmsl_srocc['gmsl_rcp26'][:,2],color='C1',alpha=0.4,linewidth=0.5)
    plt.fill_between(settings['time_extrap'],tseries_stats['extrap'][0,:],tseries_stats['extrap'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL')
    plt.plot(settings['time_extrap'],tseries_stats['extrap'][1,:],'k',linewidth=2,label='Quadratic fit extrapolated')
    plt.plot(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp26'][:,1],color='C1',linewidth=2,label='Model projection (SROCC)')
    plt.ylabel('Sea level (mm)',fontsize=9)
    plt.xlim([1990,2051])
    plt.ylim([-80,320])
    plt.grid()
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.legend(fontsize=9,loc=2)
    plt.title('SROCC RCP2.6',fontsize=9)

    plt.legend(fontsize=9,loc=2)
    plt.subplot(1,3,2)
    plt.fill_between(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp45'][:,0],gmsl_srocc['gmsl_rcp45'][:,2],color='C1',alpha=0.4,linewidth=0.5)
    plt.fill_between(settings['time_extrap'],tseries_stats['extrap'][0,:],tseries_stats['extrap'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL')
    plt.plot(settings['time_extrap'],tseries_stats['extrap'][1,:],'k',linewidth=2,label='Quadratic fit extrapolated')
    plt.plot(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp45'][:,1],color='C1',linewidth=2,label='Model projection (SROCC)')
    plt.ylabel('Sea level (mm)',fontsize=9)
    plt.xlim([1990,2051])
    plt.ylim([-80,320])
    plt.grid()
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.title('SROCC RCP4.5',fontsize=9)

    plt.subplot(1,3,3)
    plt.fill_between(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp85'][:,0],gmsl_srocc['gmsl_rcp85'][:,2],color='C1',alpha=0.4,linewidth=0.5)
    plt.fill_between(settings['time_extrap'],tseries_stats['extrap'][0,:],tseries_stats['extrap'][2,:],color='k',alpha=0.3,linewidth=0)
    plt.plot(settings['time_sat'],tseries_stats['obs'][1,:],'C0',linewidth=2,label='Observed GMSL')
    plt.plot(settings['time_extrap'],tseries_stats['extrap'][1,:],'k',linewidth=2,label='Quadratic fit extrapolated')
    plt.plot(gmsl_srocc['time'],gmsl_srocc['gmsl_rcp85'][:,1],color='C1',linewidth=2,label='Model projection (SROCC)')
    plt.ylabel('Sea level (mm)',fontsize=9)
    plt.xlim([1990,2051])
    plt.ylim([-80,320])
    plt.grid()
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.title('SROCC RCP8.5',fontsize=9)
    plt.tight_layout()
    plt.savefig(settings["dir_figures"]+'figure_5_extrap_vs_srocc.png')


# -------------------------------------------------------------------
# Ablain et al. (2018) functions
# These functions are used to generate an ensemble of GMSL
# estimates with the uncertainty structure from Ablain et al. (2018)
# -------------------------------------------------------------------
def gen_randjump(settings):
    tjump = [(1999 + 1.5 / 12),(2002 + 4.5 / 12),(2008 + 8.5 / 12),(2016 + 8.5 / 12)]
    hjump = [2,0.5,0.5,0.5]
    t_rand = np.zeros([settings['num_ens'],len(settings['time_sat'])])
    for jump in range(len(tjump)):
        rndjump = np.random.randn(settings['num_ens'])*hjump[jump]
        t_acc= settings['time_sat']>tjump[jump]
        t_rand[:,t_acc] += rndjump[:,np.newaxis]
    return(t_rand)

def gen_drift(drift,settings):
    t_rand = drift * np.random.randn(settings['num_ens'])[:,np.newaxis] * (settings['time_sat']-settings['time_sat'].mean())[np.newaxis, :]
    return(t_rand)

def gen_tpxdrift(settings):
    t_rand = np.zeros([settings['num_ens'], len(settings['time_sat'])])
    t_txa = settings['time_sat'] < (1999 + 1.5 / 12)
    t_txb = (settings['time_sat'] > (1999 + 1.5 / 12)) & (settings['time_sat']<(2002 + 4.5 / 12))
    t_rand[:, t_txa] = (np.random.randn(settings['num_ens']) * 0.7/12)[:,np.newaxis]
    t_rand[:, t_txb] = (np.random.randn(settings['num_ens']) * 0.1/12)[:,np.newaxis]
    t_rand = np.cumsum(t_rand, axis=1)
    return t_rand


def gen_autocov(sig_TPX, sig_J1, sig_J23, l_factor, settings):
    jump_0 = np.argmin(np.abs(settings['time_sat']-2002+4.5/12))
    jump_1 = np.argmin(np.abs(settings['time_sat']-2008+8.5/12))
    sig_array = np.zeros(len(settings['time_sat']))
    sig_array[:jump_0] = sig_TPX
    sig_array[jump_0:jump_1] = sig_J1
    sig_array[jump_1:] = sig_J23
    t_distance = np.abs(settings['time_sat'][:,np.newaxis] - settings['time_sat'][np.newaxis,:])
    covmat = np.exp(-0.5*(t_distance/l_factor)**2) + np.eye(len(settings['time_sat']))*1.0e-10
    covc = np.linalg.cholesky(covmat)

    t_rand = np.zeros([settings['num_ens'], len(settings['time_sat'])])
    for n in range(settings['num_ens']):
        t_rand[n, :] = sig_array*np.matmul(covc, np.random.randn(len(settings['time_sat'])))
    return t_rand


# -----------------
# Generic functions
# -----------------
def remove_seasonal(time, height):
    # Remove the mean seasonal cycle from a time series
    # use linear least squares
    acc_idx = np.isfinite(height)
    amat = np.ones([len(time[acc_idx]),6])
    amat[:, 1] = np.sin(2*np.pi*time[acc_idx])
    amat[:, 2] = np.cos(2*np.pi*time[acc_idx])
    amat[:, 3] = np.sin(4*np.pi*time[acc_idx])
    amat[:, 4] = np.cos(4*np.pi*time[acc_idx])
    amat[:, 5] = time[acc_idx] - time[acc_idx].mean()
    sol = np.linalg.lstsq(amat, height[acc_idx], rcond=None)[0]
    sol[0] = 0
    sol[5] = 0
    height_deseas = np.zeros(len(time))*np.nan
    height_deseas[acc_idx] = height[acc_idx] - np.matmul(amat, sol)
    return height_deseas

if __name__ == '__main__':
    main()
