# Extrapolating Empirical Models of Satellite-Observed Global Mean Sea Level to Estimate Future Sea Level Change
This repository contains the codes and input data to reproduce the results from the following paper:

R. S. Nerem, T. Frederikse, and B. D. Hamlington: Extrapolating Empirical Models of Satellite-Observed Global Mean Sea Level to Estimate Future Sea Level Change, Earth's Future, 2022

(c) California 

## Running the analysis
To reproduce the figures from this paper, we need the following:
- The Hector software (I've tested version 1.9 and version 2, and both should work). Download the source code from http://segal.ubi.pt/hector/ and compile the program. Then, go to the file `module_hector.py` in this repository, and change the `hector_path` variable to point to the directory where the Hector executables can be found. 
- A Python installation with the following packages installed: `SciPy`, `NumPy`, `matplotlib`, `netCDF`. 
- In the file `GMSL_projections.py`, update the `settings["dir_project"]` variable to point to the directory where this repository can be found. 

Then run `GMSL_projections.py`, either interactively in IPython or using `python GMSL_projections.py`. This will produce the figures from the paper and save them into the `Figures` directory. 

NOTE: With the standard settings, running the analysis consumes quite some memory. Reduce the parameter `settings['num_ens']` to reduce the memory usage. 

## Source data
The directory `Data` contains the source data used in the analysis, with the following files:
- `GMSL_TPJAOS_5.0_199209_202011.txt`: Global-mean sea-level curve from NASA GSFC. The data can be downloaded from PO.DAAC: Beckley, B.; Zelensky, N.P.; Holmes, S.A.;Lemoine, F.G.; Ray, R.D.; Mitchum, G.T.; Desai, S.; Brown, S.T.. 2016. Global Mean Sea Level Trend from Integrated Multi-Mission Ocean Altimeters TOPEX/Poseidon Jason-1 and OSTM/Jason-2 Version 4.2. Ver. 4.2. PO.DAAC, CA, USA. Dataset accessed 2022-02-22 at http://dx.doi.org/10.5067/GMSLM-TJ142. 
The data is described in Beckley, B. D., Callahan, P. S., Hancock, D. W., Mitchum, G. T., & Ray, R. D. (2017). On the “Cal-Mode” Correction to TOPEX Satellite Altimetry and Its Effect on the Global Mean Sea Level Time Series. Journal of Geophysical Research: Oceans, 122(11), 8371–8384. https://doi.org/10.1002/2017JC013090
- `SROCC_gmsl_26.txt` GMSL projections from the IPCC Special Report on the Ocean and Cryosphere in a Changing Climate: Oppenheimer, M., Glavovic, B., Hinkel, J., van de Wal, R. S. W., Magnan, A. K., Abd-Elgawad, A., Cai, R., Cifuentes-Jar, M., DeConto, R. M., Ghosh, T., Hay, J., Isla, F., Marzeion, B., Meyssignac, B., & Sebesvari, B. (2019). Sea Level Rise and Implications for Low Lying Islands, Coasts, and Communities. In H.-O. Pörtner, D. C. Roberts, V. Masson-Delmotte, P. Zhai, M. Tignor, E. Poloczanska, K. Mintenbeck, A. Alegría, M. Nicolai, A. Okem, J. Petzold, B. Rama, & N. M. Weyer (Eds.), IPCC Special Report on the Ocean and Cryosphere in a Changing Climate. IPCC. https://www.ipcc.ch/srocc. This file corresponds to the RCP2.6 scenario
- `SROCC_gmsl_45.txt` As above, but for the RCP4.5 scenario.
- `SROCC_gmsl_85.txt` As above, but for the RCP8.5 scenario.
- `enso_correction.mat`: Estimate of GMSL variability due to ENSO, as described in Hamlington, B. D., Frederikse, T., Nerem, R. S., Fasullo, J. T., & Adhikari, S. (2020). Investigating the Acceleration of Regional Sea‐level Rise During the Satellite Altimeter Era. Geophysical Research Letters. https://doi.org/10.1029/2019GL086528
- `global_timeseries_measures.nc`: Estimate of 20th-century GMSL changes from Frederikse, T.; Landerer, F.; Caron, L.; Adhikari, S.; Parkes, D.; Humphrey V.; Dangendorf, S.; Hogarth, P.; Zanna, L.; Cheng, L.; Wu, Y.. 2020. Reconstructed Global Mean Sea Level from GRACE and In Situ 1900 to 2018. Ver. 1.0. PO.DAAC, CA, USA. Dataset accessed 2022-02-22 at https://doi.org/10.5067/GMSLT-FJPL1. A description can be found in Frederikse, T., Landerer, F., Caron, L., Adhikari, S., Parkes, D., Humphrey, V. W., Dangendorf, S., Hogarth, P., Zanna, L., Cheng, L., & Wu, Y.-H. (2020). The causes of sea-level rise since 1900. Nature, 584(7821), 393–397. https://doi.org/10.1038/s41586-020-2591-3
- `gmsl_gia.txt`: Ensemble of estimates of GMSL changes due to GIA. See for details Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., & Blewitt, G. (2018). GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science. Geophysical Research Letters, 45(5), 2203–2212. https://doi.org/10.1002/2017GL076644
- `gmsla_le_1991_fit.txt`: Estimate of the effects of the Pinatubo eruption on GMSL. See for more details: Fasullo, J. T., Nerem, R. S., & Hamlington, B. (2016). Is the detection of accelerated sea level rise imminent? Scientific Reports, 6(1). https://doi.org/10.1038/srep31245
 
