import collections

'''
This script contains an ordered dictionary of the key words and their sample values. These values
were taken from a Fixed-Slit file.

* The information was obtained from emails in the pipeline group and the information contained in the
online dictionary located at:
https://iwjwdmsdauiwebv.stsci.edu/portal/Mashup/Clients/jwkeywords/

'''

# define dictionary
keywd_dict = collections.OrderedDict()

# Standard parameters
keywd_dict['SIMPLE']  = 'T' # Written by IDL
keywd_dict['BITPIX']  = 8 # Bits per data value, e.g. 8
keywd_dict['NAXIS']   = 0 # Number of data array dimensions,  e.g. 0
keywd_dict['EXTEND']  = 'T' # File may contain standard extensions
keywd_dict['NEXTEND'] = 3 # Number of standard extensions, e.g. 3

keywd_dict['TELESCOP']= 'JWST' # Telescope used to acquire data
keywd_dict['INSTRUME']= 'NIRSPEC' # Identifier for instrument used to acquire data
keywd_dict['RADESYS'] = 'ICRS' # Coordinate reference frame for RA and Dec

keywd_dict['DATE']    = '2016-07-12T12:43:27' # UTC date file created, e.g. 2016-07-12T12:43:27
keywd_dict['ORIGIN']  = 'STScI' # Institution responsible for creating FITS file
keywd_dict['FILENAME']= 'jwtest1001001_01101_00001_NRS1_uncal.fits' # Name of file, e.g. jwtest1001001_01101_00001_NRS1_uncal.fits
keywd_dict['FILETYPE']= 'UNCALIBRATED' # Type of data found in file
keywd_dict['DPSW_VER']= '0.1.1' # Data processing software version number, e.g. 0.1.1

# Observation identifiers
keywd_dict['DATE-OBS']= '2013-01-19' # UTC date at start of exposure, e.g. 2013-01-19
keywd_dict['TIME-OBS']= '18:23:34.230' # UTC time at start of exposure, e.g. 18:23:34.230
keywd_dict['OBS_ID']  = 'SLIT-COMBO-001' # full programmatic observation identifier
keywd_dict['VISIT_ID']= '' # visit identifier
keywd_dict['PROGRAM'] = 'test1' # program number, e.g. test1
keywd_dict['OBSERVTN']= '001' # observation number, e.g. 001
keywd_dict['VISIT']   = '001' # visit number, e.g. 001
keywd_dict['VISITGRP']= '01' # visit group identifier, e.g. 01
keywd_dict['SEQ_ID']  = 1 # parallel sequence identifier, e.g. 1
keywd_dict['ACT_ID']  = '01' # activity identifier, e.g. 01
keywd_dict['EXPOSURE']= '00001' # exposure request number, e.g. 00001

# Exposure parameters
keywd_dict['DETECTOR']= 'NRS1' # name of detector used to acquire data
keywd_dict['NINTS']   = 1 # number of integrations within exposure, e.g. 1
keywd_dict['NGROUPS'] = 20 # number of groups within integration, e.g. 20
keywd_dict['ZEROFRAM']= 'F' # boolean, T if a zero frame was read separately
keywd_dict['READPATT']= 'NRSRAPID' # readout pattern
keywd_dict['DATAPROB']= 'F' # boolean, T if science telemetry indicated any problems

# Program information
keywd_dict['TITLE']   = 'proposal_title1' # proposal title
keywd_dict['PI_NAME'] = 'UNKNOWN' # Name of principal investigator, e.g. UNKNOWN
keywd_dict['CATEGORY']= 'N/A' # program category, e.g. N/A
keywd_dict['SUBCAT']  = 'N/A' # program sub-category, e.g. N/A
keywd_dict['SCICAT']  = 'N/A' # science category assigned during TAC process, e.g. N/A
keywd_dict['CONT_ID'] = 'N/A' # continuation of the specified previous program, e.g. N/A

# Observation information
keywd_dict['TEMPLATE']= 'N/A' # proposal instruction template used, e.g. N/A
keywd_dict['OBSLABEL']= '#TODO' # proposer label for observation, e.g. #TODO

# Visit information
keywd_dict['VISITYPE']= 'PRIME' # type of visit (prime or parallel)
keywd_dict['VSTSTART']= '2013-01-19T18:27:22' # UTC visit start time, e.g. 2013-01-19T18:27:22
keywd_dict['WFSVISIT']= 'NO' # wavefront sensing and control visit indicator, options are:
#                                                    'NO', 'YES', 'SENSING_CONTROL', 'SENSING_ONLY'
keywd_dict['NEXPOSUR']= 1 # total number of exposures in visit, e.g. 1
keywd_dict['INTARGET']= 'T' # boolean, T if at least one exposure in visit is internal
keywd_dict['EXTARGET']= 'T' # boolean, T if at least one exposure in visit is external
keywd_dict['TARGOOPP']= 'F' # boolean, visit scheduled as target of opportunity

# Exposure information
keywd_dict['PNTG_SEQ']= 1 # pointing sequence number, e.g. 1
keywd_dict['EXPCOUNT']= 1 # count of the exposures in visit, e.g. 1
keywd_dict['EXP_TYPE']= 'NRS_FIXEDSLIT' # type of data in exposure, options are:
#                         'NRS_TASLIT', 'NRS_TACQ', 'NRS_TACONFIRM', 'NRS_CONFIRM', 'NRS_FIXEDSLIT', 'NRS_AUTOWAVE',
#                         'NRS_IFU', 'NRS_MSASPEC', ' NRS_AUTOFLAT', ' NRS_IMAGE', ' NRS_FOCUS', ' NRS_DARK',
#                         'NRS_LAMP', 'NRS_BOTA', 'NRS_BRIGHTOBJ'

# Target information
keywd_dict['TARGPROP']= 'UNKNOWN' # proposer's name for the target
keywd_dict['TARGNAME']= 'NGC 104' # standard astronomical catalog name for target, e.g. 'NGC 104 '
keywd_dict['TARGTYPE']= 'FIXED' # fixed target, moving target, or generic target, options are:
#                                                                      'FIXED', 'MOVING', 'GENERIC'
keywd_dict['TARG_RA'] = 0.0 # target RA computed at time of exposure, e.g. 0.0
keywd_dict['TARGURA'] = 0.0 # target RA uncertainty, e.g. 0.0
keywd_dict['TARG_DEC']= 0.0 # target DEC computed at time of exposure, e.g. 0.0
keywd_dict['TARRUDEC']= 0.0 # target Dec uncertainty, e.g. 0.0
keywd_dict['PROP_RA'] = 0.0 # proposer specified RA for the target, e.g. 0.0
keywd_dict['PROP_DEC']= 0.0 # proposer specified Dec for the target, e.g. 0.0
keywd_dict['PROPEPOC']= 2000.0 # proposer specified epoch for RA and Dec, e.g. 2000.0

# Exposure times
keywd_dict['EXPSTART']= 56311.76636840278 # UTC exposure start time (MJD), e.g. 56311.76636840278
keywd_dict['EXPMID']  = 56311.76636840278 # UTC exposure mid time (MJD), e.g. 56311.76636840278
keywd_dict['EXPEND']  = 56311.76763953704 # UTC exposure end time (MJD), e.g. 56311.76763953704

# Exposure time parameters
keywd_dict['NSAMPLES']= 1 # number of A/D samples per pixel, e.g. 1
keywd_dict['NFRAMES'] = 1 # number of frames in group, e.g. 1
keywd_dict['GROUPGAP']= 10 # number of frames dropped between groups, e.g. 10
keywd_dict['TSAMPLE'] = 10 # delta time between samples in microseconds, e.g. 10
keywd_dict['NRESET']  = 1 # number of resets between integrations, e.g. 1
keywd_dict['NRSTSTRT']= 1 # number of extra resets at start of exposure, e.g. 1
keywd_dict['TFRAME']  = 5.49132 # [seconds] time between frames, e.g. 5.49132
keywd_dict['TGROUP']  = 5.49132 # [seconds] time between groups, e.g 5.49132
keywd_dict['EFFINTTM']= 104.33508 # [seconds] effective integration time, e.g. 104.33508
keywd_dict['EFFEXPTM']= 104.33508 # [seconds] effective exposure time, e.g. 104.33508
keywd_dict['CHRGTIME']= -1.0 # [seconds] charge accumulation time, e.g. -1.0
keywd_dict['DURATION']= -1.0 # [seconds] total duration of exposuree.g. -1.0

# Subarray parameters
keywd_dict['SUBARRAY']= 'ALLSLITS' # name of subarray used, options are:
#                          '1024X16', '128X128', '128X2048', '2048X128', '2048X64', '32X32', '64X2048', '8X8',
#                          'ALLSLITS', 'BRIGHTSKY', 'FULL', 'GENERIC', 'MASK1065', 'MASK1140', 'MASK1550',
#                          'MASKLYOT', 'S1600A1', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'SLITLESSPRISM',
#                          'STRIPE', 'SUB1024A', 'SUB1024B', 'SUB128', 'SUB16', 'SUB160', 'SUB160P',
#                          'SUB1A', 'SUB1B', 'SUB2048', 'SUB256', 'SUB32', 'SUB320', 'SUB400P', 'SUB512',
#                          'SUB64', 'SUB640', 'SUB64P', 'SUB80', 'SUB96', 'SUBGRISM128', 'SUBGRISM256',
#                          'SUBGRISM64', 'SUBSTRIP256', 'SUBSTRIP96', 'SUBPRISM', 'WFSS128C', 'WFSS128R',
#                          'WFSS64C', 'WFSS64R', 'N/A']
keywd_dict['SUBSTRT1']= 1 # starting pixel number in the SIAS x direction, e.g. 1
keywd_dict['SUBSIZE1']= 2048 # number of pixels in the SIAS x direction, e.g. 2048
keywd_dict['SUBSTRT2']= 897 # starting pixel number in the SIAS y direction, e.g. 897
keywd_dict['SUBSIZE2']= 256 # number of pixels in the SIAS y direction, e.g. 256
keywd_dict['NXLIGHT'] = '#TODO' # number of light sensitive x values (columns), e.g #TODO
keywd_dict['FASTAXIS']= 2 # Direction of fast readout, options are: 1, 2, -1, -2
keywd_dict['SLOWAXIS']= 1 # Direction of slow readout, options are: 1, 2, -1, -2

# NIRSpec configuration (NIRSpec only)
keywd_dict['FILTER']  = 'F070LP' # name of the grating element used, options are:
#                                'CLEAR', 'F070LP', 'F100LP', 'F110W', 'F140X', 'F170LP', 'F290LP', 'P750L', 'NULL'
keywd_dict['GRATING'] = 'G140H' # name of grating used, options are:
#                                   'G140M', 'G235M', 'G395M', 'G140H', 'G235H', 'G395H', 'PRISM',
#                                   'MIRROR', 'NULL', 'N/A', 'ANY'
keywd_dict['GWAXTILT']= 0.35896975 # grating x tilt, e.g. 0.35896975
keywd_dict['GWAYTILT']= 0.13438272 # grating y tilt, e.g. 0.13438272
keywd_dict['FXD_SLIT']= 'NONE' # name of fixed slit aperture used, options are:
#                                   'NONE', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'S1600A1', 'NULL'
keywd_dict['MSASTATE']= 'CONFIGURED' # state of MSA, options are:
#                                       'CONFIGURED', 'LAUNCHLOCK_ALLCLOSED', 'PRIMARYPARK_ALLOPEN',
#                                       'PRIMARYPARK_ALLCLOSED', 'PRIMARYPARK_CONFIGURED'
keywd_dict['FOCUSPOS']= 0.0 # [mm] focus position for NIRSpec, e.g. 0.0

# NIRSpec MSA supporting files (NIRSpec MSA only)
keywd_dict['MSACONFG']= 'N/A' # MSA configuration file name, e.g. 'N/A'

# lamp configuration
keywd_dict['LAMP']    = 'ARGON' # internal lamp state, e.g. 'ARGON'

# Guide star information
keywd_dict['GS_ORDER']= 'N/A' # index of guide star, e.g. 'N/A'
keywd_dict['GSSTRTTM']= 'N/A' # UTC start time of guide star acquisition, e.g. 'N/A'
keywd_dict['GSENDTIM']= 'N/A' # UTC end time of guide star acquisition, e.g. 'N/A'
keywd_dict['GDSTARID']= 'N/A' # guide star identifier, e.g. 'N/A'
keywd_dict['GS_RA']   = 'N/A' # guide star right ascension, e.g. 'N/A'
keywd_dict['GS_DEC']  = 'N/A' # guide star declination, e.g. 'N/A'
keywd_dict['GSURA']   = 'N/A' # guide star right ascension uncertainty, e.g. 'N/A'
keywd_dict['GSUDEC']  = 'N/A' # guide star declination uncertainty, e.g. 'N/A'
keywd_dict['GS_MAG']  = 'N/A' # guide star magnitude in FGS detector, e.g. 'N/A'
keywd_dict['GSUMAG']  = 'N/A' # guide star magnitude uncertainty, e.g. 'N/A'
keywd_dict['PCS_MODE']= 'N/A' # Pointing Control System mode, e.g. 'N/A'
keywd_dict['GSCENTX'] = 'N/A' # guide star centroid x position, e.g. 'N/A'
keywd_dict['GSCENTY'] = 'N/A' # guide star centroid y position, e.g. 'N/A'
keywd_dict['JITTERMS']= 'N/A' # [arcsec] RMS jitter over the exposure, e.g. 'N/A'

# JWST ephemeris information
keywd_dict['COORDSYS']= 'N/A' # ephemeris coordinate system, e.g. 'N/A'
keywd_dict['EPH_TIME']= 0.0 # [sec] UTC time from ephemeris start time, e.g. 0.0
keywd_dict['JWST_X']  = 0.0 # [km] X spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_Y']  = 0.0 # [km] Y spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_Z']  = 0.0 # [km] Z spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_DX'] = 0.0 # [km/sec] X component of JWST velocity vector, e.g. 0.0
keywd_dict['JWST_DY'] = 0.0 # [km/sec] Y component of JWST velocity vector, e.g. 0.0
keywd_dict['JWST_DZ'] = 0.0 # [km/sec] Z component of JWST velocity vector, e.g. 0.0

# Spacecraft pointing information
keywd_dict['PA_V3']   = 'N/A' # [deg] position angle of V3-axis of JWST, e.g. 'N/A'
keywd_dict['RA_V1']   = 'N/A' # [deg] RA of telescope V1 axis, e.g. 'N/A'
keywd_dict['DEC_V1']  = 'N/A' # [deg] Dec of telescope V1 axis, e.g. 'N/A'

# Aperture pointing information
keywd_dict['APERNAME']= '#TODO' # mnemonic for PDB science aperture used, e.g. #TODO
keywd_dict['PA_APER'] = -999.0 # [deg] position angle of aperture used, e.g. -999.0

# WCS parameters
keywd_dict['WCSAXES'] = 3 # number of World Coordinate System axes, e.g. 3
keywd_dict['CRPIX1']  = 1024.0 # x-coordinate of the reference pixel, e.g. 1024.0
keywd_dict['CRPIX2']  = 128.0 # y-coordinate of the reference pixel, e.g. 128.0
keywd_dict['CRVAL1']  = 5.3196 # RA at the reference pixel (degrees), 5.3196
keywd_dict['CRVAL2']  = -72.98605000000001 # Dec at the reference pixel (degrees), e.g. -72.98605000000001
keywd_dict['CRVAL3']  = 2.5 # Wavelength at the reference pixel (microns), e.g. 2.5
keywd_dict['CTYPE1']  = 'RA---TAN' # first axis coordinate type
keywd_dict['CTYPE2']  = 'DEC--TAN' # second axis coordinate type
keywd_dict['CTYPE3']  = 'WAVE' # third axis coordinate type, e.g. WAVE
keywd_dict['CUNIT1']  = ['deg'] # units for first axis
keywd_dict['CUNIT2']  = ['deg'] # units for seconds axis
keywd_dict['CUNIT3']  = ['micron'] # units for third axis
keywd_dict['CDELT1']  = 0.12 # increment per pixel, axis 1, e.g. 0.12
keywd_dict['CDELT2']  = 0.12 # increment per pixel, axis 2, e.g. 0.12
keywd_dict['CDELT3']  = 0.000672 # increment per pixel, axis 3, e.g. 0.000672
keywd_dict['PC1_1']   = 1.0 # linear transformation matrix element, e.g. 1.0
keywd_dict['PC1_2']   = 0.0 # linear transformation matrix element, e.g. 0.0
keywd_dict['PC2_1']   = 0.0 # linear transformation matrix element, e.g. 0.0
keywd_dict['PC2_2']   = 1.0 # linear transformation matrix element, e.g. 1.0
keywd_dict['PC3_1']   = 1.0 # linear transformation matrix element, e.g. 1.0
keywd_dict['PC3_2']   = 0.0 # linear transformation matrix element, e.g. 0.0
keywd_dict['S_REGION']= 'N/A' # spatial extent of the observation, e.g. 'N/A'
# additional keywords needed in intput to level 2b -> email from Howard Nov 3 2016
keywd_dict['V2_REF']  = 101.1 # location of the aperture reference point in V2 (arcsec): 100-400 arcsec
keywd_dict['V3_REF']  = -202.2 # location of the aperture reference point in V3 (arcsec): -100 to -400 arcsec
keywd_dict['RA_REF']  = 156.11 # RA at the reference point (deg): 0 < RA < 360
keywd_dict['DEC_REF'] = -45.6 # Dec at the reference point (deg): -90 < Dec < +90
keywd_dict['ROLL_REF']= 5.3196 # Roll angle at the reference point (deg), e.g. 5.3196

# Velocity aberration correction
keywd_dict['DVA_RA']  = 0.0 # velocity aberration correction RA offset, e.g. 0.0
keywd_dict['DVA_DEC'] = 0.0 # velocity aberration correction Dec offset, e.g. 0.0
keywd_dict['VA_SCALE']=  1.0 # velocity aberration scale factor, e.g. 1.0

#  Time related keywords
keywd_dict['BARTDELT']= 0.0 # calculated Barycentric time correction, e.g. 0.0
keywd_dict['BSTRTIME']= 0.0 # Solar System Barycentric exposure start time, e.g. 0.0
keywd_dict['BENDTIME']= 0.0 # Solar System Barycentric exposure end time, e.g. 0.0
keywd_dict['BMIDTIME']= 0.0 # Solar System Barycentric exposure mid time, e.g. 0.0
keywd_dict['HELIDELT']= 0.0 # calculated Heliocentric time correction, e.g. 0.0
keywd_dict['HSTRTIME']= 0.0 # Heliocentric exposure start time in MJD, e.g. 0.0
keywd_dict['HENDTIME']= 0.0 # Heliocentric exposure end time in MJD, e.g. 0.0
keywd_dict['HMIDTIME']= 0.0 # Heliocentric exposure mid time in MJD, e.g. 0.0
