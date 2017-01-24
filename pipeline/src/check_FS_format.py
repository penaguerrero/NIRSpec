from __future__ import print_function, division

import os
import re
import string
from datetime import datetime

# import the scripts I wrote
import functions4fitsfiles as f4t
import hdr_keywords_dictionary as hkwd

'''
This script checks that the fits files for the Fixed Slits (FS) have the format that the pipeline
build 7 is expecting.
'''

# Define script running parameters
read_header_from_fits = False   # If False, the header will be read from a text file


# Define paths
pileline_path = string.split(os.getcwd(), sep='pipeline')[0]
build_path = os.path.join(pileline_path, 'pipeline/build7')
header_txt_path = os.path.join(pileline_path, 'pipeline/src/tests/txt_sample_header_keywds')
FS_hdr_sample = os.path.join(header_txt_path, 'FS_sample-jwtest1001001_01101_00001_NRS1_uncal_mod.txt')
fits_file = os.path.join(pileline_path, "pipeline/src/jwtest1001001_01101_00001_NRS1_uncal_mod.fits")

# read the keywords and corresponding values
if not read_header_from_fits:
    # read from text file
    keywd_dict = f4t.read_hdrtxt(FS_hdr_sample)
else:
    # read header from fits file directly
    keywd_dict = f4t.read_hdrfits(fits_file)

# create text file to log warnings
warnings_file_name = fits_file.replace(".fits", "_warnings.txt")
print ('Warnings file:  ', warnings_file_name)
tf = open(warnings_file_name, 'w')
tf.write('Formating warnings found.\n')
tf.close()

# check keyword values
for key, val in keywd_dict.iteritems():
    if key in hkwd.keywd_dict:
        print ('keywd=', key, '  val=',val, '  type=', type(val))
        # Check for specific keywords
        if key=='DPSW_VER':
            r = re.compile('\d.\d.\d')
            if r.match(val) is not None:
                print('version of ', key, ' has correct format.')
            else:
                warning = '*** Version of keyword ', key, ' does not have correct format.'
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        elif (key=='VISITGRP') or (key=='ACT_ID'):
            string_length = len(val)
            if string_length == 2:
                print('Keyword has correct number of digits.')
            else:
                warning = '*** Expected 2 digits in keyword ', key, '; got ',string_length, ' instead.'
                print(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        elif (key=='OBSERVTN') or (key=='VISIT'):
            string_length = len(val)
            if string_length == 3:
                print('Keyword has correct number of digits.')
            else:
                warning = '*** Expected 3 digits in keyword ', key, '; got ',string_length, ' instead.'
                print(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        elif (key=='EXPOSURE'):
            string_length = len(val)
            if string_length == 5:
                print('Keyword has correct number of digits.')
            else:
                warning = '*** Expected 5 digits in keyword ', key, '; got ',string_length, ' instead.'
                print(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        elif val in hkwd.keywd_dict[key]:
            print('value in list, yay!')
        else:
            print ('hkwd.keywd_dict[key] = ', hkwd.keywd_dict[key])
            #  if the number does not have a point then it is an integer
            if key!='DPSW_VER':
                r = re.compile('\d')
                if r.match(val) is not None:
                    if ('.' not in val) and ('fits' not in val) and (':' not in val) and ('-' not in val):
                        val = int(val)
                    elif ('.' in val) and ('fits' not in val) and (':' not in val) and ('-' not in val):
                        val = float(val)
            # check if the date has the right format
            elif ('-' in val) and ('T' in val) and (':' in val):
                try:
                    val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
                except:
                    ValueError
                    warning = '*** Date of keyword ', key, ' does not have the correct format.'
                    print(warning)
                    with open(warnings_file_name, "a") as tf:
                        tf.write(warning+'\n')
            elif ('-' in val) and (':' in val) and ('T' not in val):
                try:
                    val = datetime.strptime(val, '%H:%M:%S')
                except:
                    ValueError
                    warning = '*** Time of keyword ', key, ' does not have the correct format.'
                    print(warning)
                    with open(warnings_file_name, "a") as tf:
                        tf.write(warning+'\n')
                print (val)
            for type_in_list in hkwd.keywd_dict[key]:
                print(type(val), 'type_in_list=', type_in_list)
                if type(val) is type_in_list:
                    print ('type of value in list!')

# if warnings text file is empty erase it
with open(warnings_file_name, 'r') as wf:
    lines = wf.readlines()
    if len(lines) == 1:
        os.system('rm '+warnings_file_name)


print ('Script  test_FS.py  finished.')
