from __future__ import print_function, division
import os
import string
import re
from datetime import datetime

# import the scripts I wrote
import tests.functions4tests as f4t
import tests.hdr_keywords_dictionary as hkwd

'''
This script checks that the fits files for the Fixed Slits (FS) have the format that the pipeline
build 7 is expecting.
'''

### General functions

# create text file to log warnings
def create_warnings_file(fits_file):
    warnings_file_name = fits_file.replace(".fits", "_warnings.txt")
    print ('Warnings file:  ', warnings_file_name)
    tf = open(warnings_file_name, 'w')
    tf.write('Formating warnings found.\n')
    tf.close()
    return warnings_file_name

# if warnings text file is empty erase it
def check_warnings_file(warnings_file_name):
    with open(warnings_file_name, 'r') as wf:
        lines = wf.readlines()
        if len(lines) == 1:
            os.system('rm '+warnings_file_name)


### Assert functions to check specific keyword values

def format3number_assert(val):
    r = re.compile('\d.\d.\d') # check that the string has a structure like 0.1.1
    assert r.match(val) # if r.match(val) is None, the expected val format did match the received one
def check3numbers(key, val):
    warning = 'Keyword ', key, ' does not have correct format. Expecting e.g. 0.1.1, received ', val
    a = format3number_assert(val)
    if a is not None:
        print (warning)
        return warning
    else:
        print ('Format of keyword '+key+' is correct.')

def assert_len2(string_length, val_len):
    assert string_length==val_len
def check_len(key, val, val_len=2):
    string_length = len(val)
    a = assert_len2(string_length, val_len)
    warning = 'Expected ', val_len,' digits in keyword ', key, '; got ',string_length, ' instead.'
    if a is None:
        print (warning)
        return warning
    else:
        print ('Format of keyword '+key+' is correct.')

def test_assertdatetimeformat(val, check_time, check_date, check_datetime):
    if '.' in val:
        val = val.split('.')[0]
    if check_time:
        val = datetime.strptime(val, '%H:%M:%S')
    if check_date:
        val = datetime.strptime(val, '%Y-%m-%d')
    if check_datetime:
        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
    assert isinstance(val, datetime)

def test_checkdatetimeformat(val, check_time=True, check_date=False, check_datetime=False):
    warning = 'Format of date/time is incorrect.'
    a = test_assertdatetimeformat(val, check_time, check_date, check_datetime)
    if a is None:
        print (warning)
        return warning


### Unittest of keyword and format check

class NIRSpec:
    def __init__(self, file_keywd_dict, fits_file):
        self.file_keywd_dict = file_keywd_dict
        self.fits_file = fits_file

    # Will check keywords against those in hdr_keywords_dictionary.py
    def check_keywds(self, warnings_file_name):
        for hkwd_key, hkwd_val in hkwd.keywd_dict.iteritems():
            # start by making the warning for each keyword None and assigning key and val
            warning = None
            key = hkwd_key
            # Check if keyword is in the file
            if key not in file_keywd_dict:
            	warning = '*** Keyword ', key, ' is not in header.'
	    else:
            	val = file_keywd_dict[hkwd_key]
                
		# Check simple standard keyword values
                if val in hkwd_val:
                    print ('yay! the keyword value is allowed!')
                else:
                    # Check if type of value matches expected
                    valtype = type(val)
                    if valtype in hkwd_val:
                        print ('yay! Type of value is allowed!')
                    else:
                        print ('Type of value is not allowed. Expected: ', val, ' \n Instead got: ', valtype)
                        warning = 'Value of keyword ', key, ' has incorrect format.'
                    # Check for specific keywords
                    if key=='DPSW_VER':
                        warning = check3numbers(key, val)
                    elif (key=='VISITGRP') or (key=='ACT_ID'):
                        warning = check_len(key, val, val_len=2)
                    elif (key=='OBSERVTN') or (key=='VISIT'):
                        warning = check_len(key, val, val_len=3)
                    elif (key=='EXPOSURE'):
                        warning = check_len(key, val, val_len=5)
                    elif (key=='DATE') or (key=='VSTSTART'):
                        warning = test_checkdatetimeformat(val, check_date=False, check_datetime=True, check_time=False)
                    elif (key=='DATE-OBS'):
                        warning = test_checkdatetimeformat(val, check_date=True, check_datetime=False, check_time=False)
                    elif (key=='TIME-OBS'):
                        warning = test_checkdatetimeformat(val, check_date=False, check_datetime=False, check_time=True)

                    # Failed all checks, write warning
                    if warning is None:
                        warning = 'Value of keyword '+key+' is not an allowed value.'

                if warning is not None:
                    with open(warnings_file_name, "a") as tf:
                        tf.write(warning+'\n')

    def perform_test(self):
        # create text file to log warnings
        warnings_file_name = create_warnings_file(fits_file)

        # check the keywords
        self.check_keywds(warnings_file_name)

        # if warnings text file is empty erase it
        check_warnings_file(warnings_file_name)



if __name__ == '__main__':
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
        file_keywd_dict = f4t.read_hdrtxt(FS_hdr_sample)
    else:
        # read header from fits file directly
        file_keywd_dict = f4t.read_hdrfits(fits_file)

    # Perform the keyword check
    t = NIRSpec(file_keywd_dict, fits_file)
    t.perform_test()


