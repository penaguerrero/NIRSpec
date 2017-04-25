#!/usr/bin/env python

from __future__ import print_function, division
import os
import re
import numpy as np
from datetime import datetime
from astropy.io import fits
import argparse

# import the scripts I wrote
import functions4fitsfiles as f4t
import hdr_keywords_dictionary as hkwd
import sample_hdr_keywd_vals_dict as shkvd
import sample_hdr_keywd_vals_dict_MOS as shkvdMOS

'''
This script checks that the fits files for the Fixed Slits (FS) have the format that the pipeline
build 7 is expecting.

Example usage:
    The code works from the terminal.
    To create a NEW FS fits file with the updated header type:
        > python check_hdr_keywds.py blah.fits
    for an MOS fits file type:
        > python check_hdr_keywds.py blah.fits -mos

    To simply update the header of the existing FS fits file type:
        > python check_hdr_keywds.py blah.fits -u
    for an MOS fits file type
        > python check_hdr_keywds.py blah.fits -u -mos

'''

### General functions

# create text file to log warnings
def create_warnings_file(fits_file):
    warnings_file_name = fits_file.replace(".fits", "_warnings.txt")
    print ('Warnings file:  ', warnings_file_name)
    tf = open(warnings_file_name, 'w')
    tf.write('### Formating warnings found.\n')
    tf.close()
    return warnings_file_name

# if warnings text file is empty erase it
def check_warnings_file(warnings_file_name):
    with open(warnings_file_name, 'r') as wf:
        lines = wf.readlines()
        if len(lines) == 1:
            os.system('rm '+warnings_file_name)


### Functions to check specific keyword values

def check_value_type(key, val, hkwd_val):
    # Check if type of value matches expected
    valtype = type(val)
    # Check if type of value correspond to what is given, else change it
    if val is not None:
        count = 0
        for v in val:
            if v=='.':
                count += 1
        no_letters_in_string = True
        for char in val:
            if char.isalpha():
                no_letters_in_string = False
        if no_letters_in_string:
            if (count == 0) and (':' not in val) and ('-' not in val):
                if val=='':
                    warning = 'Keyword '+key+' has empty value.'
                    print (warning)
                    return warning
                val = int(val)
            if (count==1) and (':' not in val):
                val = float(val)
            valtype = type(val)
    if (valtype in hkwd_val) or (val in hkwd_val):
        print ('Keyword '+key+' has allowed value type: '+ repr(val))
        warning = None
    else:
        warning = 'Keyword '+key+' has incorrect type of value. Expected: '+str(val)+', got: '+str(valtype)
        print (warning)
    return warning


def check3numbers(key, val):
    warning = 'Keyword '+key+' does not have correct format. Expecting e.g. 0.1.1, received '+val
    r = re.compile('\d.\d.\d') # check that the string has a structure like 0.1.1
    if r.match(val) is None:
        print (warning)
        return warning
    else:
        print ('Format of keyword '+key+' is correct.')

def check_len(key, val, val_len=2):
    string_length = len(val)
    warning = 'Expected '+repr(val_len)+' digits in keyword '+key+'; got '+repr(string_length)+' instead.'
    if string_length == val_len:
        print ('Format of keyword '+key+' is correct.')
    else:
        print (warning)
        return warning

def check_datetimeformat(key, val, check_time, check_date, check_datetime):
    warning = 'Keyword '+key+ ' has incorrect format.'
    if '.' in val:
        vlist = val.split(':')
        v = float(vlist[-1])
        v = str(int(np.round(v, decimals=0)))
        val = vlist[0]+':'+vlist[1]+':'+v
    if check_time:
        val = datetime.strptime(val, '%H:%M:%S')
    if check_date:
        val = datetime.strptime(val, '%Y-%m-%d')
    if check_datetime:
        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
    if isinstance(val, datetime):
        print ('Format of keyword '+key+' is correct.')
    else:
        print (warning)
        return warning


### keyword and format check

class NIRSpec_hdr_format_check:
    def __init__(self, file_keywd_dict, fits_file, only_update, mosdata):
        self.file_keywd_dict = file_keywd_dict
        self.fits_file = fits_file
        self.only_update = only_update
        self.mosdata = mosdata
        self.warnings_list = []
        self.missing_keywds = []

    # Will check keywords against those in hdr_keywords_dictionary.py
    def check_keywds(self, warnings_file_name):
        for hkwd_key, hkwd_val in hkwd.keywd_dict.iteritems():
            # start by making the warning for each keyword None and assigning key and val
            key = hkwd_key

            # Check if keyword is in the file
            if key not in file_keywd_dict:
                self.missing_keywds.append(key)
                warning = '*** Keyword '+key+' is not in header.'
                val = None
            else:
                val = file_keywd_dict[hkwd_key]

                # Check simple standard keyword values
                if val in hkwd_val:
                    print ('Keyword '+key+' has allowed value: '+ val)
                    warning = None
                else:
                    warning = check_value_type(key, val, hkwd_val)

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
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=True,
                                                   check_time=False)
                elif key=='DATE-OBS':
                    warning = check_datetimeformat(key, val, check_date=True, check_datetime=False,
                                                   check_time=False)
                elif key=='TIME-OBS':
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=False,
                                                   check_time=True)

            if warning is not None:
                self.warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')


    def add_keywds(self, only_update):
        '''
        This function adds the missing keywords from the hdr_keywords_dictionary.py (hkwd) file and gives
        the fake values taken from the dictionary sample_hdr_keywd_vals_dict.py (shkvd).
        Args:
            only_update: If false a copy of the original fits file will be created with the
                         updated header.
        '''
        # create name for updated fits file
        updated_fitsfile = fits_file
        if not only_update:
            updated_fitsfile = fits_file.replace('.fits', '_updatedHDR.fits')
            os.system('cp '+fits_file+' '+updated_fitsfile)
        # add missimg keywords
        print ('Path of updated file: ', updated_fitsfile)
        for i, key in enumerate(self.missing_keywds):
            # get the index of the keyword previous to the one you want to add
            prev_key_idx = hkwd.keywd_dict.keys().index(key) - 1
            # add the keyword in the right place from the right dictionary
            new_value = shkvd.keywd_dict[key]
            after_key = shkvd.keywd_dict.keys()[prev_key_idx]
            if self.mosdata:
                new_value = shkvdMOS.keywd_dict[key]
                after_key = shkvdMOS.keywd_dict.keys()[prev_key_idx]
            fits.setval(updated_fitsfile, key, value=new_value, after=after_key)
        print ('\n New header: ')
        hdulist = fits.open(updated_fitsfile)
        hdr = hdulist[0].header
        print (repr(hdr))

    def perform_check(self):
        # create text file to log warnings
        print('')
        warnings_file_name = create_warnings_file(fits_file)

        # check the keywords
        print('\n   Starting keyword check...')
        self.check_keywds(warnings_file_name)

        # if warnings text file is empty erase it
        check_warnings_file(warnings_file_name)

        # if the warnings are just the misssing keywords and an empty value on VISIT_ID erase it
        expected_warnings = ['Keyword VISIT_ID has empty value.', '*** Keyword V2_REF is not in header.',
                             '*** Keyword V3_REF is not in header.', '*** Keyword RA_REF is not in header.',
                             '*** Keyword DEC_REF is not in header.', '*** Keyword ROLL_REF is not in header.']
        if self.warnings_list == expected_warnings:
            with open(warnings_file_name, 'r') as wf:
                lines = wf.readlines()
                os.system('rm '+warnings_file_name)

        # create new file with updated header
        print('\n   Adding keywords...')
        self.add_keywds(self.only_update)



if __name__ == '__main__':

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("-u",
                        dest="only_update",
                        action='store_true',
                        default=False,
                        help='Use if NOT wanting to create a new file with updated header.')
    parser.add_argument("-mos",
                        dest="mosdata",
                        action='store_true',
                        default=False,
                        help='Use this for processing MOS data.')
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    only_update = args.only_update
    mosdata = args.mosdata

    # Define paths - THIS IS FOR TESTING THE CODE ONLY
    #pileline_path = os.path.abspath(os.curdir).split('pipeline')[0]
    #build_path = os.path.join(pileline_path, 'pipeline/build7')
    #fits_file = os.path.join(pileline_path, "pipeline/src/jwtest1001001_01101_00001_NRS1_uncal_mod.fits")
    #only_update = False

    # read the keywords and corresponding values from fits file directly
    file_keywd_dict = f4t.read_hdrfits(fits_file)

    # Perform the keyword check
    t = NIRSpec_hdr_format_check(file_keywd_dict, fits_file, only_update, mosdata)
    t.perform_check()

    print ('\n * Script  check_hdr_keywds.py  finished * \n')


