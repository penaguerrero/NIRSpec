#!/usr/bin/env python

from __future__ import print_function, division

import os
import re
from datetime import datetime

import numpy as np
import pytest

# import the scripts I wrote
import functions4fitsfiles as f4t
import hdr_keywords_dictionary as hkwd

'''
This script checks that the fits files for the Fixed Slits (FS) have the format that the pipeline
build 7 is expecting.

To test and get an html report run code with:
    > pytest --html=report.html
'''


# Define paths
pileline_path = os.path.abspath(os.curdir).split('pipeline')[0]
fits_file = os.path.join(pileline_path, "pipeline/src/jwtest1001001_01101_00001_NRS1_uncal_mod.fits")



# Create a fixture that reads the file and obtains the header
@pytest.fixture(scope='module')   # the scope means this function info will be shared within a module
def get_hdr():

    # read header from fits file directly
    file_keywd_dict = f4t.read_hdrfits(fits_file)

    dictionary = {'file_keywd_dict':file_keywd_dict}

    return dictionary



# Check for the existence of keywords and verify correct format of corresponding value
def test_NIRSpec_hdr_keywds(get_hdr):

    # read hdr
    file_keywd_dict = get_hdr['file_keywd_dict']

    # check keywords against those in hdr_keywords_dictionary.py
    for hkwd_key, hkwd_val in hkwd.keywd_dict.iteritems():

        # start by making the warning for each keyword None and assigning key and val
        key = hkwd_key

        # Check if keyword is in the file
        assert [k==key for k in file_keywd_dict.keys()]

        if key in file_keywd_dict:
            val = file_keywd_dict[hkwd_key]
        else:
            val=None

        # Check simple standard keyword values
        assert [v==val for v in hkwd_val]

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
        assert [v==valtype for v in hkwd_val]

        # Check for specific keywords
        if key=='DPSW_VER':
            r = re.compile('\d.\d.\d')
            assert r.match(val)
        elif (key=='VISITGRP') or (key=='ACT_ID'):
            string_length = len(val)
            assert string_length==2
        elif (key=='OBSERVTN') or (key=='VISIT'):
            string_length = len(val)
            assert string_length==3
        elif (key=='EXPOSURE'):
            string_length = len(val)
            assert string_length==5
        elif (key=='DATE') or (key=='VSTSTART'):
            val = str(datetime.strptime(val, '%Y-%m-%dT%H:%M:%S'))
            assert isinstance(val, str)
        elif (key=='DATE-OBS'):
            val = datetime.strptime(val, '%Y-%m-%d')
            assert isinstance(val, datetime)
        elif (key=='TIME-OBS'):
            vlist = val.split(':')
            v = float(vlist[-1])
            v = str(int(np.round(v, decimals=0)))
            val = vlist[0]+':'+vlist[1]+':'+v
            val = datetime.strptime(val, '%H:%M:%S')
            assert isinstance(val, datetime)



