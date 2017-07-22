from __future__ import print_function, division
import numpy as np
import os
from glob import glob

from .. import core_utils
from ..D_extract_2d import extract_2d_utils
from ..auxiliary_code import CV3_testdata_used4build7


"""
This script compares pipeline WCS info with ESA results for FIXED SLIT.

slit = 'S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1'
det = 'NRS1' or 'NRS2'
subarray_origin = [SUBSTRT1, SUBSTRT2] from original image; needed since SUBSTRT in
                    world_coordinates file is not in full frame reference
"""

def do_idl_match(arrA, arrB):
    """
    This function does the same that the IDL match function does. It finds the elements common
    in both arrays and it returns two arrays with the index of those elements in each array.
    (The arrays do not need to be the same length)
    Args:
        arrA: numpy array
        arrB: numpy array

    Returns:
        subA: numpy array of index of arrA from elements also present in arrB
        subB: numpy array of index of arrB from elements also present in arrA
    """
    # First turn the arrays into sets
    setA, setB = set(arrA), set(arrB)
    # Find the intersection with Python built-in method
    interAB = setA.intersection(setB)
    # Find the index corresponding to the intersection elements and return them as arrays
    subA, subB = [], []
    for i, ai in enumerate(arrA):
        if ai in interAB:
            subA.append(i)
    for i, bi in enumerate(arrB):
        if bi in interAB:
            subB.append(i)
    return np.array(subA), np.array(subB)


def compare_wcs(infile_hdr, infile_name, data_extension, esaroot):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

    Args:
        infile_hdr: list, header of the output fits file from the 2d_extract step
        infile_name: str, name of the output fits file from the 2d_extract step (with full path)
        data_extension: int, number of the data extension from the infile
        esaroot: str, full path of where to find all ESA intermediary products to make comparisons for the tests

    Returns:

    """

    # Run compute_world_coordinates.py in order to produce the necessary file
    os.system("compute_world_coordinates('"+infile_name+"')")

    # Read the resulting file
    cwc_fname = infile_name.replace(".fits", "_world_coordinates.fits")
    fdata = core_utils.get_filedata(cwc_fname, ext=data_extension)
    pwave = fdata[0,:,:]
    pdy = fdata[3,:,:]
    pskyx = fdata[1,:,:]
    pskyy = fdata[2,:,:]
    # get the origin of the subwindow and the grating from the extract_2d file header
    det = extract_2d_utils.find_DETECTOR(infile_hdr)
    # the science extensions are always the same in the extract_2d file
    sci_data_list, px0_list, py0_list, sltname_list = [], [], [], []
    grat = core_utils.get_keywd_val(infile_name, "GRATING", ext=0)
    filt = core_utils.get_keywd_val(infile_name, "FILTER", ext=0)
    print ("Grating = ", grat, "   Filter = ", filt)
    sci_exts = core_utils.get_sci_extensions(infile_name)
    for se in sci_exts:
        sci_data_list.append(core_utils.get_filedata(infile_name, ext=se))
        px0_se = core_utils.get_keywd_val(infile_name, "SLTSTRT1", ext=se)+core_utils.get_keywd_val(infile_name, "SUBSTRT1", ext=0)-1
        py0_se = core_utils.get_keywd_val(infile_name, "SLTSTRT2", ext=se)+core_utils.get_keywd_val(infile_name, "SUBSTRT2", ext=0)-1
        px0_list.append(px0_se)
        py0_list.append(py0_se)
        sltname = core_utils.get_keywd_val(infile_name, "SLTNAME", ext=se)
        sltname_list.append(sltname)
    n_p = np.shape(pwave)
    npx = n_p[0]
    npy = n_p[1]
    px = np.arange(1, npx+1)+np.array(px0_list)
    py = np.arange(1, npy+1)+np.array(py0_list)
    print  ("Pipeline subwindow corner pixel ID: ", px0_list, py0_list)

    # get the corresponding ESA file to the input file
    if det == "NRS1":
        file4detector = 0
    elif det == "NRS2":
        file4detector = 1
    NID = CV3_testdata_used4build7.CV3_testdata_dict["FS"][grat]["NID"]
    CV3filename = CV3_testdata_used4build7.CV3_testdata_dict["FS"][grat]["CV3filename"][file4detector]
    if len(NID) != 1:
        # Find the NID of the file used
        initial_infile_name = infile_name.replace(glob("_rate*", ""))
        for i, fi in enumerate(CV3_testdata_used4build7.CV3_testdata_dict["FS"][grat]["level1Bfilenames"]):
            for j, f in enumerate(fi):
                if initial_infile_name in f:
                    NID = CV3_testdata_used4build7.CV3_testdata_dict["FS"][grat]["NID"][i]
                    CV3filename = CV3_testdata_used4build7.CV3_testdata_dict["FS"][grat]["CV3filename"][i][j]
    Vnumber = CV3filename.split("_")[0].replace("NRS", "")
    ESA_dir_name = Vnumber+"_"+NID+"_JLAB88"

    esafile_dir = esaroot+"RegressionTestData_CV3_March2017_FixedSlit/"+ESA_dir_name+"/"+ESA_dir_name+"_trace_SLIT"

    # read in ESA data
    for sltname in sltname_list:
        # change the format of the string to match the ESA trace
        sltname = sltname.split("S")[1]
        if sltname[-1] == "A1":
            sltname = "A_"+sltname.split("A")[0]+"_1"
        elif sltname[-1] == "A2":
            sltname = "A_"+sltname.split("A")[0]+"_2"
        elif sltname[-1] == "A":
            sltname = "A_"+sltname.split("A")[0]+"_"
        elif sltname[-1] == "B":
            sltname = "B_"+sltname.split("A")[0]+"_"
        esa_file = esafile_dir+"/Trace_SLIT_"+sltname+ESA_dir_name+".fits"
    if det == "NRS1":
        eflux = core_utils.get_filedata(esa_file, ext=1)
        ewave = core_utils.get_filedata(esa_file, ext=4)
        edy = core_utils.get_filedata(esa_file, ext=5)
    else:
        eflux = core_utils.get_filedata(esa_file, ext=6)
        ewave = core_utils.get_filedata(esa_file, ext=9)
        edy = core_utils.get_filedata(esa_file, ext=10)
    n_p = np.shape(eflux)
    nex = n_p[0]
    ney = n_p[1]

    # get the origin of the subwindow
    CRVAL1 = core_utils.get_keywd_val(esa_file, "CRVAL1", ext=0)
    CRPIX1 = core_utils.get_keywd_val(esa_file, "CRPIX1", ext=0)
    CRVAL2 = core_utils.get_keywd_val(esa_file, "CRVAL1", ext=0)
    CRPIX2 = core_utils.get_keywd_val(esa_file, "CRPIX1", ext=0)
    print ("CRVAL1 = ", CRVAL1, "    CRPIX1 = ", CRPIX1)
    print ("CRVAL2 = ", CRVAL1, "    CRPIX2 = ", CRPIX1)
    if det == "NRS1":
        ex0 = CRVAL1 - CRPIX1 + 1
        ey0 = CRVAL2 - CRPIX2 + 1
    else:
        ex0 = 2048 - CRPIX1 + CRVAL1
        ey0 = 2048 - CRVAL2 + CRPIX2
    print ("ESA subwindow corner pixel ID: ", ex0, ey0)
    ex = np.arange(nex) + ex0
    ey = np.arange(ney) + ey0

    # match up the correct elements in each data set
    subpx, subex = do_idl_match(px, ex)
    subpy, subey = do_idl_match(py, ey)
    countx, county = len(subpx), len(subpy)
    print ("Matched elements in the 2D spectra: ", countx, county)
    # flatten the arrays
    imp = subpy.flatten()
    ime = subey.flatten()
    print (ime[0:9])

    # get the difference between the two in units of resels
    # don't include pixels where one or the other solution is 0 or NaN
    ig, igy = np.array([]), np.array([])
    for i, imp_i in enumerate(imp):
        if pwave[imp_i] != 0.0:
            if ewave[ime[i]] != 0.0:
                if np.isfinite(pwave[imp_i]):
                    if np.isfinite(ewave[ime[i]]):
                        ig = np.append(ig, i)
        # and now in the y direction
        if pdy[imp_i] != 0.0:
            if edy[ime[i]] != 0.0:
                if np.isfinite(pdy[imp_i]):
                    if np.isfinite(edy[ime[i]]):
                        igy = np.append(igy, i)
    delwav1, delwav2 = [], []
    for imp_i in imp:
        for i in ig:
            d = pwave[imp_i[i]]
            delwav1.append(d)
    for ime_i in ime:
        for i in ig:
            d = ewave[ime_i[i]]
            delwav2.append(d)
    delwave = np.array(delwav1) - 6 - np.array(delwav2)



