from __future__ import print_function, division
import numpy as np
import os
import subprocess
import sys
from astropy.io import fits


"""
This script compares pipeline WCS info with ESA results for Multi-Object Spectroscopy (MOS) data.

slit = [shutter, ID]  --> e.g. [3, 19924]
det = 'NRS1' or 'NRS2'

"""

def get_sci_extensions(fits_file_name):
    """
    This functions obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    sci_list = []
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_list.append(ext)
    return sci_list

def compare_wcs(infile_hdr, infile_name, data_extension, esaroot, auxiliary_code_path=None):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

    Args:
        infile_hdr: list, header of the output fits file from the 2d_extract step
        infile_name: str, name of the output fits file from the 2d_extract step (with full path)
        data_extension: int, number of the data extension from the infile
        esaroot: str, full path of where to find all ESA intermediary products to make comparisons for the tests
        auxiliary_code_path: str, path where to find the auxiliary code. If not set the code will assume
                            it is in the the auxiliary code directory

    Returns:

    """

    # Run compute_world_coordinates.py in order to produce the necessary file
    # !!! note that the code expects to be in the build environment !!!
    if auxiliary_code_path is None:
        auxiliary_code_path = "./"
    #print ("auxiliary_code_path=", auxiliary_code_path)
    compute_world_coords_code = os.path.join(auxiliary_code_path, "compute_world_coordinates.py")
    run_compute_world_coords = subprocess.Popen(compute_world_coords_code+" mos "+infile_name, shell=True,
                                                stdout=subprocess.PIPE)
    run_compute_world_coords.wait()

    # The world coordinate file was created but it needs to be renamed
    basenameinfile_name = os.path.basename(infile_name)
    fileID = basenameinfile_name.split("_")[0]   # obtain the id of the file
    working_dir_path = os.getcwd()
    wcoordfile = working_dir_path+"/"+fileID+"_world_coordinates.fits"
    # to move file to location of infile
    #cwc_fname = infile_name.replace(".fits", "_world_coordinates.fits")
    # to rename file within the working directory
    cwc_fname = basenameinfile_name.replace(".fits", "_world_coordinates.fits")
    os.system("mv "+wcoordfile+" "+cwc_fname)
    #print (cwc_fname)

    # Read the resulting file and get the desired arrays
    fdata = fits.getdata(cwc_fname, data_extension)
    pwave = fdata[0,:,:]
    pdy = fdata[3,:,:]
    pskyx = fdata[1,:,:]
    pskyy = fdata[2,:,:]

    # reverse the order of the rows
    pwave = pwave[::-1]
    pdy = pdy[::-1]

    # get the origin of the subwindow and the grating from the extract_2d file header
    det = infile_hdr["DETECTOR"]

    # the science extensions are always the same in the extract_2d file
    sci_data_list, px0_list, py0_list = [], [], []
    grat = fits.getval(infile_name, "GRATING", 0)
    filt = fits.getval(infile_name, "FILTER", 0)
    print ("Grating:", grat, "    Filter:", filt)
    sci_exts = get_sci_extensions(infile_name)
    for se in sci_exts:
        sci_data_list.append(fits.getdata(infile_name, se))
        px0_se = fits.getval(infile_name, "SLTSTRT1", se)+fits.getval(infile_name, "SUBSTRT1", 0)-1
        py0_se = fits.getval(infile_name, "SLTSTRT2", se)+fits.getval(infile_name, "SUBSTRT2", 0)-1
        px0_list.append(px0_se)
        py0_list.append(py0_se)
    n_p = np.shape(pwave)
    #print ("n_p=", n_p)
    npx = n_p[0]
    npy = n_p[1]
    px = np.arange(npx)+np.array(px0_list)
    py = np.arange(npy)+np.array(py0_list)
    print  ("Pipeline subwindow corner pixel ID: ", px0_list, py0_list)
    #print  ("px, py: ", px, py)

    # get the corresponding ESA file to the input file
    # to do this, the script needs the python dictionary of the CV3 data
    sys.path.append(auxiliary_code_path)
    import CV3_testdata_used4build7
    if det == "NRS1":
        file4detector = 0
    elif det == "NRS2":
        file4detector = 1
    for NID, nid_dict_key in CV3_testdata_used4build7.CV3_testdata_dict["MOS"]["NID"].items():
        if nid_dict_key["grism"] == grat:
            if nid_dict_key["filter"] == filt:
                CV3filename = nid_dict_key["CV3filename"][file4detector]
                print ("NID of ESA file:", NID)
                print("CV3filename =", CV3filename)
    # read in ESA data
    # the ESA direcoty names use/follow their name conventions
    ESA_dir_name = CV3filename.split("_")[0].replace("NRS", "")+"_"+NID+"_JLAB88"
    esafile_directory = esaroot+"/RegressionTestData_CV3_March2017_MOS/"+ESA_dir_name+"/"+ESA_dir_name+"_trace_MOS"

    #esafile = os.path.join(esafile_directory, )
    # This fixed path is just to test that the code works
    #esafile = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/build7/test_data/MOS_CV3/complete_pipeline_testset/V9621500100101_short_msa.fits"
    #print ("Using this ESA file: \n", esafile)
    #if det == "NRS1":
    #    eflux =



if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    auxiliary_code_path = pipeline_path+"/src/pytests/calwebb_spec2_pytests/auxiliary_code"
    infile_name = pipeline_path+"/build7/test_data/MOS_CV3/complete_pipeline_testset/jwtest1010001_01101_00001_NRS1_uncal_rate_short_assign_wcs_extract_2d.fits"
    data_extension = 1
    infile_hdr = fits.getheader(infile_name, 0)
    esaroot = pipeline_path+"/build7/test_data/ESA_intermediary_products"

    # Run the principal function of the script
    compare_wcs(infile_hdr, infile_name, data_extension, esaroot, auxiliary_code_path=auxiliary_code_path)