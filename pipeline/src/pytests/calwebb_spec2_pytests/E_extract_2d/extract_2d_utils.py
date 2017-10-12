"""
This file contains the functions which will be used to test the extract_2d step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


### VERIFICATION FUNCTIONS

def s_ext2d_exists(output_hdul):
    """
    This function checks that the keyword S_EXTR2D was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_EXTR2D" in output_hdul
    return result


### VALIDATION FUNCTIONS

# Functions to compare WCS pipeline info with ESA results
# these functions are a translation into Python from James' IDL
# routines compare_wcs_fs.pro and compare_wcs_mos.pro

def check_FS_true(output_hdul):
    """
    This function checks if the fits file is a Fixed Slit.
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be Fixed Slit
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_FIXEDSLIT":
            result = True
    return result

def find_which_slit(output_hdul):
    """
    This function determines which Fixed Slit was used
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        s: string, slit used or is None if not in the list
    """
    # the order of this list corresponds to the
    slits = ["S200A1", "S200A2", "S200B1", "S400A1", "S1600A1"]
    if "FXD_SLIT" in output_hdul:
        for i, s in enumerate(slits):
            if "FXD_SLIT" == s:
                return i+1, s


def check_MOS_true(output_hdul):
    """
    This function checks if the fits file is Multi-Object Spectroscopy (MOS).
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be MOS data
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_MSASPEC":
            result = True
    return result


def check_IFU_true(output_hdul):
    """
    This function checks if the fits file is IFU data.
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be IFU data
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_IFU":
            result = True
    return result


def find_DETECTOR(output_hdul):
    """
    This function determines which detector was used
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        det: string, either NRS1 or NRS2
    """
    if "DETECTOR" in output_hdul:
        det = output_hdul["DETECTOR"]
        return det

