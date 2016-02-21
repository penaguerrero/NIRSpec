from __future__ import print_function, division
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import os


'''
This script crates the plots in DeltaV2-DeltaV3 space, that compare the 3 tests ran.
     TEST1 - Average positions P1 and P2, transform to V2-V3 space, and compare to average
             reference positions (V2-V3 space)
     TEST2 - Transform individual positions P1 and P2 to V2-V3 space, average V2-V3 space
             positions, and compare to average reference positions.
     TEST3 - Transform P1 and P2 individually to V2-V3 space and compare star by star and
             position by position.
'''


def v2v3mean_plot(case, meanV2, meanV3, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plot in V2-V3 space of the 3 tests: averaging in pixel space, averaging on sky,
     and no averaging.
    Args:
        case               -- string, for example '491Scene1_rapid_real_bgFrac0.3'
        meanV2             -- list of 3 numpy array of mean values of V2 for Tests 1, 2, and 3
        meanV3             -- list of 3 numpy array of mean values of V3 for Tests 1, 2, and 3
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # check if the plot is for an Nk set
    basename = case
    if not isinstance(meanV2, float):
        basename = case+'_'+str(len(meanV2[0]))+'samples'

    # Make the plot of MEANS
    plot_title = 'Mean Residual Values'
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.suptitle(plot_title, fontsize=18, y=0.96)
    plt.title(basename)
    plt.xlabel(r'$\Delta$V2')
    plt.ylabel(r'$\Delta$V3')
    xmin, xmax = -0.05, 0.05
    ymin, ymax = -0.05, 0.05
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, xmin, xmax*2, colors='k', linestyles='dashed')
    plt.vlines(0.0, ymin, ymax*2, colors='k', linestyles='dashed')
    plt.plot(meanV2[0], meanV3[0], 'b^', ms=10, alpha=0.7, label='Avg in Pixel Space')
    plt.plot(meanV2[1], meanV3[1], 'go', ms=10, alpha=0.7, label='Avg in Sky')
    plt.plot(meanV2[2], meanV3[2], 'r*', ms=13, alpha=0.7, label='No Avg')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if isinstance(meanV2, float):
        textinfig = r'V2$\mu$ = %0.2f    V3$\mu$ = %0.2f' % (meanV2, meanV3)
        ax1.annotate(textinfig, xy=(1.02, 0.35), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'means_'+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'means_'+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')


def v2theta_plot(case, meanV2, theta, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plot in V2-theta space of the 3 tests: averaging in pixel space, averaging on sky,
     and no averaging.
    Args:
        case               -- string, for example '491Scene1_rapid_real_bgFrac0.3'
        meanV2             -- list of 3 numpy array of theta values of V2 for Tests 1, 2, and 3
        theta              -- list of 3 numpy array of theta values for Tests 1, 2, and 3
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # check if the plot is for an Nk set
    basename = case
    if not isinstance(meanV2, float):
        basename = case+'_'+str(len(meanV2[0]))+'samples'

    # Make the plot of V2-THETA
    plot_title = r'Residual Mean Calculated Angle, $\theta$'
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.suptitle(plot_title, fontsize=18, y=0.96)
    plt.title(basename)
    plt.xlabel(r'$\Delta$V2')
    plt.ylabel(r'$\theta$')
    xmin, xmax = -0.01, 0.01
    ymin, ymax = -40.0, 40.0
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, xmin, xmax*2, colors='k', linestyles='dashed')
    plt.vlines(0.0, ymin, ymax*2, colors='k', linestyles='dashed')
    plt.plot(meanV2[0], theta[0], 'b^', ms=10, alpha=0.7, label='Avg in Pixel Space')
    plt.plot(meanV2[1], theta[1], 'go', ms=10, alpha=0.7, label='Avg in Sky')
    plt.plot(meanV2[2], theta[2], 'r*', ms=13, alpha=0.7, label='No Avg')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if isinstance(meanV2, float):
        textinfig = r'V2$\mean$ = %0.2f    $\theta$ = %0.2f' % (meanV2, theta)
        ax1.annotate(textinfig, xy=(1.02, 0.35), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'thetaV2_'+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'thetaV2_'+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')


def v3theta_plot(case, meanV3, theta, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plot in V3-theta space of the 3 tests: averaging in pixel space, averaging on sky,
     and no averaging.
    Args:
        case               -- string, for example '491Scene1_rapid_real_bgFrac0.3'
        meanV3             -- list of 3 numpy array of theta values of V3 for Tests 1, 2, and 3
        theta              -- list of 3 numpy array of theta for Tests 1, 2, and 3
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # check if the plot is for an Nk set
    basename = case
    if not isinstance(meanV3, float):
        basename = case+'_'+str(len(meanV3[0]))+'samples'

    # Make the plot of V3-THETA
    plot_title = r'Residual Mean Calculated Angle, $\theta$'
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.suptitle(plot_title, fontsize=18, y=0.96)
    plt.title(basename)
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\Delta$V3')
    xmin, xmax = -40.0, 40.0
    ymin, ymax = -0.02, 0.02
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, xmin, xmax*2, colors='k', linestyles='dashed')
    plt.vlines(0.0, ymin, ymax*2, colors='k', linestyles='dashed')
    plt.plot(theta[0], meanV3[0], 'b^', ms=10, alpha=0.7, label='Avg in Pixel Space')
    plt.plot(theta[1], meanV3[1], 'go', ms=10, alpha=0.7, label='Avg in Sky')
    plt.plot(theta[2], meanV3[2], 'r*', ms=13, alpha=0.7, label='No Avg')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if isinstance(meanV2, float):
        textinfig = r'V3$\mean$ = %0.2f    $\theta$ = %0.2f' % (meanV3, theta)
        ax1.annotate(textinfig, xy=(1.02, 0.35), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'thetaV3_'+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'thetaV3_'+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')


def theta_plot(case, theta, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plot of theta for the 3 tests: averaging in pixel space, averaging on sky,
     and no averaging.
    Args:
        case               -- string, for example '491Scene1_rapid_real_bgFrac0.3'
        theta              -- list of 3 numpy array of theta for Tests 1, 2, and 3
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # check if the plot is for an Nk set
    basename = case
    if not isinstance(meanV3, float):
        basename = case+'_'+str(len(meanV3[0]))+'samples'

    # Make the plot of THETA
    plot_title = r'Residual Mean Calculated Angle, $\theta$'
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.suptitle(plot_title, fontsize=18, y=0.96)
    plt.title(basename)
    plt.xlabel('Sample Number')
    plt.ylabel(r'$\theta$')
    xmin, xmax = -100.0, 5100.0
    ymin, ymax = -40.0, 40.0
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, -1000, 6000, colors='k', linestyles='dashed')
    plt.vlines(0.0, -1000, 6000, colors='k', linestyles='dashed')
    plt.plot(theta[0], 'b^', ms=10, alpha=0.7, label='Avg in Pixel Space')
    plt.plot(theta[1], 'go', ms=10, alpha=0.7, label='Avg in Sky')
    plt.plot(theta[2], 'r*', ms=13, alpha=0.7, label='No Avg')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if isinstance(meanV2, float):
        textinfig = r'V3$\mean$ = %0.2f    $\theta$ = %0.2f' % (meanV3, theta)
        ax1.annotate(textinfig, xy=(1.02, 0.35), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'theta_'+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'theta_'+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')


def stddev_plot(case, sigmaV2, sigmaV3, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plot in sigmaV2-sigmaV3 space of the 3 tests: averaging in pixel space, averaging on sky,
     and no averaging.
    Args:
        case               -- string, for example '491Scene1_rapid_real_bgFrac0.3'
        sigmaV2            -- list of 3 numpy array of sigma values of V2 for Tests 1, 2, and 3
        sigmaV3            -- list of 3 numpy array of sigma values of V3 for Tests 1, 2, and 3
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # check if the plot is for an Nk set
    basename = case
    if not isinstance(sigmaV2, float):
        basename = case+'_'+str(len(sigmaV2[0]))+'samples'

    # Make the plot of V3-THETA
    plot_title = r'Standard Deviations, $\sigma$'
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.suptitle(plot_title, fontsize=18, y=0.96)
    plt.title(basename)
    plt.xlabel(r'V2$\sigma$')
    plt.ylabel(r'V3$\sigma$')
    xmin, xmax = -0.005, 0.025
    ymin, ymax = -0.005, 0.025
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, xmin, xmax*2, colors='k', linestyles='dashed')
    plt.vlines(0.0, ymin, ymax*2, colors='k', linestyles='dashed')
    plt.plot(sigmaV2[0], sigmaV3[0], 'b^', ms=10, alpha=0.7, label='Avg in Pixel Space')
    plt.plot(sigmaV2[1], sigmaV3[1], 'go', ms=10, alpha=0.7, label='Avg in Sky')
    plt.plot(sigmaV2[2], sigmaV3[2], 'r*', ms=13, alpha=0.7, label='No Avg')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if isinstance(meanV2, float):
        textinfig = r'V2$\sigma$ = %0.2f    V3$\sigma$ = %0.2f' % (sigmaV2, sigmaV3)
        ax1.annotate(textinfig, xy=(1.02, 0.35), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'sigma_'+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'sigma_'+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')


#######################################################################################################################


if __name__ == '__main__':


    #### Set parameters

    centroid_windows = [3, 5, 7]
    Nsigma2plot = 2
    case = '491Scene1_rapid_real_bgFrac0.3'
    save_plot = False
    show_plot = True


    ######################################################

    # general path to text files
    gen_path = os.path.abspath('../resultsXrandomstars')

    # Loop over centroid_windows
    for cwin in centroid_windows:
        # load the data fom the 3 tests
        test_files_list = glob(os.path.join(gen_path, 'TEST*'+case+'*_Nsigma'+repr(Nsigma2plot)+'*'+repr(cwin)+'.txt'))
        #          0        1        2          3        4       5         6         7         8
        # data = sample, sigmaV2, sigmaV3, sigmaTheta, meanV2, meanV3, meanTheta, LastIter, RejStars
        dataT1 = np.loadtxt(test_files_list[0], comments='#', unpack=True)
        dataT2 = np.loadtxt(test_files_list[1], comments='#', unpack=True)
        dataT3 = np.loadtxt(test_files_list[2], comments='#', unpack=True)

        # compact variables
        meanV2 = [dataT1[4], dataT2[4], dataT3[4]]
        meanV3 = [dataT1[5], dataT2[5], dataT3[5]]
        sigmaV2 = [dataT1[1], dataT2[1], dataT3[1]]
        sigmaV3 = [dataT1[2], dataT2[2], dataT3[2]]
        theta = [dataT1[6], dataT2[6], dataT3[6]]
        cwincase = case+'_CentroidWindow'+repr(cwin)

        # Means plot
        v2v3mean_plot(cwincase, meanV2, meanV3, save_plot=save_plot, show_plot=show_plot, destination=None)

        # Thetas plot
        v2theta_plot(cwincase, meanV2, theta, save_plot=save_plot, show_plot=show_plot, destination=None)
        v3theta_plot(cwincase, meanV3, theta, save_plot=save_plot, show_plot=show_plot, destination=None)
        theta_plot(cwincase, theta, save_plot=save_plot, show_plot=show_plot, destination=None)

        # Standard Deviation plot
        stddev_plot(cwincase, sigmaV2, sigmaV3, save_plot=save_plot, show_plot=show_plot, destination=None)


