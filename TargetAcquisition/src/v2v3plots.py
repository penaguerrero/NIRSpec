from __future__ import print_function, division
import matplotlib.pyplot as plt
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


def make_v2v3plots(case, T_diffVs, LS_res, save_plot=False, show_plot=False, destination=None):
    """
    This function creates the plots in V2-V3 space.
    Args:
        case               -- string, for example 'Scene2_rapid_real_bgFrac'
        T_diffVs           -- list of 6 lists corresponding to differences (true-measured) of V2 and V3 for centroid windows 3, 5, and 7
        LS_res             -- list of 6 lists corresponding to sigmas and means of centroid windows 3, 5, and 7
        save_plot          -- True or False
        show_plot          -- True or False
        destination        -- string, destination directory
    Returns:

    """
    # Set the paths
    results_path = os.path.abspath('../plots4presentationIST')

    # unfold variables
    diffV2_3, diffV3_3, diffV2_5, diffV3_5, diffV2_7, diffV3_7 = T_diffVs
    LSsigmas_3, LSsigmas_5, LSsigmas_7, LSdeltas_3, LSdeltas_5, LSdeltas_7 = LS_res

    # check if the plot is for an Nk set
    basename = ''
    if isinstance(LSsigmas_3[0], list):
        basename = '_'+str(len(LSsigmas_3[0]))+'set'

    # Make the plot of MEANS
    plot_title = 'Mean Residual Values for '+case
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.title(plot_title)
    plt.xlabel(r'$\Delta$V2')
    plt.ylabel(r'$\Delta$V3')
    plt.xlim(0.0, 0.4)
    plt.ylim(0.0, 0.1)
    plt.plot(LSdeltas_3[0], LSdeltas_3[1], 'b^', ms=10, alpha=0.7, label='Centroid window=3')
    plt.plot(LSdeltas_5[0], LSdeltas_5[1], 'go', ms=10, alpha=0.7, label='Centroid window=5')
    plt.plot(LSdeltas_7[0], LSdeltas_7[1], 'r*', ms=13, alpha=0.7, label='Centroid window=7')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if not isinstance(LSsigmas_3[0], list):
        textinfig3 = r'V2$\mu3$ = %0.2f    V3$\mu3$ = %0.2f' % (LSdeltas_3[0], LSdeltas_3[1])
        textinfig5 = r'V2$\mu5$ = %0.2f    V3$\mu5$ = %0.2f' % (LSdeltas_5[0], LSdeltas_5[1])
        textinfig7 = r'V2$\mu7$ = %0.2f    V3$\mu7$ = %0.2f' % (LSdeltas_7[0], LSdeltas_7[1])
        ax1.annotate(textinfig3, xy=(1.02, 0.35), xycoords='axes fraction' )
        ax1.annotate(textinfig5, xy=(1.02, 0.32), xycoords='axes fraction' )
        ax1.annotate(textinfig7, xy=(1.02, 0.29), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'means_'+case+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'means_'+case+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')

    # Make the plot of STANDARD DEVIATIONS
    plot_title = 'Standard Deviation of Residuals for '+case
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.title(plot_title)
    plt.xlabel(r'$\Delta$V2')
    plt.ylabel(r'$\Delta$V3')
    plt.xlim(0.0, 0.35)
    plt.ylim(0.0, 0.25)
    plt.plot(LSsigmas_3[0], LSsigmas_3[1], 'b^', ms=10, alpha=0.7, label='Centroid window=3')
    plt.plot(LSsigmas_5[0], LSsigmas_5[1], 'go', ms=10, alpha=0.7, label='Centroid window=5')
    plt.plot(LSsigmas_7[0], LSsigmas_3[1], 'r*', ms=13, alpha=0.7, label='Centroid window=7')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if not isinstance(LSsigmas_3[0], list):
        textinfig3 = r'V2$\sigma3$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_3[0], LSsigmas_3[1])
        textinfig5 = r'V2$\sigma5$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_5[0], LSsigmas_5[1])
        textinfig7 = r'V2$\sigma7$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_7[0], LSsigmas_7[1])
        ax1.annotate(textinfig3, xy=(1.02, 0.35), xycoords='axes fraction' )
        ax1.annotate(textinfig5, xy=(1.02, 0.32), xycoords='axes fraction' )
        ax1.annotate(textinfig7, xy=(1.02, 0.29), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'stddevs_'+case+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'stddevs_'+case+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')

    # Make the plot showing distribution
    plot_title = 'Distribution of Residuals for '+case
    fig1 = plt.figure(1, figsize=(12, 10))
    ax1 = fig1.add_subplot(111)
    plt.title(plot_title)
    plt.xlabel(r'$\Delta$V2')
    plt.ylabel(r'$\Delta$V3')
    xmin, xmax = -0.7, 0.7
    ymin, ymax = -0.4, 0.4
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.hlines(0.0, xmin, xmax*2, colors='k', linestyles='dashed')
    plt.vlines(0.0, ymin, ymax*2, colors='k', linestyles='dashed')
    plt.plot(diffV2_3, diffV3_3, 'b^', ms=6, alpha=0.7, label='Centroid window=3')
    plt.plot(diffV2_5, diffV3_5, 'go', ms=4, alpha=0.7, label='Centroid window=5')
    plt.plot(diffV2_7, diffV3_7, 'r*', ms=6, alpha=0.7, label='Centroid window=7')
    # Shrink current axis by 10%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))   # put legend out of the plot box
    if not isinstance(LSsigmas_3[0], list):
        textinfig3 = r'V2$\sigma3$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_3[0], LSsigmas_3[1])
        textinfig5 = r'V2$\sigma5$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_5[0], LSsigmas_5[1])
        textinfig7 = r'V2$\sigma7$ = %0.2f    V3$\sigma3$ = %0.2f' % (LSsigmas_7[0], LSsigmas_7[1])
        textinfig3b = r'V2$\mu3$ = %0.2f    V3$\mu3$ = %0.2f' % (LSdeltas_3[0], LSdeltas_3[1])
        textinfig5b = r'V2$\mu5$ = %0.2f    V3$\mu5$ = %0.2f' % (LSdeltas_5[0], LSdeltas_5[1])
        textinfig7b = r'V2$\mu7$ = %0.2f    V3$\mu7$ = %0.2f' % (LSdeltas_7[0], LSdeltas_7[1])
        ax1.annotate(textinfig3, xy=(1.02, 0.35), xycoords='axes fraction' )
        ax1.annotate(textinfig5, xy=(1.02, 0.32), xycoords='axes fraction' )
        ax1.annotate(textinfig7, xy=(1.02, 0.29), xycoords='axes fraction' )
        ax1.annotate(textinfig3b, xy=(1.02, 0.23), xycoords='axes fraction' )
        ax1.annotate(textinfig5b, xy=(1.02, 0.20), xycoords='axes fraction' )
        ax1.annotate(textinfig7b, xy=(1.02, 0.17), xycoords='axes fraction' )
    if save_plot:
        if destination is not None:
            fig_name = os.path.join(destination, 'residuals_'+case+basename+'.jpg')
        else:
            fig_name = os.path.join(results_path, 'residuals_'+case+basename+'.jpg')
        fig1.savefig(fig_name)
        print ("\n Plot saved: ", fig_name)
    if show_plot:
        plt.show()
    else:
        plt.close('all')
