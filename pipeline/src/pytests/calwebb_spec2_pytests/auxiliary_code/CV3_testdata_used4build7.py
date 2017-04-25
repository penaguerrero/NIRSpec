"""
This script contains a dictionary of the CV3 data used for testing build 7, as well
as the corresponding ESA intermediary products.
"""

CV3_testdata_dict = {}

# Fixed Slit test data
# ALLSLITS test data, this is in subarray mode 2048x256
CV3_testdata_dict["FS"]["PRISM"] = {
    "filter" : "CLEAR",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FFB",
    "NID" : 9852,
    "CV3filename" : ["NRSSLIT-COMBO-077_1_491_SE_2013-01-20T01h09m52.fits",
                     "NRSSLIT-COMBO-077_1_492_SE_2013-01-20T01h10m04.fits"],
    "level1Bfilenames" : ["jwtest1019001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1020001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["PRISM"] = {
    "filter" : "F100LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FFB",
    "NID" : 9845,
    "CV3filename" : ["NRSSLIT-COMBO-070_1_491_SE_2013-01-20T00h44m57.fits",
                     "NRSSLIT-COMBO-070_1_492_SE_2013-01-20T00h45m16.fits"],
    "level1Bfilenames" : ["jwtest1017001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1018001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G140M"] = {
    "filter" : "F070LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : ["CLS/FFV", "CLS/ARGON"],
    "NID" : [9814, 9812],
    "CV3filename" : [["NRSSLIT-COMBO-039_1_491_SE_2013-01-19T22h20m46.fits",
                     "NRSSLIT-COMBO-039_1_492_SE_2013-01-19T22h21m04.fits"],
                     ["NRSSLIT-COMBO-037_1_491_SE_2013-01-19T22h14m57.fits",
                      "NRSSLIT-COMBO-037_1_492_SE_2013-01-19T22h15m20.fits"]],
    "level1Bfilenames" : [["jwtest1011001_01101_00001_NRS1_uncal_mod.fits",
                           "jwtest1012001_01101_00001_NRS2_uncal_mod.fits"],
                          ["jwtest1009001_01101_00001_NRS1_uncal_mod.fits",
                           "jwtest1010001_01101_00001_NRS2_uncal_mod.fits"]]
}
CV3_testdata_dict["FS"]["G140H"] = {
    "filter" : "F070LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/ARGON",
    "NID" : 9776,
    "CV3filename" : ["NRSSLIT-COMBO-001_1_491_SE_2013-01-19T18h27m13.fits",
                     "NRSSLIT-COMBO-001_1_492_SE_2013-01-19T18h27m30.fits"],
    "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1002001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G140H"] = {
    "filter" : "F100LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF1",
    "NID" : 9791,
    "CV3filename" : ["NRSSLIT-COMBO-016_1_491_SE_2013-01-19T20h12m14.fits",
                     "NRSSLIT-COMBO-016_1_492_SE_2013-01-19T20h14m05.fits"],
    "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1004001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G235M"] = {
    "filter" : "F170LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF2",
    "NID" : 9830,
    "CV3filename" : ["NRSSLIT-COMBO-055_1_491_SE_2013-01-19T23h33m01.fits",
                     "NRSSLIT-COMBO-055_1_492_SE_2013-01-19T23h33m37.fits"],
    "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1014001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G235H"] = {
    "filter" : "F290LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF3",
    "NID" : 9799,
    "CV3filename" : ["NRSSLIT-COMBO-024_1_491_SE_2013-01-19T20h54m11.fits",
                     "NRSSLIT-COMBO-024_1_492_SE_2013-01-19T20h54m51.fits"],
    "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1006001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G395H"] = {
    "filter" : "F290LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF3",
    "NID" : 9801,
    "CV3filename" : ["NRSSLIT-COMBO-026_1_491_SE_2013-01-19T21h12m50.fits",
                     "NRSSLIT-COMBO-026_1_492_SE_2013-01-19T21h14m03.fits"],
    "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1008001_01101_00001_NRS2_uncal_mod.fits"]
}


# MOS test data
CV3_testdata_dict["MOS"]["G140H"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE1",
    "NID" : 41598,
    "CV3filename" : ["NRSV96214001001P0000000002105_1_491_SE_2016-01-24T01h59m01.fits ",
                     "NRSV96214001001P0000000002105_1_492_SE_2016-01-24T01h59m01.fits"],
    "level1Bfilenames" : ["jwtest1015001_01101_00001_NRS1_uncal.fits", "jwtest1016001_01101_00001_NRS1_uncal.fits"],
    "MSA_config" : "V9621400100101"
}
"""
CV3_testdata_dict["MOS"]["G140H"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE1",
    "NID" : 37286,
    "CV3filename" : ["NRSV00300060001P0000000002103_1_491_SE_2016-01-06T04h04m49.fits",
                     "NRSV00300060001P0000000002103_1_492_SE_2016-01-06T04h04m49.fits"],
    "level1Bfilenames" : ["", ""],
    "MSA_config" : "V0030006000101"
}
"""
CV3_testdata_dict["MOS"]["G140M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE1",
    "NID" : 39547,
    "CV3filename" : ["NRSV84600010001P0000000002101_4_491_SE_2016-01-17T17h34m08.fits",
                     "NRSV84600010001P0000000002101_4_492_SE_2016-01-17T17h34m08.fits"],
    "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal.fits", "jwtest1002001_01101_00001_NRS2_uncal.fits"],
    "MSA_config" : "V8460001000101"
}
CV3_testdata_dict["MOS"]["G235M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE2",
    "NID" : 39553,
    "CV3filename" : ["NRSV84600011001P0000000002101_2_491_SE_2016-01-17T18h18m48.fits",
                     "NRSV84600011001P0000000002101_2_492_SE_2016-01-17T18h18m48.fits"],
    "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits", "jwtest1006001_01101_00001_NRS2_uncal.fits"],
    "MSA_config" : "V8460001100101"
}
CV3_testdata_dict["MOS"]["G395M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE3",
    "NID" : 41543,
    "CV3filename" : ["NRSV96215001001P0000000002103_1_491_SE_2016-01-24T01h25m07.fits",
                     "NRSV96215001001P0000000002103_1_492_SE_2016-01-24T01h25m07.fits"],
    "level1Bfilenames" : ["jwtest1010001_01101_00001_NRS1_uncal.fits", "jwtest1009001_01101_00001_NRS2_uncal.fits"],
    "MSA_config" : "V9621500100101"
}
CV3_testdata_dict["MOS"]["PRISM"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE4",
    "NID" : 37328,
    "CV3filename" : ["NRSV00300060001P000000000210T_1_491_SE_2016-01-06T06h27m34.fits",
                     "NRSV00300060001P000000000210T_1_492_SE_2016-01-06T06h27m34.fits"],
    "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal.fits", "jwtest1014001_01101_00001_NRS1_uncal.fits"],
    "MSA_config" : "V0030006000104"
}


# IRS2 data has a raw physical size of 2048x3200
# In some places IRS2 considered subarray data only because it is not 2048x2048.
"""
CV3_testdata_dict["IRS2"]["G140H"] = {
    "filter" : "OPAQUE",
    "mode" : "MOS",
    "CAA_lamp" : "LINE1",
    "NID" : 30055,
    "CV3filename" : ["NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits",
                     "NRSSMOS-MOD-G1H-02-5344031756_1_492_SE_2015-12-10T03h25m56.fits"],
    "level1Bfilenames" : ["", ""],
    "MSA_config" : ""
}
"""
CV3_testdata_dict["IRS2"]["G235H"] = {
    "filter" : "OPAQUE",
    "mode" : "MOS",
    "CAA_lamp" : "LINE2",
    "NID" : 35373,
    "CV3filename" : ["NRSV00300010001P0000000002109_1_491_SE_2016-01-02T19h18m49.fits",
                     "NRSV00300010001P0000000002109_1_492_SE_2016-01-02T19h18m49.fits"],
    "level1Bfilenames" : ["jwtest1011001_01101_00001_NRS1_uncal.fits", "jwtest1012001_01101_00001_NRS1_uncal.fits"],
    "MSA_config" : "V0030001000101"
}


# IFU test data
CV3_testdata_dict["IFU"]["PRISM"] = {
    "filter" : "CLEAR",
    "CAA_lamp" : "NO_LAMP",
    "NID" : [37668, 37669],
    "CV3filename" : [["NRSSIMA-QUAL-04-B-6007022859_1_491_SE_2016-01-07T02h37m13.fits",
                      "NRSSIMA-QUAL-04-B-6007022859_1_492_SE_2016-01-07T02h37m13.fits"],
                     ["NRSSIMA-QUAL-04-2-6007023323_1_491_SE_2016-01-07T02h41m22.fits",
                      "NRSSIMA-QUAL-04-2-6007023323_1_492_SE_2016-01-07T02h41m22.fits"]],
    "level1Bfilenames" : [["jwtest1001001_01101_00001_NRS1_uncal.fits", "jwtest1002001_01101_00001_NRS2_uncal.fits"],
                          ["jwtest1003001_01101_00001_NRS2_uncal.fits", "jwtest1004001_01101_00001_NRS2_uncal.fits"]],
    "notes" : ["no external source, use for subtraction for NID37669,msa_config=ARDCLOSED",
               "very bright external source, use with NID37668,msa_config=ARDCLOSED"],
}
CV3_testdata_dict["IFU"]["G140M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE1",
    "NID" : 30192,
    "CV3filename" : ["NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06.fits",
                     "NRSSMOS-MOD-G1M-17-5344175105_1_492_SE_2015-12-10T18h00m05.fits"],
    "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits", "jwtest1006001_01101_00001_NRS1_uncal.fits"],
    "notes" : ""
}
CV3_testdata_dict["IFU"]["G235M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE2",
    "NID" : 30227,
    "CV3filename" : ["NRSSMOS-MOD-G2M-17-5344211451_1_491_SE_2015-12-10T21h23m36.fits",
                     "NRSSMOS-MOD-G2M-17-5344211451_1_492_SE_2015-12-10T21h23m37.fits"],
    "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal.fits", "jwtest1008001_01101_00001_NRS1_uncal.fits"],
    "notes" : ""
}
CV3_testdata_dict["IFU"]["G395M"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE3",
    "NID" : 30273,
    "CV3filename" : ["NRSSMOS-MOD-G3M-17-5345014854_1_491_SE_2015-12-11T01h57m10.fits",
                     "NRSSMOS-MOD-G3M-17-5345014854_1_492_SE_2015-12-11T01h57m10.fits"],
    "level1Bfilenames" : ["jwtest1009001_01101_00001_NRS1_uncal.fits", "jwtest1010001_01101_00001_NRS1_uncal.fits"],
    "notes" : ""
}


# DARKS
# the IRS2 data was taken with readout pattern NRSIRS2RAPID
CV3_testdata_dict["DARK"]["G140H"] = {
    "filter" : "OPAQUE",
    "readpatt": "NRSIRS2RAPID",
    "mode" : "IRS2",
    "CAA_lamp" : "NO_LAMP",
    #          1      2      3      4      5      6      7      8      9     10     11      12     13
    "NID" : [30487, 30489, 30491, 30494, 30496, 30500, 30506, 30507, 30508, 30513, 30514, 30577, 30555,
             30560, 30562, 30564, 30565, 30567, 30568, 30573, 30576, 30624, 30625, 30626, 30628, 30631],
    #          14     15     16     17     18     19     20     21     22     23     24     25     26
    "CV3filename" : [["NRSDET-DARK-IRS2-5345184144_1_491_SE_2015-12-11T19h03m38.fits", #1
                     "NRSDET-DARK-IRS2-5345184144_1_492_SE_2015-12-11T19h03m36.fits"],
                     ["NRSDET-DARK-IRS2-5345184144_2_491_SE_2015-12-11T19h20m47.fits", #2
                      "NRSDET-DARK-IRS2-5345184144_2_492_SE_2015-12-11T19h20m46.fits"],
                     ["NRSDET-DARK-IRS2-5345184144_3_491_SE_2015-12-11T19h38m37.fits", #3
                      "NRSDET-DARK-IRS2-5345184144_3_492_SE_2015-12-11T19h38m37.fits"],
                     ["NRSDET-DARK-IRS2-5345184144_4_491_SE_2015-12-11T19h55m00.fits", #4
                      "NRSDET-DARK-IRS2-5345184144_4_492_SE_2015-12-11T19h54m59.fits"],
                     ["NRSDET-DARK-IRS2-5345184144_5_491_SE_2015-12-11T20h12m29.fits", #5
                      "NRSDET-DARK-IRS2-5345184144_5_492_SE_2015-12-11T20h12m29.fits"],
                     ["NRSDET-DARK-IRS2-5345184144_6_491_SE_2015-12-11T20h32m39.fits", #6
                      "NRSDET-DARK-IRS2-5345184144_6_492_SE_2015-12-11T20h32m40.fits"],
                     ["NRSDET-DARK-IRS2-5345214649_1_491_SE_2015-12-11T22h09m20.fits", #7
                      "NRSDET-DARK-IRS2-5345214649_1_492_SE_2015-12-11T22h09m19.fits"],
                     ["NRSDET-DARK-IRS2-5345214649_2_491_SE_2015-12-11T22h26m19.fits", #8
                      "NRSDET-DARK-IRS2-5345214649_2_492_SE_2015-12-11T22h26m20.fits"],
                     ["NRSDET-DARK-IRS2-5345214649_3_491_SE_2015-12-11T22h43m19.fits", #9
                      "NRSDET-DARK-IRS2-5345214649_3_492_SE_2015-12-11T22h43m20.fits"],
                     ["NRSDET-DARK-IRS2-5345214649_4_491_SE_2015-12-11T23h00m49.fits", #10
                      "NRSDET-DARK-IRS2-5345214649_4_492_SE_2015-12-11T23h00m50.fits"],
                     ["NRSDET-DARK-IRS2-5345214649_5_491_SE_2015-12-11T23h19m32.fits", #11
                      "NRSDET-DARK-IRS2-5345214649_5_492_SE_2015-12-11T23h19m31.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_10_491_SE_2015-12-12T10h55m20.fits",#12
                      "NRSDET-DARK-IRS2-5346073733_10_492_SE_2015-12-12T10h55m20.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_1_491_SE_2015-12-12T07h59m04.fits", #13
                      "NRSDET-DARK-IRS2-5346073733_1_492_SE_2015-12-12T07h59m04.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_2_491_SE_2015-12-12T08h19m25.fits", #14
                      "NRSDET-DARK-IRS2-5346073733_2_492_SE_2015-12-12T08h19m24.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_3_491_SE_2015-12-12T08h36m35.fits", #15
                      "NRSDET-DARK-IRS2-5346073733_3_492_SE_2015-12-12T08h36m34.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_4_491_SE_2015-12-12T08h54m35.fits", #16
                      "NRSDET-DARK-IRS2-5346073733_4_492_SE_2015-12-12T08h54m35.fit"],
                     ["NRSDET-DARK-IRS2-5346073733_5_491_SE_2015-12-12T09h11m25.fits", #17
                      "NRSDET-DARK-IRS2-5346073733_5_492_SE_2015-12-12T09h11m24.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_6_491_SE_2015-12-12T09h30m45.fits", #18
                      "NRSDET-DARK-IRS2-5346073733_6_492_SE_2015-12-12T09h30m44.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_7_491_SE_2015-12-12T09h47m05.fits", #19
                      "NRSDET-DARK-IRS2-5346073733_7_492_SE_2015-12-12T09h47m05.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_8_491_SE_2015-12-12T10h04m02.fits", #20
                      "NRSDET-DARK-IRS2-5346073733_8_492_SE_2015-12-12T10h04m03.fits"],
                     ["NRSDET-DARK-IRS2-5346073733_9_491_SE_2015-12-12T10h21m33.fits", #21
                      "NRSDET-DARK-IRS2-5346073733_9_492_SE_2015-12-12T10h21m33.fits"],
                     ["NRSDET-DARK-IRS2-5346161822_1_491_SE_2015-12-12T16h48m34.fits", #22
                      "NRSDET-DARK-IRS2-5346161822_1_492_SE_2015-12-12T16h48m33.fits"],
                     ["NRSDET-DARK-IRS2-5346161822_2_491_SE_2015-12-12T17h05m54.fits", #23
                      "NRSDET-DARK-IRS2-5346161822_2_492_SE_2015-12-12T17h05m54.fits"],
                     ["NRSDET-DARK-IRS2-5346161822_3_491_SE_2015-12-12T17h21m54.fits", #24
                      "NRSDET-DARK-IRS2-5346161822_3_492_SE_2015-12-12T17h21m54.fits"],
                     ["NRSDET-DARK-IRS2-5346161822_4_491_SE_2015-12-12T17h39m43.fits", #25
                      "NRSDET-DARK-IRS2-5346161822_4_492_SE_2015-12-12T17h39m44.fits"],
                     ["NRSDET-DARK-IRS2-5346161822_5_491_SE_2015-12-12T18h58m53.fits", #26
                      "NRSDET-DARK-IRS2-5346161822_5_492_SE_2015-12-12T18h58m53.fits"]],
    "level1Bfilenames" : [["jwtest1001001_01101_00001_NRS1_uncal.fits", "jwtest1002001_01101_00001_NRS1_uncal.fits"], #1
                          ["jwtest1003001_01101_00001_NRS1_uncal.fits", "jwtest1004001_01101_00001_NRS1_uncal.fits"], #2
                          ["jwtest1005001_01101_00001_NRS1_uncal.fits", "jwtest1006001_01101_00001_NRS1_uncal.fits"], #3
                          ["jwtest1007001_01101_00001_NRS1_uncal.fits", "jwtest1008001_01101_00001_NRS1_uncal.fits"], #4
                          ["jwtest1009001_01101_00001_NRS1_uncal.fits", "jwtest1010001_01101_00001_NRS1_uncal.fits"], #5
                          ["jwtest1011001_01101_00001_NRS1_uncal.fits", "jwtest1012001_01101_00001_NRS1_uncal.fits"], #6
                          ["jwtest1013001_01101_00001_NRS1_uncal.fits", "jwtest1014001_01101_00001_NRS1_uncal.fits"], #7
                          ["jwtest1015001_01101_00001_NRS1_uncal.fits", "jwtest1016001_01101_00001_NRS1_uncal.fits"], #8
                          ["jwtest1017001_01101_00001_NRS1_uncal.fits", "jwtest1018001_01101_00001_NRS1_uncal.fits"], #9
                          ["jwtest1019001_01101_00001_NRS1_uncal.fits", "jwtest1020001_01101_00001_NRS1_uncal.fits"], #10
                          ["jwtest1021001_01101_00001_NRS1_uncal.fits", "jwtest1022001_01101_00001_NRS1_uncal.fits"], #11
                          ["jwtest1024001_01101_00001_NRS1_uncal.fits", "jwtest1023001_01101_00001_NRS1_uncal.fits"], #12
                          ["jwtest1025001_01101_00001_NRS1_uncal.fits", "jwtest1026001_01101_00001_NRS1_uncal.fits"], #13
                          ["jwtest1027001_01101_00001_NRS1_uncal.fits", "jwtest1028001_01101_00001_NRS1_uncal.fits"], #14
                          ["jwtest1029001_01101_00001_NRS1_uncal.fits", "jwtest1030001_01101_00001_NRS1_uncal.fits"], #15
                          ["jwtest1031001_01101_00001_NRS1_uncal.fits", "jwtest1032001_01101_00001_NRS1_uncal.fits"], #16
                          ["jwtest1033001_01101_00001_NRS1_uncal.fits", "jwtest1034001_01101_00001_NRS1_uncal.fits"], #17
                          ["jwtest1035001_01101_00001_NRS1_uncal.fits", "jwtest1036001_01101_00001_NRS1_uncal.fits"], #18
                          ["jwtest1037001_01101_00001_NRS1_uncal.fits", "jwtest1038001_01101_00001_NRS1_uncal.fits"], #19
                          ["jwtest1039001_01101_00001_NRS1_uncal.fits", "jwtest1040001_01101_00001_NRS1_uncal.fits"], #20
                          ["jwtest1041001_01101_00001_NRS1_uncal.fits", "jwtest1042001_01101_00001_NRS1_uncal.fits"], #21
                          ["jwtest1043001_01101_00001_NRS1_uncal.fits", "jwtest1044001_01101_00001_NRS1_uncal.fits"], #22
                          ["jwtest1045001_01101_00001_NRS1_uncal.fits", "jwtest1046001_01101_00001_NRS1_uncal.fits"], #23
                          ["jwtest1047001_01101_00001_NRS1_uncal.fits", "jwtest1048001_01101_00001_NRS1_uncal.fits"], #24
                          ["jwtest1049001_01101_00001_NRS1_uncal.fits", "jwtest1050001_01101_00001_NRS1_uncal.fits"], #25
                          ["jwtest1051001_01101_00001_NRS1_uncal.fits", "jwtest1052001_01101_00001_NRS1_uncal.fits"]],#26
    "MSA_config" : "ADRCLOSED"
}
CV3_testdata_dict["DARK"]["G140H"] = {
    "filter" : "OPAQUE",
    "readpatt": "NRSRAPID",
    "mode" : "MOS",
    "CAA_lamp" : "NO_LAMP",
    #          1      2      3      4      5      6      7      8      9     10     11      12     13
    "NID" : [30419, 30421, 30422, 30423, 30424, 30425, 30426, 30427, 30429, 30430, 30409, 30431, 30432,
             30433, 30434, 30435, 30441, 30410, 30411, 30413, 30414, 30415, 30416, 30417, 30418],
    #          14     15     16     17     18     19     20     21     22     23     24     25
    "CV3filename" : [["NRSDET-DARK-TRAD-5345123434_10_491_SE_2015-12-11T13h10m09.fits", #1
                     "NRSDET-DARK-TRAD-5345123434_10_492_SE_2015-12-11T13h10m09.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_11_491_SE_2015-12-11T13h15m19.fits", #2
                      "NRSDET-DARK-TRAD-5345123434_11_492_SE_2015-12-11T13h15m19.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_12_491_SE_2015-12-11T13h17m39.fits", #3
                      "NRSDET-DARK-TRAD-5345123434_12_492_SE_2015-12-11T13h17m39.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_13_491_SE_2015-12-11T13h19m49.fits", #4
                      "NRSDET-DARK-TRAD-5345123434_13_492_SE_2015-12-11T13h19m49.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_14_491_SE_2015-12-11T13h23m29.fits", #5
                      "NRSDET-DARK-TRAD-5345123434_14_492_SE_2015-12-11T13h23m29.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_15_491_SE_2015-12-11T13h26m19.fits", #6
                      "NRSDET-DARK-TRAD-5345123434_15_492_SE_2015-12-11T13h26m19.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_16_491_SE_2015-12-11T13h29m19.fits", #7
                      "NRSDET-DARK-TRAD-5345123434_16_492_SE_2015-12-11T13h29m19.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_17_491_SE_2015-12-11T13h31m39.fits", #8
                      "NRSDET-DARK-TRAD-5345123434_17_492_SE_2015-12-11T13h31m39.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_18_491_SE_2015-12-11T13h34m09.fits", #9
                      "NRSDET-DARK-TRAD-5345123434_18_492_SE_2015-12-11T13h34m09.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_19_491_SE_2015-12-11T13h36m29.fits", #10
                      "NRSDET-DARK-TRAD-5345123434_19_492_SE_2015-12-11T13h36m29.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_1_491_SE_2015-12-11T12h49m00.fits", #11
                      "NRSDET-DARK-TRAD-5345123434_1_492_SE_2015-12-11T12h49m00.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_20_491_SE_2015-12-11T13h39m02.fits",#12
                      "NRSDET-DARK-TRAD-5345123434_20_492_SE_2015-12-11T13h39m02.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_21_491_SE_2015-12-11T13h41m52.fits", #13
                      "NRSDET-DARK-TRAD-5345123434_21_492_SE_2015-12-11T13h41m51.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_22_491_SE_2015-12-11T13h43m51.fits", #14
                      "NRSDET-DARK-TRAD-5345123434_22_492_SE_2015-12-11T13h43m51.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_23_491_SE_2015-12-11T13h46m11.fits", #15
                      "NRSDET-DARK-TRAD-5345123434_23_492_SE_2015-12-11T13h46m11.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_24_491_SE_2015-12-11T13h48m42.fits", #16
                      "NRSDET-DARK-TRAD-5345123434_24_492_SE_2015-12-11T13h48m42.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_25_491_SE_2015-12-11T15h07m32.fits", #17
                      "NRSDET-DARK-TRAD-5345123434_25_492_SE_2015-12-11T15h07m32.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_2_491_SE_2015-12-11T12h51m29.fits", #18
                      "NRSDET-DARK-TRAD-5345123434_2_492_SE_2015-12-11T12h51m29.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_3_491_SE_2015-12-11T12h53m39.fits", #19
                      "NRSDET-DARK-TRAD-5345123434_3_492_SE_2015-12-11T12h53m39.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_4_491_SE_2015-12-11T12h55m59.fits", #20
                      "NRSDET-DARK-TRAD-5345123434_4_492_SE_2015-12-11T12h55m59.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_5_491_SE_2015-12-11T12h58m29.fits", #21
                      "NRSDET-DARK-TRAD-5345123434_5_492_SE_2015-12-11T12h58m29.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_6_491_SE_2015-12-11T13h00m39.fits", #22
                      "NRSDET-DARK-TRAD-5345123434_6_492_SE_2015-12-11T13h00m39.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_7_491_SE_2015-12-11T13h02m59.fits", #23
                      "NRSDET-DARK-TRAD-5345123434_7_492_SE_2015-12-11T13h02m59.fit"],
                     ["NRSDET-DARK-TRAD-5345123434_8_491_SE_2015-12-11T13h05m19.fits", #24
                      "NRSDET-DARK-TRAD-5345123434_8_492_SE_2015-12-11T13h05m19.fits"],
                     ["NRSDET-DARK-TRAD-5345123434_9_491_SE_2015-12-11T13h07m39.fits", #25
                      "NRSDET-DARK-TRAD-5345123434_9_492_SE_2015-12-11T13h07m39.fit"]],
    "level1Bfilenames" : [["jwtest1001001_01101_00001_NRS1_uncal.fits", "jwtest1002001_01101_00001_NRS1_uncal.fits"], #1
                          ["jwtest1003001_01101_00001_NRS1_uncal.fits", "jwtest1004001_01101_00001_NRS1_uncal.fits"], #2
                          ["jwtest1005001_01101_00001_NRS1_uncal.fits", "jwtest1006001_01101_00001_NRS1_uncal.fits"], #3
                          ["jwtest1007001_01101_00001_NRS1_uncal.fits", "jwtest1008001_01101_00001_NRS1_uncal.fits"], #4
                          ["jwtest1009001_01101_00001_NRS1_uncal.fits", "jwtest1010001_01101_00001_NRS1_uncal.fits"], #5
                          ["jwtest1011001_01101_00001_NRS1_uncal.fits", "jwtest1012001_01101_00001_NRS1_uncal.fits"], #6
                          ["jwtest1013001_01101_00001_NRS1_uncal.fits", "jwtest1014001_01101_00001_NRS1_uncal.fits"], #7
                          ["jwtest1015001_01101_00001_NRS1_uncal.fits", "jwtest1016001_01101_00001_NRS1_uncal.fits"], #8
                          ["jwtest1017001_01101_00001_NRS1_uncal.fits", "jwtest1018001_01101_00001_NRS1_uncal.fits"], #9
                          ["jwtest1019001_01101_00001_NRS1_uncal.fits", "jwtest1020001_01101_00001_NRS1_uncal.fits"], #10
                          ["jwtest1021001_01101_00001_NRS1_uncal.fits", "jwtest1022001_01101_00001_NRS1_uncal.fits"], #11
                          ["jwtest1023001_01101_00001_NRS1_uncal.fits", "jwtest1024001_01101_00001_NRS1_uncal.fits"], #12
                          ["jwtest1025001_01101_00001_NRS1_uncal.fits", "jwtest1026001_01101_00001_NRS1_uncal.fits"], #13
                          ["jwtest1027001_01101_00001_NRS1_uncal.fits", "jwtest1028001_01101_00001_NRS1_uncal.fits"], #14
                          ["jwtest1029001_01101_00001_NRS1_uncal.fits", "jwtest1030001_01101_00001_NRS1_uncal.fits"], #15
                          ["jwtest1031001_01101_00001_NRS1_uncal.fits", "jwtest1032001_01101_00001_NRS1_uncal.fits"], #16
                          ["jwtest1033001_01101_00001_NRS1_uncal.fits", "jwtest1034001_01101_00001_NRS1_uncal.fits"], #17
                          ["jwtest1035001_01101_00001_NRS1_uncal.fits", "jwtest1036001_01101_00001_NRS1_uncal.fits"], #18
                          ["jwtest1037001_01101_00001_NRS1_uncal.fits", "jwtest1038001_01101_00001_NRS1_uncal.fits"], #19
                          ["jwtest1039001_01101_00001_NRS1_uncal.fits", "jwtest1040001_01101_00001_NRS1_uncal.fits"], #20
                          ["jwtest1041001_01101_00001_NRS1_uncal.fits", "jwtest1042001_01101_00001_NRS1_uncal.fits"], #21
                          ["jwtest1043001_01101_00001_NRS1_uncal.fits", "jwtest1044001_01101_00001_NRS1_uncal.fits"], #22
                          ["jwtest1045001_01101_00001_NRS1_uncal.fits", "jwtest1046001_01101_00001_NRS1_uncal.fits"], #23
                          ["jwtest1047001_01101_00001_NRS1_uncal.fits", "jwtest1048001_01101_00001_NRS1_uncal.fits"], #24
                          ["jwtest1049001_01101_00001_NRS1_uncal.fits", "jwtest1050001_01101_00001_NRS1_uncal.fits"]],#25
    "MSA_config" : "ADRCLOSED"
}
