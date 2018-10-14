
from __future__ import print_function

import os, sys, unittest, importlib
sys.path.insert(0,".")

modules = [
    "test_sym_u",
    "testcol",
    "test_cImageD11",
    "gv_general.test_gv_general",
    "testcolumnfile",
    "test_put_incr",
    "testscale",
    "test_ubito",
    "test_uncomputegv",
    "test_transform",
    "twinprob.test_twin",
    "testcolfile2db.testcolfile2db",
    "testconnectedpixels.testconnectedpixels",
    "testlabelimage.testlabelimage",

]

if "all" in sys.argv:
    modules += ["index_demos.generate_gv",
                "test_index_unknown.test_index_unknown",
                "peaksearchtiftest.make_test_data"
    ]
else:
    print ("Add \"all\" to command line to run all tests")
HERE = os.getcwd()
print( "HERE",HERE )

for M in modules:
    os.chdir(HERE)
    print("Running suite for ",M)
    if M.find(".")>-1:
        path = M.split(".")
        for direc in path[:-1]:
            os.chdir( direc )
        M = path[-1]
    MOD = importlib.import_module(M)
    mySuite = unittest.loader.findTestCases( MOD )
    runner = unittest.TextTestRunner()
    try:
        runner.run(mySuite)
    except:
        raise
