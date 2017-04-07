

import os, sys, unittest, importlib
sys.path.insert(0,".")

modules = [
    "gv_general.test_gv_general",
    "testcolumnfile",
    "test_put_incr",
    "testscale",
    "test_ubito",
    "test_uncomputegv",
    "test_transform",
    "index_demos.generate_gv",
    "peaksearchtiftest.make_test_data",
    "test_index_unknown.test_index_unknown",
    "testcolfile2db.testcolfile2db",
    "testconnectedpixels.testconnectedpixels",
    "testlabelimage.testlabelimage",
       ]

HERE = os.getcwd()
print "HERE",HERE

for M in modules:
    os.chdir(HERE)
    print(M)
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
