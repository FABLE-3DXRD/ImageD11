
import os, sys

"""
Usage : python install_ImageD11_from_git.py [path] [checkout=ImageD11]

  - runs a git clone of the ImageD11 source into folder path/checkout
  - runs setup.py build_ext --inplace  (e.g. an editable install)
  - inserts the code path into the front of sys.path
  
  To be called in ipython via:
  
  In python do:
  
"""


def run_ImageD11_from_git( path, checkout ):
    if not os.path.exists( path ):
        print("Creating")
        os.makedirs( path )
    assert os.path.isdir( path )
    code_path = os.path.join(path, checkout)
    if not os.path.exists(code_path):
        os.system('git clone https://github.com/FABLE-3DXRD/ImageD11 ' + code_path)
        assert os.path.exists(code_path),'failed to checkout from git'
    bld = os.path.join(code_path, "build")
    if not os.path.exists(bld):
        os.system('cd '+code_path+' && python setup.py build_ext --inplace')
        assert os.path.exists(bld), 'failed to compile'
    print("# Setting path via: ")
    sys.path.insert(0, code_path)
    print("sys.path.insert(0,",code_path,")")
    import ImageD11, ImageD11.cImageD11
    print("# Running from:", ImageD11.__file__)
    return code_path
    
    
def guess_ImageD11_code_path():  # This should be in silx somewhere?
    """ Locates:
    scripts can hold a local installation of ImageD11 if you need one
    """
    path_items = os.getcwd().split('/')
    if 'visitor' in path_items:
        vi = path_items.index( 'visitor' )
        experiment, session = path_items[ vi + 1], path_items[ vi + 3 ]
        return os.path.join( "/data", "visitor", experiment, "id11", session, "SCRIPTS" )
    return os.environ['HOME']+'/Code'


def setup_ImageD11_from_git( path=None, checkout="ImageD11" ):
    if path is None:
        return run_ImageD11_from_git( guess_ImageD11_code_path(), checkout )
    else:
        return run_ImageD11_from_git( path, checkout )

