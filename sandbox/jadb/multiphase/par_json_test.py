# exec(open('/data/id11/nanoscope/install_ImageD11_from_git.py').read())

# Add ImageD11 local
# PYTHONPATH = setup_ImageD11_from_git( os.path.join( os.environ['HOME'],'Code'), 'ImageD11' )

from ImageD11.parameters import parameters, AnalysisSchema
from ImageD11.unitcell import Phases
from ImageD11.columnfile import columnfile


def main():
    # read a json file into AnalysisSchema
    json_pars = AnalysisSchema('pars.json')

    # print the json dict
    print(json_pars.pars_dict)

    # print the geometry pars object
    print(json_pars.geometry_pars_obj.get_parameters())

    # print the phase pars object dict
    print(json_pars.phase_pars_obj_dict)

    # print a specific phase pars object
    print(json_pars.phase_pars_obj_dict["ferrite"].get_parameters())

    # print a combined pars dict
    print(json_pars.xfab_pars_dict)

    # read a json file directly with parameters.loadparameters
    pars = parameters()
    pars.loadparameters(filename='pars.json')

    # print the pars
    print(pars.get_parameters())

    # load a columnfile
    cf = columnfile('peaks.flt')

    # print the header parameters
    print(cf.parameters.get_parameters())

    # print the lattice (should be 'I' for the header parameters)
    print(cf.parameters.get('cell_lattice_[P,A,B,C,I,F,R]'))

    # load the cf parameters from disk
    cf.parameters.loadparameters('pars.json')

    # print the lattice (should be 229 for the JSON parameters)
    print(cf.parameters.get('cell_lattice_[P,A,B,C,I,F,R]'))

    # import the phases
    phases = Phases('pars.json')

    # print the unitcells dict
    print(phases.unitcells)


if __name__ == '__main__':
    main()
