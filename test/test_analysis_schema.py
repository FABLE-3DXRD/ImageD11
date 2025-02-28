import os

import numpy as np
import unittest

import ImageD11.parameters
import ImageD11.unitcell


class TestFromDefault(unittest.TestCase):
    def test_from_default(self):
        """Test that we can create default parameters which get saved to disk"""
        asc = ImageD11.parameters.AnalysisSchema.from_default(detector='eiger')
        asc.save('pars.json')
        # make sure files are saved to disk
        self.assertTrue(os.path.exists('geometry.par'))
        self.assertTrue(os.path.exists('pars.json'))
        self.assertTrue(os.path.exists('CeO2.par'))
        # make sure the content of each file matches our asc dicts
        old_pars_dict = ImageD11.parameters.read_par_file('geometry.par').get_parameters()
        new_pars_dict = asc.get_xfab_pars_dict()
        for key, value in old_pars_dict.items():
            if key in new_pars_dict.keys():
                if isinstance(value, (float, int)):
                    self.assertTrue(np.allclose(value, new_pars_dict[key]))
                else:
                    self.assertEqual(value, new_pars_dict[key])
        old_pars_dict = ImageD11.parameters.read_par_file('CeO2.par').get_parameters()
        new_pars_dict = asc.phase_pars_obj_dict['CeO2'].get_parameters()
        for key, value in old_pars_dict.items():
            if key in new_pars_dict.keys():
                if isinstance(value, (float, int)):
                    self.assertTrue(np.allclose(value, new_pars_dict[key]))
                else:
                    self.assertEqual(value, new_pars_dict[key])

    @classmethod
    def tearDownClass(cls):
        """clean up files after being made"""
        os.remove('geometry.par')
        os.remove('pars.json')
        os.remove('CeO2.par')


class TestAddPhase(unittest.TestCase):
    def test_add_phase(self):
        """Test that we can add a phase"""
        asc = ImageD11.parameters.AnalysisSchema.from_default(detector='eiger')
        new_ucell = ImageD11.unitcell.unitcell([1., 2., 3., 4., 5., 6], 1)
        asc.add_phase_from_unitcell('test_phase', new_ucell)
        # save to disk
        asc.save('pars.json')
        # make sure files are saved to disk
        self.assertTrue(os.path.exists('test_phase.par'))
        # now if we import pars.json as a Phases object, can we access the new unitcell?
        phases = ImageD11.unitcell.Phases('pars.json')
        self.assertTrue(np.allclose(phases.unitcells['test_phase'].lattice_parameters, new_ucell.lattice_parameters))

    @classmethod
    def tearDownClass(cls):
        """clean up files after being made"""
        os.remove('geometry.par')
        os.remove('pars.json')
        os.remove('CeO2.par')
        os.remove('test_phase.par')


class TestToOldParsFile(unittest.TestCase):
    def test_to_old(self):
        """Test that we can make an old-style parameter file"""
        asc = ImageD11.parameters.AnalysisSchema.from_default(detector='eiger')
        asc.to_old_pars_file('oldpars.par', 'CeO2')
        self.assertTrue(os.path.exists('oldpars.par'))
        # contents in the old pars file should match the asc
        old_pars_dict = ImageD11.parameters.read_par_file('oldpars.par').get_parameters()
        new_pars_dict = asc.get_xfab_pars_dict(phase_name='CeO2')
        for key, value in old_pars_dict.items():
            if key in new_pars_dict.keys():
                if isinstance(value, (float, int)):
                    self.assertTrue(np.allclose(value, new_pars_dict[key]))
                else:
                    self.assertEqual(value, new_pars_dict[key])

    @classmethod
    def tearDownClass(cls):
        """clean up files after being made"""
        os.remove('oldpars.par')


class TestParsFollowJson(unittest.TestCase):
    def test(self):
        # make an old pars file
        asc = ImageD11.parameters.AnalysisSchema.from_default(detector='eiger')
        asc.to_old_pars_file('oldpars.par', 'CeO2')
        # load it in
        asc = ImageD11.parameters.AnalysisSchema.from_old_pars_file('oldpars.par', phase_name='CeO2')
        os.mkdir('new_folder')
        asc.save(os.path.join('new_folder', 'pars.json'))
        self.assertTrue(os.path.exists(os.path.join('new_folder', 'pars.json')))
        self.assertTrue(os.path.exists(os.path.join('new_folder', 'geometry.par')))
        self.assertTrue(os.path.exists(os.path.join('new_folder', 'CeO2.par')))

    @classmethod
    def tearDownClass(cls):
        """clean up files after being made"""
        os.remove(os.path.join('new_folder', 'pars.json'))
        os.remove(os.path.join('new_folder', 'geometry.par'))
        os.remove(os.path.join('new_folder', 'CeO2.par'))
        os.rmdir('new_folder')
