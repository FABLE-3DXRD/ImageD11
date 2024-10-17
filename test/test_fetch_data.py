"""Unit tests for ImageD11/fetch_data.py"""
import os
import shutil
import unittest

import six
from six.moves import urllib

import numpy as np

from ImageD11 import fetch_data
from ImageD11.sinograms import dataset as id11dset


if six.PY2:  # python 2/3 compatibility
    FileNotFoundError = IOError


class TestDownloadURL(unittest.TestCase):
    def setUp(self):
        valid_dsetname = 'Si_cube_S3DXRD_nt_moves_dty'
        valid_base_url = fetch_data.dataset_base_urls[valid_dsetname]
        valid_parfilename = fetch_data.dataset_filenames[valid_dsetname]['geomfile']
        self.valid_url = os.path.join(valid_base_url, valid_parfilename)
        self.valid_dest = './geometry.par'
        if os.path.exists(self.valid_dest):
            raise IOError('Intended download destination already exists! Aborting test...')
        self.invalid_folder = './djasqowdonqwfioqfioqnwd'  # unlikely to collide?
        if os.path.exists(self.invalid_folder):
            raise IOError('Intended not-working path actually exists! Aborting test...')

    # this is just testing urllib itself, which is weird on py2 vs 3, so don't do this for now
    # def test_invalid_url(self):
    #     invalid_url = "broken url"
    #
    #     with self.assertRaises(ValueError):
    #         fetch_data.download_url(invalid_url, self.valid_dest)

    def test_invalid_path(self):
        invalid_path = os.path.join(self.invalid_folder, 'pars.par') # parent folder doesn't exist

        with self.assertRaises(IOError):
            fetch_data.download_url(self.valid_url, invalid_path)

    def test_valid_download(self):

        # download the file
        fetch_data.download_url(self.valid_url, self.valid_dest)

        # compare contents

        test_file_contents = """chi 0.0
distance 148420.9411459722
fit_tolerance 0.1
min_bin_prob 1e-05
no_bins 1000000
o11 -1
o12 0
o21 0
o22 -1
omegasign 1.0
t_x 0
t_y 0
t_z 0
tilt_x 0.0
tilt_y -0.0018224188997752666
tilt_z 0.002220060473496065
wavelength 0.1927
wedge 0.0
weight_hist_intensities 0
y_center 1028.0096389590503
y_size 75.0
z_center 1111.1810324237968
z_size 75.0
"""

        self.assertEqual(test_file_contents, open(self.valid_dest).read())

    def tearDown(self):
        if os.path.exists(self.valid_dest):
            os.remove(self.valid_dest)


class TestGetDataset(unittest.TestCase):
    def setUp(self):
        self.valid_dsetname = 'Si_cube_S3DXRD_nt_moves_dty'
        self.valid_base_url = fetch_data.dataset_base_urls[self.valid_dsetname]
        self.valid_folder = './test_fetch_data'
        if os.path.exists(self.valid_folder):
            raise IOError('Intended download destination already exists! Aborting test...')
        self.invalid_folder = './djasqowdonqwfioqfioqnwd'  # unlikely to collide?
        if os.path.exists(self.invalid_folder):
            raise IOError('Intended not-working path actually exists! Aborting test...')

    def test_invalid_dsetname(self):
        invalid_dsetname = '128108nsjfjsdojaspixpkmsx'
        with self.assertRaises(ValueError):
            fetch_data._get_dataset(invalid_dsetname, self.valid_folder, allow_download=False)

    def test_valid_but_missing_folder(self):
        with self.assertRaises(FileNotFoundError):
            fetch_data._get_dataset(self.valid_dsetname, self.valid_folder, allow_download=False)

    def test_dataset_exists(self):
        # manually make the dataset
        # then check that the function can find it

        os.mkdir(self.valid_folder)
        raw_data_root_dir = os.path.join(self.valid_folder, 'raw')
        processed_data_root_dir = os.path.join(self.valid_folder, 'processed')

        ds = id11dset.DataSet(dataroot=raw_data_root_dir,
                              analysisroot=processed_data_root_dir,
                              sample=fetch_data.dataset_metadata[self.valid_dsetname]['sample'],
                              dset=fetch_data.dataset_metadata[self.valid_dsetname]['dataset'])

        # set a custom flag so we know we're saving it to disk and reloading it
        ds.omega = np.array([[0., 1], [2, 3]])
        ds.dty = np.array([[0., 1], [2, 3]])
        ds.pbpfile = 'Hello! I am testing this dataset function.'
        ds.save()

        ds_returned = fetch_data._get_dataset(self.valid_dsetname, self.valid_folder, allow_download=False)
        self.assertEqual(ds_returned.pbpfile, ds.pbpfile)

    def test_sparse_exists(self):
        # manually make a dataset
        # manually download a sparse file for it
        # download a par file too
        # check that function makes a dataset with the sparse file and par file included

        os.mkdir(self.valid_folder)
        raw_data_root_dir = os.path.join(self.valid_folder, 'raw')
        processed_data_root_dir = os.path.join(self.valid_folder, 'processed')

        ds = id11dset.DataSet(dataroot=raw_data_root_dir,
                              analysisroot=processed_data_root_dir,
                              sample=fetch_data.dataset_metadata[self.valid_dsetname]['sample'],
                              dset=fetch_data.dataset_metadata[self.valid_dsetname]['dataset'])

        os.makedirs(ds.analysispath)

        # download the sparse file
        urllib.request.urlretrieve(self.valid_base_url + fetch_data.dataset_filenames[self.valid_dsetname]['sparsefile'], ds.sparsefile)
        # download the parfile
        urllib.request.urlretrieve(
            self.valid_base_url + fetch_data.dataset_filenames[self.valid_dsetname]['sparsefile'], os.path.join(processed_data_root_dir, fetch_data.dataset_filenames[self.valid_dsetname]['parfile']))

        ds_returned = fetch_data._get_dataset(self.valid_dsetname, self.valid_folder, allow_download=False)
        # this dataset should have the right sparse and par file paths
        # should also have 41 scans
        self.assertEqual(ds_returned.sparsefile, os.path.abspath(ds.sparsefile))
        self.assertEqual(ds_returned.parfile, os.path.abspath(os.path.join(processed_data_root_dir, fetch_data.dataset_filenames[self.valid_dsetname]['parfile'])))
        self.assertEqual(len(ds_returned.scans), 41)
        # should also appear on disk
        self.assertTrue(os.path.exists(ds_returned.dsfile))

    def test_download(self):
        os.mkdir(self.valid_folder)
        ds_returned = fetch_data._get_dataset(self.valid_dsetname, self.valid_folder, allow_download=True)
        # check all files downloaded

        # where should the files have gone?
        parfile = os.path.abspath(os.path.join(self.valid_folder, 'processed', fetch_data.dataset_filenames[self.valid_dsetname]['parfile']))
        e2dxfile = os.path.abspath(os.path.join(self.valid_folder, 'processed', fetch_data.dataset_filenames[self.valid_dsetname]['e2dxfile']))
        e2dyfile = os.path.abspath(os.path.join(self.valid_folder, 'processed', fetch_data.dataset_filenames[self.valid_dsetname]['e2dyfile']))

        raw_data_root_dir = os.path.join(self.valid_folder, 'raw')
        processed_data_root_dir = os.path.join(self.valid_folder, 'processed')

        # where should the dataset and sparse files have gone?
        # use dataset logic for this
        ds = id11dset.DataSet(dataroot=raw_data_root_dir,
                              analysisroot=processed_data_root_dir,
                              sample=fetch_data.dataset_metadata[self.valid_dsetname]['sample'],
                              dset=fetch_data.dataset_metadata[self.valid_dsetname]['dataset'])

        self.assertTrue(os.path.exists(parfile))
        self.assertTrue(os.path.exists(e2dxfile))
        self.assertTrue(os.path.exists(e2dyfile))
        self.assertTrue(os.path.exists(ds.sparsefile))
        self.assertTrue(os.path.exists(ds.dsfile))

        # does the returned dataset object have these attributes?

        self.assertEqual(ds_returned.parfile, parfile)
        self.assertEqual(ds_returned.e2dxfile, e2dxfile)
        self.assertEqual(ds_returned.e2dyfile, e2dyfile)
        self.assertEqual(ds_returned.sparsefile, os.path.abspath(ds.sparsefile))
        self.assertEqual(ds_returned.dsfile, os.path.abspath(ds.dsfile))

        # does the returned dataset object have 41 scans? (from sparse)
        self.assertEqual(len(ds_returned.scans), 41)

    def tearDown(self):
        if os.path.exists(self.valid_folder):
            shutil.rmtree(self.valid_folder)
