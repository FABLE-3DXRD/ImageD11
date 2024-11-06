# interact with the sample environment
# 1) ADMET tensile rig for now
# 2) ...
# Haixing Fang, haixing.fang@esrf.fr
# May 28, 2024

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
import xml.etree.ElementTree as ET


    
class admet:
    '''
    read data from ADMET stress rig, which is normally save as .mtxData file
    example:
    filename = '/home/esrf/haixing0a/Downloads/ADMET_testing_20240327/ma5607_HEA2.mtwData'
    stress = admet(filename)
    '''
    def __init__(self, filename = None, outname = None):
        self.filename = filename
        self.header = "Position(mm), Load(N)"
        self.outname = outname
        if filename is not None:
            self.data, self.key_names = self.read_data()
        else:
            self.data = None
            self.key_names = None
    
    def read_data(self):
        if os.path.exists(self.filename):
            sample = self.parse_xml_file(self.filename)
            sample['Sample@TIME']=np.array(sample['Sample@TIME'])/1000   # [s]
            sample['Sample@POSN']=np.array(sample['Sample@POSN'])*1000   # [mm]
            sample['Sample@LOAD']=np.array(sample['Sample@LOAD'])*10     # [N]
            key_names = []
            for key in sample.keys():
                key_names.append(key)
            return sample, key_names
        else:
            print('{} is not found'.format(self.filename))
            return None
    
    def save_data(self, output_name = 'test.dat'):
        if self.outname is None:
            self.outname = output_name
            
        # Write header and data to the file
        with open(self.outname, 'w') as file:
            file.write(self.header + '\n')
            for row in zip(self.data['Sample@POSN'], self.data['Sample@LOAD']):
                file.write(','.join(map(str, row)) + '\n')
        print('Data has been written to {}'.format(self.outname))            
        
    def plot_data(self):
        if self.data is not None:
            plt.figure()
            plt.plot(self.data['Sample@POSN'], self.data['Sample@LOAD'], 'o')
            plt.xlabel('Position (mm)',fontsize=20)
            plt.ylabel('Load (N)',fontsize=20)
            plt.tick_params(axis='both',labelsize=16)
            plt.tight_layout()
            plt.show()
    
    def parse_xml_file(self, xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()

        data, sample = self.read_xml_data(root)
        return sample   

    def read_xml_data(self, element):
        data_values = {}
        sample_values = {}

        if element.text:
            data_values[element.tag] = element.text.strip()

        for child in element:
            child_data, child_sample_values = self.read_xml_data(child)

            if child.tag in data_values:
                if isinstance(data_values[child.tag], list):
                    data_values[child.tag].append(child_data)
                else:
                    data_values[child.tag] = [data_values[child.tag], child_data]
            else:
                data_values[child.tag] = child_data

            if child.tag == 'Sample':
                if child.attrib:
                    for attr_key, attr_value in child.attrib.items():
                        attr_tag = child.tag + "@" + attr_key
                        if attr_tag in sample_values:
                            sample_values[attr_tag].append(float(attr_value))
                        else:
                            sample_values[attr_tag] = [float(attr_value)]
            elif child_sample_values:
                sample_values.update(child_sample_values)

        return data_values, sample_values

