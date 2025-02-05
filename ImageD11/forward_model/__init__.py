"""
forward model for box-beam 3DXRD and s3DXRD using point-focused beam
The aim is to exploit forward model for the following purposes:
1) facilitate analysis for (s)3DXRD by filtering and cleaning peaks for indexed grains, calculating completeness
2) bridging the combined analysis between DCT and (s)3DXRD for the same sample but measured in different geometries/instruments
3) forward model based reconstruction for both near-field and far-field geometries
4) parameters conversion among ImageD11 par, DCT parameters, poni, which all can be converted to a dictionary and compatible with fwd-DCT parameters
5) forward projector and back projector, aiming to improving/refining the grain and strain maps (under development)
6) other auxillary tools for multi-modality analysis, including interacting with PCT, sample environment, orientation conversion, deformation analysis (schmidt factor calc) etc.
Developed since Oct 2023

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

__version__ = "1.0.1"
__author__ = 'Haixing Fang',
__author_email__ = 'haixing.fang@esrf.fr'