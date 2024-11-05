# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:51:26 2020

@author: Arsenic
"""
import numpy as np
from math import cos,sin, pi,atan2,asin

def degToRad(deg):
    rad = np.array(deg)*pi/180
    return rad

def radToDeg(rad):
    deg = np.array(rad)*180/pi
    return deg

def rot_matrix(eulers):
    "macro base in crystal frame"
    "RD = rotated_matrix[0,:]"
    "TD = rotated_matrix[1,:]"
    "ND = rotated_matrix[2,:]"
    eulers = degToRad(eulers)
    rotated_matrix = np.zeros((3,3))
    rotated_matrix[0,0] = cos(eulers[0])*cos(eulers[2])-sin(eulers[0])*sin(eulers[2])*cos(eulers[1])
    rotated_matrix[0,1] = -cos(eulers[0])*sin(eulers[2])-sin(eulers[0])*cos(eulers[2])*cos(eulers[1])
    rotated_matrix[0,2] = sin(eulers[0])*sin(eulers[1])
    rotated_matrix[1,0] = sin(eulers[0])*cos(eulers[2])+cos(eulers[0])*sin(eulers[2])*cos(eulers[1])
    rotated_matrix[1,1] = -sin(eulers[0])*sin(eulers[2])+cos(eulers[0])*cos(eulers[2])*cos(eulers[1])
    rotated_matrix[1,2] = -cos(eulers[0])*sin(eulers[1])
    rotated_matrix[2,0] = sin(eulers[2])*sin(eulers[1])
    rotated_matrix[2,1] = cos(eulers[2])*sin(eulers[1])
    rotated_matrix[2,2] = cos(eulers[1])
    return rotated_matrix

def gen_directions(lattice):
    "generates the possible slip directions"
    if lattice == 'BCC':
        directions = np.array([[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1]])
    if lattice == 'FCC':
        directions = np.array([[0,1,1], [0,-1,1], [1,0,1], [-1,0,1],[-1,1,0],[1,1,0]])
    return directions

def gen_planes(plane_type):
    if plane_type == '110':
        planes = np.array([[1,1,0],[1,0,1],[0,1,1],[-1,1,0],[-1,0,1],[0,-1,1]])
        return planes
    if plane_type == '112':
        planes = np.array([[1,1,2],[1,2,1],[2,1,1],[1,1,-2],[1,-2,1],[-2,1,1],[-1,1,2],[2,-1,1],[1,2,-1],[1,-1,2],[2,1,-1],[-1,2,1]])
        return planes
    if plane_type == '123':
        planes = np.array([[1,2,3],[3,1,2],[2,3,1],[-1,2,3],[3,-1,2],[2,3,-1],[1,-2,3],[3,1,-2],[2,-3,1],[1,2,-3],[-3,1,2],[2,-3,1],[3,2,1],[2,1,3],[1,3,2],[-3,2,1],[2,1,-3],[1,-3,2],[3,-2,1],[-2,1,3],[1,3,-2],[3,2,-1],[2,-1,3],[-1,3,2]])
        return planes
    if plane_type == '134':
        planes = np.array([[1,3,4],[4,1,3],[3,4,1],[-1,3,4],[4,-1,3],[3,4,-1],[1,-3,4],[4,1,-3],[-3,4,1],[1,3,-4],[-4,1,3],[3,-4,1],[4,3,1],[3,1,4],[1,4,3],[-4,3,1],[3,1,-4],[1,-4,3],[4,-3,1],[-3,1,4],[1,4,-3],[4,3,-1],[3,-1,4],[-1,4,3]])
        return planes
    if plane_type == '111':
        planes = np.array([[-1,1,1],[1,1,1],[-1,-1,1],[1,-1,1]])
        return planes
    else:
        print("unregistered plane type")

def calc_sfs(eulers,slipPlane,lattice):
    rotated_matrix = rot_matrix(eulers)
    planes = gen_planes(slipPlane)
    schmidFact = []
    dirs = gen_directions(lattice)
    for plane in range(planes.shape[0]):
        for direction in range(dirs.shape[0]):
            if np.dot(planes[plane],dirs[direction]) == 0:
                schmidFact.append((np.abs(np.multiply(np.dot(rotated_matrix[1,:],planes[plane]),np.dot(rotated_matrix[1,:],dirs[direction]))))/(np.multiply(np.linalg.norm(planes[plane]),np.linalg.norm(dirs[direction]))))
    return schmidFact
    
def plane_trace_components(plane,eulers):
    "calculates the trace angle of given plane in given crystal onto sample surface"
    "opening convention from left to right, convention '2' from ImageJ"
    matrix = rot_matrix(eulers)
    trace_x = np.dot(np.cross(plane, matrix[2,:]), matrix[1,:])/np.linalg.norm(np.cross(plane, matrix[2,:]))  # divide by norm of vector
    trace_y = -np.dot(np.cross(plane, matrix[2,:]), matrix[0,:])/np.linalg.norm(np.cross(plane, matrix[2,:]))  # divide by norm of vector
    return trace_x, trace_y


def project_direction(direction,eulers):
    "calculates the normed components of a given direction in sample reference frame"
    matrix = rot_matrix(eulers)
    b_x = np.dot(matrix[0,:],direction)
    b_y = np.dot(matrix[1,:],direction)
    return b_x,b_y

def get_gamma_angle(plane,direction,eulers):
    "calculates the gamma angle associated to a slip system"
    b_x, b_y = project_direction(direction,eulers)
    plane_trace_x, plane_trace_y = plane_trace_components(plane,eulers)
    longi = plane_trace_x*b_x + plane_trace_y*b_y
    transv = -(plane_trace_y)*b_x + (plane_trace_x)*b_y
    gamma = atan2(longi,transv)
    return gamma

def random_quats(n):
    "uniform quats population"
    u1,u2,u3 = np.random.rand(3,n)
    return np.array([np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
                     np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
                     np.sqrt(u1)     * np.sin(2 * np.pi * u3),
                     np.sqrt(u1)     * np.cos(2 * np.pi * u3)]).T

def quats_to_euler(quats):
    "conversion from quaternions to euler angles"
    "/!\ resulting euler angles in radians"
    numQuats = quats.shape[0]
    eulers = np.zeros((numQuats,3),dtype=float)
    for quat in range(numQuats):
        eulers[quat,0] = atan2(2*(quats[quat,0]*quats[quat,1]+quats[quat,2]*quats[quat,3]),1-2*(quats[quat,1]*quats[quat,1]+quats[quat,2]*quats[quat,2]))
        eulers[quat,1] = asin(2*(quats[quat,0]*quats[quat,2]-quats[quat,3]*quats[quat,1]))
        eulers[quat,2]= atan2(2*(quats[quat,0]*quats[quat,3]+quats[quat,1]*quats[quat,2]),1-2*(quats[quat,2]*quats[quat,2]+quats[quat,3]*quats[quat,3]))
    return eulers