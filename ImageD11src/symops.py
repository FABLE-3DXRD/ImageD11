
from __future__ import print_function



# Systematic absence checking
# by Gavin Vaughan, April 2011

# Cell centering

def lattice_centre(h,k,l,ctype):
	if ctype == "P":
		return False
	elif ctype == "A":
		return (k+l)%2 != 0
	elif ctype == "B":
		return (h+l)%2 != 0
	elif ctype == "C":
		return (h+k)%2 != 0
	elif ctype == "I":
		return (h+k+l)%2 != 0
	elif ctype == "F":
		return (h+k)%2!=0 or (h+l)%2!=0 or (k+l)%2!=0
	elif ctype == "R":
		return (-h+k+l)%3 != 0
	else:
		return False

# rotations

def rotation_axis(h,k,l,rtype,axis):
	return False

# mirrors

def mirror_plane(h,k,l,axis):
	return False

# 1st type of screw axes

def screw_axis(h,k,l,stype,axis):
	if stype== '21' or stype == '42' or stype == '63':
		mod=2
	elif stype == '31' or stype == '32' or stype == '62' or stype == '64':
		mod=3
	elif stype == '41' or stype == '43':
		mod=4
	elif stype == '61' or stype == '65':
		mod=6
	else:
		raise Exception(stype+" is not a valid screw axis\n")   
	if axis == 1:
		if k != 0 and l != 0: return False
		return h%mod != 0
	if axis == 2:
		if h != 0 and l != 0: return False
		return k%mod != 0
	if axis == 3:
		if h != 0 and k != 0: return False
		return l%mod != 0
	return False
	

# Glide Planes 

def glide_plane(h,k,l,gtype,axis):
	if axis == 1:
		if h != 0: return False
	elif axis == 2:
		if k != 0: return False
	elif axis == 3:
		if l != 0: return False
	return glidetest(gtype,axis)
		
		
def glidetest(gtype,axis):
	if gtype == 'a':
		return h%2 != 0	
	if gtype == 'b':
		return k%2 != 0			
	if gtype == 'c':
		return l%2 != 0
	if gtype == 'n' :
		if axis == 1: return (k+l)%2 != 0
		if axis == 2: return (h+l)%2 != 0
		if axis == 3: return (h+k)%2 != 0
	if gtype == 'd':
		if axis == 1: return (k+l)%4 != 0
		if axis == 2: return (h+l)%4 != 0
		if axis == 3: return (h+k)%4 != 0
	return False
		
def test_absence(h, k, l, sg):
# break down the SG name into its symmetry components
# e.g. P21c -> P 1 21/c 1 
# To be unambiguous, it would be best if they were entered in the explicit way
# P21c could mean e.g. P 2 1 c or P 1 21 c or P 21 c 1 etc etc
# So use an unambiguous string like "P 1 21/c 1"
	symmop=sgstring.split()
	if len(symmop) != 4:		
		print("Only read %d symmops"%symmop)	
		raise Exception("You must supply at least 4 space separated symmetry operations"\
			 " and optional compound operations separated by '/'")
	
	is_absent=lattice_centre(h,k,l,symmop[0])
	if is_absent: return True
	
	for i in range(1,4):	
		if "/" in symmop[i]:
			op1, op2 = symmop[i].split("/")
			is_absent = checkop(h,k,l,op1,i)
			if is_absent: return True
			is_absent = checkop(h,k,l,op2,i)
			if is_absent: return True
		else:
			is_absent = checkop(h,k,l,symmop[i],i)
			if is_absent: return True
		
	return False
			
		
def checkop(h,k,l,op,axis):
# could do something clever, but this is just brute force
	rots=['1', '-1', '2', '-2', '3', '-3', '4', '-4', '6', '-6']
	mirror=['m']
	screws=['21', '31', '41', '61', '32', '42', '62', '43', '63', '64', '65']
	glides=['a', 'b', 'c', 'n', 'd']	
	
	if op in rots:
		return False
	if op in mirror:
		return False
	if op in screws:
		return screw_axis(h,k,l,op,axis)
	if op in glides:
		return glide_plane(h,k,l,op,axis)
	raise Exception(op+" is not a valid symmetry operation")
	

#never mind inversion for the absence check

if __name__ == "__main__":  
	import sys
	if len(sys.argv) != 8:
		print("Usage %s h k l sg_with_spaces"%(sys.argv[0]))
		print(len(sys.argv))
		sys.exit()
	try:
		h=int(sys.argv[1])
	except:
		print("Sorry %s is not an integer\n" % sys.argv[1])
		sys.exit()
	try:
		k=int(sys.argv[2])
	except:
		print("Sorry %s is not an integer\n" % sys.argv[2])
		sys.exit()
	try:
		l=int(sys.argv[3])
	except:
		print("Sorry %s is not an integer\n" % sys.argv[3])
		sys.exit()
		
	
	sgstring=sys.argv[4]+" "+sys.argv[5]+" "+sys.argv[6]+" "+sys.argv[7]
	absent=test_absence(h,k,l,sgstring)
	if absent is True:
		print("The ",h,k,l,"reflection is absent in space group "+sgstring)
	else:
		print("The ",h,k,l,"reflection is present in space group "+sgstring)
		
		
