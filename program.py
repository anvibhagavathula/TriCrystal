
from helpers import * 

#--------------------------------------------------------------------------------------------------------

# WELCOME MESSAGE 
print ('\n\n*************************** WELCOME TO TriCRYSTAL ***************************\n\n * TriCRYSTAL: Commensurate and incommensurate crystal structures of layered materials\n (c) 2020 Johnson Chemistry Group, Dalhousie University\n')
print (' Description: TriCrystal builds commensurate and incommensurate structures of layered materials.\n Current version reads CIF files and writes the new structure to a QUANTUM ESPRESSO input file.\n Additional information such as the bond distance between atoms, lattice vectors in Bohr and Angstrom, and a simple 3D plot of each layer is also provided.\n\n Authored by: T. Kabengele \n Contact: tilas.kabengele@dal.ca \n\n using: Python 3.7.6 \n with libraries: numpy, matplotlib \n imported library: crystals \n (c) Copyright 2020, Laurent P. René de Cotret \n Install via: pip install crystals\n\n')
print ('* TriCRYSTAL--', datetime.now(),'\n\n')
print ('**************************************************************************\n \n')

# USER INPUTTED CRYSTAL CIF FILE 
my_crystal = Crystal.from_cif(input('***Input cif file*** \n'))

# USER INPUTTED SUPERCELL PARAMETERS AND ROTATION ANGLE 
print ('\n***Rotation parameters*** ')
m = int(input('Enter m '))
n = int(input('Enter n '))

# USER INPUTTED VACUUM CELL PARAMETER 
print('\n***Cell Vacuuum Parameter')
z = float(input('Enter z '))

#--------------------------------------------------------------------------------------------------------

dir1 = os.path.join("/Users/anvitabhagavathula/Desktop/TriCrystal/", 'periodic_table.csv')
colnames = ['number', 'symbol']
periodic_table = pandas.read_csv(dir1, usecols=colnames)
number = periodic_table.number.tolist()
symbol = periodic_table.symbol.tolist()

# CALCULATING LATTICE PARAMETERS 
a, b, c, alpha, beta, gamma = my_crystal.lattice_parameters

# INITIALISING BOTTOM AND MIDDLE LAYERS 
tt_mid,tt_bot,elt_mid,elt_bot = bulk(my_crystal)

# CALCULATING LATTICE VECTORS 
a1, a2, a3 = my_crystal.lattice_vectors
uc = a1,a2,a3
uc = np.array(uc)

#--------------------------------------------------------------------------------------------------------

print ('\n\nIntializing atoms...\n\n')

print ('Initial TOP atoms..')
for i in range(0,len(elt_mid)):
    s = np.array(tt_mid[i])
    s = np.dot(s, np.linalg.inv(uc))
    print ('Atom No.',i+1, ' ',elt_mid[i], ' ',s)

print ('\nInitial BOTTOM atoms..')
for j in range(0,len(elt_bot)):
    s = np.array(tt_bot[j])
    s = np.dot(s, np.linalg.inv(uc))
    print ('Atom No.',len(elt_mid)+j+1, ' ',elt_bot[j], ' ',s)

# SELECTING ZERO-TH ATOMS FROM MIDDLE AND BOTTOM LAYERS 
print ('\nSelect zeroeth MIDDLE atom')
zeroeth1 = int(input('Enter Atom No. '))
print ('\nSelect zeroeth BOTTOM atom')
zeroeth2 = int(input('Enter Atom No. '))

# CORRECTING ARRAY INDICES FOR ZEROETH ATOMS 
idx1 = zeroeth1-1
idx2 = (zeroeth2-1)-len(elt_mid)
print ('\nZeroeth TOP (angstrom)', elt_mid[idx1], tt_mid[idx1])
print ('\nZeroeth BOTTOM (angstrom)', elt_bot[idx2], tt_bot[idx2])

# FINDING THE BOND LENGTH 
lengths = pdist(tt_bot, 'euclidean')
bond_distance = round(min(lengths),3)
print ('\nBond distance = ', bond_distance)

# PRINTING LATTICE PARAMETERS 
print ('\nLattice Vectors (Angstrom)')
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a1[0],a1[1],a1[2]))
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a2[0],a2[1],a2[2]))
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a3[0],a3[1],a3[2]))

a2b = 1.8897259886
print ('\nLattice Vectors (Bohr)')
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a1[0]*a2b,a1[1]*a2b,a1[2]*a2b))
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a2[0]*a2b,a2[1]*a2b,a2[2]*a2b))
print (' ','{:12.6f} {:12.6f} {:12.6f}'.format(a3[0]*a2b,a3[1]*a2b,a3[2]*a2b))

# CALCULATING ROTATION ANGLE 
A = np.array([1,0,0])
B = np.array([np.cos(np.deg2rad(60)),np.sin(np.deg2rad(60)),0])
V = m*A + n*B
rotation_angle = np.arccos((np.dot(A,V.T))/(np.linalg.norm(A)*np.linalg.norm(V)))
rotation_angle = np.rad2deg(rotation_angle)

# CALCULATING ROTATION MATRIX 
theta = np.deg2rad(rotation_angle)
phi = np.deg2rad(60) - 2*theta
R = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
print ('\nRotation angle theta (degrees) = ', rotation_angle)
print ('\nMoire angle gamma (degrees) = ',np.rad2deg(phi))

## CALCULATING INTERLAYER SCALING
inter_scale = scale_interlayer(z)

#--------------------------------------------------------------------------------------------------------

print ("\n\nCALCULATING ATOMIC POSITIONS...")
for i in range(1,2):
        print ('\n\nPlease wait...\n\n')
        sys.stdout.flush()

print ("&control")
print ( " title='crystal',")
print ( " prefix='crystal',")
print ( " pseudo_dir='/users/abhagava/scratch/pseudos/B86bPBE/',")
print ( " calculation='scf',\n etot_conv_thr=1.0D-5,\n forc_conv_thr=1.0D-4,\n/")
print ( "&system")
print (" ibrav=0,")
print (" nat=42,")

nat = ntype(my_crystal)
print (" ntyp= ",nat,",",sep="")
print (" ecutwfc=90.0,")
print (" ecutrho=900.0,")
print (" vdw_corr=\'xdm\',\n/")
print ("&electrons")
print (" conv_thr = 1d-8\n/\n&ions\n/\n&cell\n/")
print ("ATOMIC_SPECIES")

for atm in my_crystal.chemical_composition:
        ELE = str(atm)
        ele = ELE.lower()
        mass_number = Element(ele).mass
        print (ele, "     ",mass_number," ",ele,'.UPF',sep="")


print('\nATOMIC_POSITIONS crystal')

#--------------------------------------------------------------------------------------------------------

# LOOPS FOR BOTTOM LAYER 

bt0 = time.time()
bt1 = time.time()


tt1 = []
tt2 = []
atoms_bot = []
symb1 = []
symb2 = []
k = 0
for atm in tt_bot:
    u = elt_bot[k]
    k = k + 1
    for i in range(0,(m+n)*15):
        ttx = atm + i*a1
        tt1.append(ttx)
        symb1.append(u)
k = 0
for atm in tt1:
    u = symb1[k]
    k = k + 1
    for i in range(1,(n+m)*15):
        tty = atm + i*a2
        tt2.append(tty)
        symb2.append(u)
symb_bot = symb1 + symb2
atoms_bot = list(tt1) + list (tt2)

elbt1 = time.time() - bt1

### Initializing new unit cell ###

# new cell parameters
newa1b, newa2b, v1b, v2b, v3b, v4b = newcell(my_crystal,tt_bot[idx2],m,n)
unitcell = newa1b,newa2b,a3
unitcell = np.array(unitcell)

# polygon boundary
boundary_bot,p1b,p2b,p3b,p4b = poly(v1b,v2b,v3b,v4b)

# center atom
org = central(boundary_bot)
destinations = MultiPoint(atoms_bot)
nearest_geoms = nearest_points(org, destinations)
origin = np.array([nearest_geoms[1].x, nearest_geoms[1].y, 0])

ex = 1
Rb = np.array([[np.cos(0), -np.sin(0), 0], [np.sin(0), np.cos(0), 0], [0, 0, 1]])
v1r,v2r,v3r,v4r = rotcell(v1b*ex,v2b*ex,v3b*ex,v4b*ex,origin,Rb)
boundary_bot,p1b,p2b,p3b,p4b = poly(v1r,v2r,v3r,v4r)

# list of atomic numbers from symbols
symb_num_bot = []
for i in range(0,len(symb_bot)):
    idx = symbol.index(symb_bot[i])
    symb_num_bot.append(number[idx])

# number of types of atomic species
typ = ntype(my_crystal)

# Initializing check for which atoms lie within the new unit cell
bot = []
symbot = []
atoms_bot = atoms_bot - origin
supx,supy,supz = atoms_bot.T

for i in range(0,len(atoms_bot)):
    num = symb_num_bot[i]
    if inpoly(atoms_bot[i],boundary_bot) == True:
        bt =  supx[i], supy[i], supz[i]
        bot.append(bt)
        symbot.append(symb_bot[i])

botl = bot
bot = np.array(bot)
bot_frac = np.dot(bot, (np.linalg.inv(unitcell)))

for i in range(1,m+n):
    bot_frac[bot_frac<0] += 1
bot_frac[bot_frac>1] += -1

i = 0
sim = []
for atm1 in bot_frac:
    count = 0
    for atm2 in bot_frac:
        if (round(atm1[0],2) == round(atm2[0],2) and round(atm1[1],2) == round(atm2[1],2) and round(atm1[2],2)== round(atm2[2],2)) == True:
            count+=1
        if count > 1:
            sim.append(i)
    i+=1
sim = np.array(sim)
sim = np.unique(sim)

if len(sim) >= 1:
    bot_frac = np.delete(bot_frac, sim[1:], 0)
    symbot = np.delete(symbot, sim[1:], 0)

elbt = time.time() - bt0

#--------------------------------------------------------------------------------------------------------

# LOOPS FOR MIDDLE LAYER 
md0 = time.time()

tt1 = []
tt2 = []
atoms_mid = []
symb1 = []
symb2 = []
k = 0
for atm in tt_mid:
    u = elt_mid[k]
    k = k + 1
    for i in range(0,(m+n)*15):
        ttx = atm + i*a1
        tt1.append(ttx)
        symb1.append(u)
k = 0
for atm in tt1:
    u = symb1[k]
    k = k + 1
    for i in range(1,(n+m)*15):
        tty = atm + i*a2
        tt2.append(tty)
        symb2.append(u)
symb_mid = symb1 + symb2
atoms_mid = list(tt1) + list (tt2)

### Initializing new unit cell ###

# new cell parameters
newa1t, newa2t, v1t, v2t, v3t, v4t = newcell(my_crystal,tt_mid[idx1],m,n)
unitcell2 = newa1t,newa2t,a3
unitcell2 = np.array(unitcell2)

# polygon boundary
boundary_mid,p1t,p2t,p3t,p4t = poly(v1t,v2t,v3t,v4t)

# center atom
org = central(boundary_mid)
destinations = MultiPoint(atoms_mid)
nearest_geoms = nearest_points(org, destinations)
origin = np.array([nearest_geoms[1].x, nearest_geoms[1].y, 0])

ex = 1
Rm = Rb
v1r,v2r,v3r,v4r = rotcell(v1t*ex,v2t*ex,v3t*ex,v4t*ex,origin,Rm)
boundary_mid,p1t,p2t,p3t,p4t = poly(v1r,v2r,v3r,v4r)

# list of atomic numbers from symbols
symb_num_mid = []
for i in range(0,len(symb_mid)):
    idx = symbol.index(symb_mid[i])
    symb_num_mid.append(number[idx])

# number of types of atomic species
typ = ntype(my_crystal)

# Initializing check for which atoms lie within the new unit cell
mid = []
symmid = []
atoms_mid = np.dot((np.array(atoms_mid)-origin), R)
supx,supy,supz = atoms_mid.T


for i in range(0,len(atoms_mid)):
    num = symb_num_mid[i]
    if inpoly(atoms_mid[i],boundary_mid) == True:
        tt =  supx[i], supy[i], supz[i]
        mid.append(tt)
        symmid.append(symb_mid[i])

midl = mid
mid = np.array(mid)
mid_frac = np.dot(mid, (np.linalg.inv(unitcell2)))

# removing negative frac coordinates and those greater than 1
for i in range(1,m+n):
    mid_frac[mid_frac<0] += 1
mid_frac[mid_frac>1] += -1

i = 0
sim = []
for atm1 in mid_frac:
    count = 0
    for atm2 in mid_frac:
        if (round(atm1[0],2) == round(atm2[0],2) and round(atm1[1],2) == round(atm2[1],2) and round(atm1[2],2)== round(atm2[2],2)) == True:
            count+=1
        if count > 1:
            sim.append(i)
    i+=1
sim = np.array(sim)
sim = np.unique(sim)

if len(sim) >= 1:
    mid_frac = np.delete(mid_frac, sim[1:], 0)
    symtop = np.delete(symmid, sim[1:], 0)

elmd = time.time() - md0

#--------------------------------------------------------------------------------------------------------

# WRITING REMAINING SCF.IN FILE 

z_posn = [] 

for atm in mid_frac: 
    z_posn.append(atm[2])

i = 0
nat_bot=0
for atm in bot_frac:
    print ('{:2} {:12.6f} {:12.6f} {:12.6f}'.format(symbot[i], atm[0], atm[1], z_posn[i]-inter_scale))
    i+=1
    nat_bot+=1

i = 0
nat_mid=0
for atm in mid_frac:
    print ('{:2} {:12.6f} {:12.6f} {:12.6f}'.format(symmid[i], atm[0], atm[1], z_posn[i]))
    i+=1
    nat_mid+=1

i = 0
nat_top=0
for atm in bot_frac:
    print ('{:2} {:12.6f} {:12.6f} {:12.6f}'.format(symbot[i], atm[0], atm[1], z_posn[i]+inter_scale)) 
    i+=1
    nat_top+=1

print ('\nK_POINTS automatic')
print ('6 6 1 1 1 1')

old_a = np.array([a*0.5*np.sqrt(3), -0.5*b, 0])
old_b = np.array([0, b, 0])
old_c = np.array([0, 0, c])

newa1 = -m*old_a + n*old_b
newa2 = -m*old_b + n*old_a

uc1 = np.around(newa1b*1.8897259886, decimals=12)
uc2 = np.around(newa2b*1.8897259886, decimals=12)
uc3 = np.around(a3*1.8897259886, decimals=12)
uc1 = list(uc1)
uc2 = list(uc2)
uc3 = list(uc3)

print ('\nCELL_PARAMETERS bohr')
print ('   ','{:17.12f} {:17.12f} {:17.12f}'.format(uc1[0],uc1[1],uc1[2]))
print ('   ','{:17.12f} {:17.12f} {:17.12f}'.format(uc2[0],uc2[1],uc2[2]))
print ('   ','{:17.12f} {:17.12f} {:17.12f}'.format(uc3[0],uc3[1],z))