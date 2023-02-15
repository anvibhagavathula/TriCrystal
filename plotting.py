from program import * 

#############################################################
##################### PLOTTING RESULTS ######################

pt1 = time.time()

## figure settings ##
fig = plt.figure(figsize = (8, 8))
ax = fig.add_subplot(111, projection='3d')
fig.set_facecolor('white')
ax.set_facecolor('white')
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
plt.axis('off')
plt.grid(b=None)


# Initializing bottom coordinates for new unitcell plotting
coord = [p1b, p2b, p3b, p4b]
coord.append(coord[0])
xb, yb = zip(*coord)

# plotting boundary at lowest z position
zb = lowest(my_crystal)
plt.plot(xb,yb,zb)

# plotting bottom layer atoms
bot = bot.T
supx,supy,supz = list(bot)
ax.scatter(supx, supy, supz, s=5*num, c='tab:green', alpha=0.4)
ax.scatter(supx, supy, supz, c='tab:blue')

# Adding bonds to bottom layer
neigh = NearestNeighbors(radius=bond_distance+0.1).fit(botl)
for atm in botl:
    rng = neigh.radius_neighbors([list(atm)])
    dis = np.asarray(rng[0][0])
    idx = np.asarray(rng[1][0])
    if inpoly(atm,boundary_bot) == True:
        for i in idx:
            atz = []
            atz.append(atm)
            atz.append(list(botl[i]))
            xb, yb, zb = zip(*atz)
            num = symb_num_bot[i]
            plt.plot(xb,yb,zb, '-bo',alpha=1)


# Initializing TOP coordinates for new unitcell plotting
coord = [p1t, p2t, p3t, p4t]
coord.append(coord[0])
xm, ym = zip(*coord)

# plotting the boundary at hiest z position
zm = highest(my_crystal)
plt.plot(xm,ym,zm)

# plotting TOP layer atoms
mid = mid.T
supx,supy,supz = list(mid)
ax.scatter(supx, supy, supz, s=5*num, c='tab:red', alpha=0.4)
ax.scatter(supx, supy, supz, c='tab:blue')

# Adding bonds to TOP layer
neigh = NearestNeighbors(radius=bond_distance+0.1).fit(midl)
for atm in midl:
    rng = neigh.radius_neighbors([list(atm)])
    dis = np.asarray(rng[0][0])
    idx = np.asarray(rng[1][0])
    if inpoly(atm,boundary_mid) == True:
        for i in idx:
            atz = []
            atz.append(atm)
            atz.append(list(midl[i]))
            xt, yt, zt = zip(*atz)
            num = symb_num_mid[i]
            plt.plot(xt,yt,zt, '-ro',alpha=1)

# TOP COORDINATES 
# Initializing bottom coordinates for new unitcell plotting
coord = [p1b, p2b, p3b, p4b]
coord.append(coord[0])
xb, yb = zip(*coord)

# plotting boundary at lowest z position
zb = lowest(my_crystal)+1.0
plt.plot(xb,yb,zb)

# plotting bottom layer atoms
bot = bot.T
supx,supy,supz = list(bot)
ax.scatter(supx, supy, supz, s=5*num, c='tab:green', alpha=0.4)
ax.scatter(supx, supy, supz, c='tab:blue')

# Adding bonds to bottom layer
neigh = NearestNeighbors(radius=bond_distance+0.1).fit(botl)
for atm in botl:
    rng = neigh.radius_neighbors([list(atm)])
    dis = np.asarray(rng[0][0])
    idx = np.asarray(rng[1][0])
    if inpoly(atm,boundary_bot) == True:
        for i in idx:
            atz = []
            atz.append(atm)
            atz.append(list(botl[i]))
            xb, yb, zb = zip(*atz)
            num = symb_num_bot[i]
            plt.plot(xb,yb,zb, '-bo',alpha=1)

plt.show()

elplt = time.time() - pt1


j2 = 90 + np.rad2deg(phi/2)
mk = 1/(2*(np.abs(np.sin((phi/2)))))


########################## SUMMARY REPORT ###############################

print ("\n********************* SUMMARY REPORT ***********************")
#print ('\nRotation angle (deg) = ', np.round(rotation_angle,3))
#print ('Relative Rotation (deg) = ',np.round(np.rad2deg(phi),3))
print ('\nHermann moire rotation = ', j2)
print ('Hermann moire constant = ', mk)
print ('\nTop atoms(rotated) = ',len(top_frac))
print ('Bottom atoms  = ',len(bot_frac))
print ('\nTotal atoms \n=', len(bot_frac)+len(top_frac))
#print ('\n Gamma = ', j2 + np.round(rotation_angle,3))
#print ( '\n lattice vectors = ',1.8897259886*a, 1.8897259886*b, 1.8897259886*c)
#print ('\n Erin method lattice vectors = ',1.8897259886*old_a,1.8897259886*old_b)
#print ('time for replication', elbt1)
#print ('time for plotting', elplt)
#print ('\ntotal time top layer',eltp)
#print ('total time bottom layer',elbt)
print ('\n*************************** Done!! **************************\n')