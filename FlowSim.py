#Joshua Parker
#Program to simulate gas molecules attaching to a MOF at dynamic temperatures,
# and graphing the changing state.
#Thanks goes to Professor Rob Owen for giving me examples of how to use
# matplotlib to create the graphs and save them as files.
#Thanks goes to Professor Stephen Fitzgerald for the overall guidance and
# suggested structure and formulas used.
import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import scipy.constants as sci
 
import time
import random

#Grid dimensions
gridX = 40 
gridY = 40
#MOF dimensions
mofX = 20 
mofY = 20
#Numbers of gas molecules
numH2 = 60
numHD = 0
numD2 = 60 

temperature = 1 #Ranges from 1 to 201, should correspond to Kelvin units. NOT ZERO, something is divided by temperature
movements = 10 #This is how many steps aka frames there are.
#As of 05/14/2017, it takes ~2 seconds for 10 movements and ~184 seconds for 1000. Other parameters such as number of gases and mof sites affect it, but largest 
#(known)impact is that of movements.
beta = 200/movements #amount by which temperature increases with each movement.

#Formula for relating temperature and binding energy of a gas molecule to the probability of it escaping the MOF
def escapeRate(E,T):
	n = np.exp(-E/T)
	return n
#This should be an object which represents an arbitrary gas molecule.
#Each gas molecules needs to know what type of gas it is, its coordinates,
#and whether it's adsorbed to a MOF site.
class Gas:
	def __init__(self,type,a,b,dir):
		self.type = type
		self.adsorbed = False
		self.energy = 0
		self.x = a
		self.y = b
		self.oldX = a
		self.oldY = b
		self.direction = dir
		if type == "h2":
			self.energy = 10 #These are the bindingn energies of the different types of gases. 
		elif type == "hd":
			self.energy = 5 #These are the bindingn energies of the different types of gases. 
		elif type == "d2":
			self.energy = 30 #These are the bindingn energies of the different types of gases. 
		#1 is up, 2 is right, 3 is down, 4 is left
	def move(self,a,b):
		self.oldX = self.x
		self.oldY = self.y
		self.x = a
		self.y = b
	def reverse(self):
		self.direction=random.randint(1,4)
		
#Object which represents a MOF SITE
class Mof:
	def __init__(self,a,b):
		self.type = "open" # open means the site is unoccupied, closed means that it is occupied
		self.external = False #External refers to whether this site is on the exterior of the mof sample, aka whether it's adjacent to space to which gas can escape.
		self.x = a
		self.y = b
		self.neighborX = -10 #These are used as part of checking whether a site is external, and if so what the coordinates are of its adjacent free space.
		self.neighborY = -10 #-10 is just the default, indicating that it's not external, and is changed when the mof site is switched to external = true
	def fill(self):
		self.type = "closed"
	def clear(self):
		self.type = "open"
	def neighbor(self,nx,ny):
		self.external = True
		self.neighborX = nx
		self.neighborY = ny


#Passing in the list of mofs and list of gas molecules, specific coordinates, return a list of all objects with those coordinates.
def location(mofList,gasList,x,y):
	arr = []
	gas = False
	mof = False
	for m in mofList:
		if m.x == x and m.y == y:
			arr.append(m)
			if mof:
				print("two mofs on same space??")
			mof = True
	for g in gasList:
		if g.x == x and g.y == y:
			arr.append(g)
			if gas:
				print("two gases on same space??")
			gas = True

	return arr

#This moves all gas molecules at each time step. It's the mass of moving parts most likely to screw up when changing things, so be careful.
def step(mofList,gasList,t):
	T = temperature + beta*t
	for g in gasList:
		#Probability of Escape = e^(-E/kT)
		#Since it's just relative, we can do it as proportions
		#We want temperature to change with time.
		#desorptionRate = E/T where E = Gas-specific adsorptionRate

		#adsorptionRate is 1-desorptionRate
		if not g.adsorbed:
			#case 1: no other objects have coordinates above current, and is within grid.
			#case 2: intended coordinates out of bounds, instead have it stop and change direction to random direction.
			#case 3: intended coordinates are taken by another gas molecule, instead have it stop and change direction to random direction.
			#case 4: intended coordinates are open mof, and it gets adsorbed.
			#case 5: intended coordinates are open mof, and it does not get adsorbed, instead have it stop and change direction to random direction.
			#case 6: intended coordinates are closed mof, instead have it stop and change direction to random direction.
				#fringecase: it is surrounded by other objects, have it just stay as is.

			# move - if within boundaries, if len(loc)==0
			# reverse - if beyond boundaries, if loc()==mof and permit==0, if loc().type=="closed", if len(loc())==1 and type(loc()[0])==Gas
			# adsorb
			# stay
			permit = random.random() 
			if g.direction==1:
				loc = location(mofList,gasList,g.x,g.y-1)
				if len(loc) == 0 and g.y > 0: #move up if open space ahead
					g.y = g.y-1
				elif g.y == 0: #reverse direction if boundary ahead
					g.reverse()
				elif len(loc) == 1:	#if there is only 1 object in space ahead
					if type(loc[0]) == Gas: #If that object is a gas, reverse
						g.reverse()
					elif type(loc[0]) == Mof: #if that object is a mof:
						if loc[0].type == "closed" or permit <= escapeRate(g.energy,T): #if that mof is closed, or permit is denied, reverse
							g.reverse()
						elif loc[0].type == "open" and permit > escapeRate(g.energy,T): #if that mof is open, and permit is allowed, move in, fill mof, set gas to adsorbed=True
							g.y = g.y-1
							loc[0].fill()
							g.adsorbed = True
				elif len(loc) == 2:	#If the space is occupied by 2 things ( hopefully gas + mof ) then reverse
					g.reverse()
				else:
					print("wtf1")

			elif g.direction==2:
				loc = location(mofList,gasList,g.x+1,g.y)
				if len(loc)==0 and g.x < gridX-1: #move right
					g.x = g.x+1
				elif g.x == gridX-1:
					g.reverse()
				elif len(loc) == 1:
					if type(loc[0]) == Gas:
						g.reverse()
					elif type(loc[0]) == Mof:
						if loc[0].type == "closed" or permit <= escapeRate(g.energy,T):
							g.reverse()
						elif loc[0].type == "open" and permit > escapeRate(g.energy,T):
							g.x = g.x+1
							loc[0].fill()
							g.adsorbed = True
				elif len(loc) == 2:
					g.reverse()
				else:
					print("wtf2")

			elif g.direction==3:
				loc = location(mofList,gasList,g.x,g.y+1)
				if len(loc)==0 and g.y < gridY-1: #move down
					g.y = g.y+1
				elif g.y == gridY-1:
					g.reverse()
				elif len(loc) == 1:
					if type(loc[0]) == Gas:
						g.reverse()
					elif type(loc[0]) == Mof:
						if loc[0].type == "closed" or permit <= escapeRate(g.energy,T):
							g.reverse()
						elif loc[0].type == "open" and permit > escapeRate(g.energy,T):
							g.y = g.y+1
							loc[0].fill()
							g.adsorbed = True
				elif len(loc) == 2:
					g.reverse()
				else:
					print("wtf3")

			elif g.direction==4:
				loc = location(mofList,gasList,g.x-1,g.y)
				if len(loc)==0 and g.x > 0: #move left
					g.x = g.x-1
				elif g.x == 0:
					g.reverse()
				elif len(loc) == 1:
					if type(loc[0]) == Gas:
						g.reverse()
					elif type(loc[0]) == Mof:
						if loc[0].type == "closed" or permit <= escapeRate(g.energy,T):
							g.reverse()
						elif loc[0].type == "open" and permit > escapeRate(g.energy,T):
							g.x = g.x-1
							loc[0].fill()
							g.adsorbed = True
				elif len(loc) == 2:
					g.reverse()
				else:
					print("wtf4")

			else:
				print("Huh??")
		else: #If the gas molecule is adsorbed to a mof site:
			#move randomly
			desorbed = False
			loc = location(mofList,gasList,g.x,g.y)
			if loc[0].type != "closed": #This is merely an error checker to make sure that all gas molecules which are adsorbed are in fact on a mof site
				print("adsorbed gas is not on a closed mof?")
				print(type(loc[0]))
			if loc[0].external == True: #If the gas molecule is on an exterior site, it has a chance to escape.
				chance = random.random() 
				adjacencyX = loc[0].neighborX
				adjacencyY = loc[0].neighborY
				if adjacencyX  == -10 or adjacencyY == -10: #This is another error checker, making sure that the adjacency coordinates are not the defaults
				#Otherwise we could be sending gas molecules to oblivion if that wasn't working.
					print("adj wrong")
				elif len(location(mofList,gasList, adjacencyX, adjacencyY))==0: #This means that the adjacent spot is empty space
					if escapeRate(g.energy,T) > chance:
						#Desorption
						g.x = adjacencyX
						g.y = adjacencyY
						loc[0].type = "open"
						g.adsorbed = False
						desorbed = True
			#choose random direction
				# if the space in that direction is closed or out of bounds, reroll.
			if not desorbed:
				moved = False
				#while not moved:
				rng = random.randint(1,4)
				if rng == 1: #try to move up
					target = location(mofList,gasList,g.x,g.y+1)
					if len(target) == 1:
						if type(target[0]) == Mof and target[0].type == "open":
							loc[0].clear()
							g.y = g.y+1
							target[0].fill()
							moved = True

				elif rng == 2: #try to move right
					target = location(mofList,gasList,g.x+1,g.y)
					if len(target) == 1:
						if type(target[0]) == Mof and target[0].type == "open":
							loc[0].clear()
							g.x = g.x+1
							target[0].fill()
							moved = True

				elif rng == 3: #try to move down
					target = location(mofList,gasList,g.x,g.y-1)
					if len(target) == 1:
						if type(target[0]) == Mof and target[0].type == "open":
							loc[0].clear()
							g.y = g.y-1
							target[0].fill()
							moved = True

				elif rng == 4: #try to move left
					target = location(mofList,gasList,g.x-1,g.y)
					if len(target) == 1:
						if type(target[0]) == Mof and target[0].type == "open":
							loc[0].clear()
							g.x = g.x-1
							target[0].fill()
							moved = True

#Initial setup method, creates all of the mof sites. Pass in the list where you want to store the mof objects, and the dictionary for setup.
def setMof(mofList):
	for m in range((gridX-mofX)//2, (gridX+mofX)//2):
		for o in range((gridY-mofY)//2, (gridY+mofY)//2):
			mof = Mof(m,o)
			if m == (gridX-mofX)//2:
				mof.neighbor((gridX-mofX)//2-1,o)
			elif m == (gridX+mofX)//2-1:
				mof.neighbor((gridX+mofX)//2,o)
			elif o == (gridY-mofY)//2:
				mof.neighbor(m,(gridY-mofY)//2-1)
			elif o == (gridY+mofY)//2-1:
				mof.neighbor(m,(gridY+mofY)//2)
			mofList.append(mof)

#Initial setup method, creates all of the gas molecules. Makes sure that they all start in free space.
def seed(mofList,gasList):
	h2 = 0
	hd = 0
	d2 = 0
	while h2 < numH2:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 0:
			g = Gas("h2",x,y,d)
			gasList.append(g)
			h2+=1
	while hd < numHD:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 0:
			g = Gas("hd",x,y,d)
			gasList.append(g)
			hd+=1
	while d2 < numD2:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 0:
			g = Gas("d2",x,y,d)
			gasList.append(g)
			d2+=1
#Initial setup method, creates all of the gas molcules the same way as seed, but makes sure that they all start off on a mof site.
def seedWithin(mofList,gasList):
	h2 = 0
	hd = 0
	d2 = 0
	while h2 < numH2:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 1:
			if loc[0].type == "open":
				g = Gas("h2",x,y,d)
				loc[0].type = "closed"
				g.adsorbed = True
				gasList.append(g)
				h2+=1
	while hd < numHD:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 1:
			if loc[0].type == "open":
				g = Gas("hd",x,y,d)
				loc[0].type = "closed"
				g.adsorbed = True
				gasList.append(g)
				hd+=1
	while d2 < numD2:
		x = random.randint(0,gridX-1)
		y = random.randint(0,gridY-1)
		d = random.randint(1,4)
		loc = location(mofList,gasList,x,y)
		if len(loc) == 1:
			if loc[0].type == "open":
				g = Gas("d2",x,y,d)
				loc[0].type = "closed"
				g.adsorbed = True
				gasList.append(g)
				d2+=1	
#This method handles the visuals. It takes in the list of mofs and list of gases, as well as the time step, and plots it, then saves it in a sub directory of the
# current directory called "testframes". You'll definitely want to make a folder for the images alone. 
def draw(mofList,gasList,stepNumber):

#BEWARE: the latter portion of this method relies on a bunch of documentation for matplotlib, look it up and save past iterations before changing things, 
# or it could get messier than it already is.

	freeH2 = 0
	freeHD = 0
	freeD2 = 0
	trappedH2 = 0
	trappedHD = 0
	trappedD2 = 0
	#Coordinates of the open mof sites
	mX = []
	mY = []
	#Coordinates of H2 molecules floating around
	h2X = []
	h2Y = []
	#Coordinates of HD molecules floating around
	hdX = []
	hdY = []
	#Coordinates of D2 molecules floating around
	d2X = []
	d2Y = []
	#Coordinates of trapped H2 molecules
	th2X = []
	th2Y = []
	#Coordinates of trapped HD molecules
	thdX = []
	thdY = []
	#Coordinates of trapped D2 molecules
	td2X = []
	td2Y = []


	for g in gasList:
		if g.type=="h2" and g.adsorbed==False:
			freeH2 += 1
			h2X.append(g.x)
			h2Y.append(g.y)
		elif g.type=="h2":
			trappedH2 += 1
			th2X.append(g.x)
			th2Y.append(g.y)
		elif g.type=="hd" and g.adsorbed==False:
			freeHD += 1
			hdX.append(g.x)
			hdY.append(g.y)
		elif g.type=="hd":
			trappedHD += 1
			thdX.append(g.x)
			thdY.append(g.y)
		elif g.type=="d2" and g.adsorbed==False:
			freeD2 += 1
			d2X.append(g.x)
			d2Y.append(g.y)
		elif g.type=="d2":
			trappedD2 += 1
			td2X.append(g.x)
			td2Y.append(g.y)
		else:
			print("We hit an unknown gas??")
	for m in mofList:
		if m.type=="open":
			mX.append(m.x)
			mY.append(m.y)

	if freeH2 + trappedH2 != numH2:
		print("H2 is not conserved")
	if freeHD + trappedHD != numHD:
		print("HD is not conserved")
	if freeD2 + trappedD2 != numD2:
		print("D2 is not conserved")

	#******************** This is where it gets really messy ******************************

	#This following code sets the axis ranges, plots the scatter plot,
	# gives it axis labels, and then saves and shows it.
	fig,(ax,barplot) = plt.subplots(1,2,gridspec_kw = {'width_ratios':[5,1]})
	ax.set_xlim(0,gridX)
	ax.set_ylim(0,gridY)
	#plot the mof markers
	ax.scatter(mX,mY,marker='s',facecolors='none', edgecolors='g',label="Mof Site")
	ax.scatter(th2X,th2Y,marker='s',facecolors='none', edgecolors='g')
	ax.scatter(td2X,td2Y,marker='s',facecolors='none', edgecolors='g')
	blueStar, = ax.plot(th2X,th2Y,"b*")
	blackStar, = ax.plot(td2X,td2Y,"r*")
	#plot the H2 markers
	ax.scatter(h2X, h2Y, c='b', marker='o',label="H2")
	#plot the D2 markers
	ax.scatter(d2X, d2Y, c='r', marker='o',label="D2")
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(0.98, 0.5))


	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)

	#plt.annotate('Free H2 molecules: %s' %(freeH2),xy=(25,-1),annotation_clip=False)
	#plt.annotate('trapped H2 molecules: %s' %(trappedH2),xy=(25,-3),annotation_clip=False)
	#plt.annotate('Free D2 molecules: %s' %(freeD2),xy=(5,-1),annotation_clip=False)
	#plt.annotate('Trapped D2 molecules: %s' %(trappedD2),xy=(5,-3),annotation_clip=False)
	# We want a bar chart, where x axis is free H2, trapped H2, free D2, trapped D2, 
	# and y axis is numbers
	
	#I'm not sure why I need two of each, but this was the only way I could find to make them properly display. the barplot documentation is not as in-depth as I'd like.
	objects = ('H2','H2','D2','D2')
	
	plt.bar(1,trappedH2,color='b')
	plt.bar(2,trappedD2,color='r')
	barplot.set_ylim(0,100)
	barplot.set_xticklabels(objects,size=10)
	barplot.yaxis.set_visible(False)

	T = temperature + beta*stepNumber
	ax.annotate('Temperature: %s' %T,xy=(12,-2),annotation_clip=False)
	#This is the line where we choose where to save the file and what to call it.
	plt.savefig( "testFrames/SimStep%04i.png"%stepNumber)
	#plt.show() will display the graph every step, and the program will pause until you close it. You'll usually want it commented out/off.
	#plt.show()
	#plt.close(fig) lets things move a bit faster, look up documentation if you're curious, otherwise it's pretty irrelevant, just leave it be.
	plt.close(fig)
	

#What we need: we need a list of gas objects, a list of mof objects.
#We instantiate the objects, then step through them, and draw at end of each step.
#step() will need both lists of objects, and then have them move.
#draw() will need both lists of objects, plot the resulting image, and save it to a folder in the current directory.

#This is the main function which we call first
def main():
	startTime = time.process_time() #Because i like seeing how long it takes. Especially if you change things to take a long time, you can get an idea of how long to wait.
	mofList = [] #This will store all mof objects
	gasList = [] #This will store all gas objects

	setMof(mofList) #This initializes the MOF sample
	seedWithin(mofList,gasList) #This is interchangeable with seed (or should be), it creates the gas molecules.
	#The only difference is whether the gas molecules start inside or outside of the mof.
	draw(mofList,gasList,0) #Shows us the state right after initialization, before any movements.
	for t in range(1, movements+1): 
		step(mofList,gasList,t)
		draw(mofList,gasList,t)

	for g in gasList: #This is an error checker I found useful when making adjustments, checks to make sure that all of the gas molecules are within boundaries.
		if g.x >= gridX or g.x < 0:
			print("lost to the void")
		if g.y >= gridY or g.y < 0:
			print("lost to the abyss")
			print(g.y)

	endTime = time.process_time()
	print("Time elapsed: %s" %(endTime-startTime))
	
if __name__ == '__main__':
	main()