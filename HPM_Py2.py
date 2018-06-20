# This program will calculate the proper motions of stars and calculate their   #
# average PM to detect outliers which will be HPM candidates. The program will  #
# read in data from files obtained from the WSA to calculate the proper motion  #
# of stars in a region  using the coordinates from region data from ds9 and the #
# WSA database. Data will be output to a file.                                  #

#modules
import csv
from itertools import islice
import re
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import Counter

#it just werks

def restart(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate):

	print "\n\n\nEnd of program:"
	print "\n\t1. restart program\n\t2. redraw graph\n\t3. exit the program"
	input = int(raw_input("> "))

	if input == 1:
		start() #restarts the whole program

	elif input == 2:
		proper_motions(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate)
		#goes back to the proper motion calculations

	elif input == 3:
		exit() #exits the program

	else:
		print "Did not understand input. Please use a valid integer number as the input."
		restart(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate)
		#if you messed up the input then you'll go back to the start of this function to try again


def proper_motions(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate):

	print "\n####################################################################"
	print "#                         Proper Motions                           #"
	print "####################################################################"

	#TRAILING ~ INDICATES IMPORTANT BITS

	RA_outlierU = []; pm_RA_component2 = []; pm_mag_clipped = []
	DEC_outlierU = []; pm_DEC_component2 = []; pm_RA_clipped = []
	RA_outlierF = []; pm_mag2 = []; pm_DEC_clipped = []
	DEC_outlierF = []; JminusK2 = []
	pm_RA_outlier = []; HminusK2 = []
	pm_RA_outlier_error = []; JminusH2 = []
	pm_DEC_outlier = []; JminusK_error2 = []
	pm_DEC_outlier_error = []; JminusH_error2 = []
	pm_mag_outlier = []; HminusK_error2 = []
	pm_mag_outlier_error = []; RA_U = []
	JminusK_outlier = []; RA_F = []
	JminusH_outlier = []; DEC_U = []
	HminusK_outlier = []; DEC_F = []
	JminusK_outlier_error = []
	JminusH_outlier_error = []
	HminusK_outlier_error = []

	#converts final lists to numpy arrays to make calculations in one step rather than a loop
	RA1_matched_array = np.asfarray(RA1_matched3); DEC1_matched_array = np.asfarray(DEC1_matched3)
	RA2_matched_array = np.asfarray(RA2_matched3); DEC2_matched_array = np.asfarray(DEC2_matched3)

	d_RA = RA1_matched_array - RA2_matched_array #difference between U09B8 and Froebrich RA -- units of degrees
	d_DEC = DEC1_matched_array - DEC2_matched_array #difference between U09B8 and Froebrich DEC
	d_time = epoch1 - epoch2 #time difference between U09B8 and Froebrich epochs -- units of years

	pm_RA = ((d_RA * 3600) / d_time)  #finally calculating the proper motion components!
	pm_DEC = ((d_DEC * 3600)/ d_time) #units should be arcsec/yr

	pm_RA_mean = np.median(pm_RA) #THIS IS USED AS THE OFFSET OF THE COORDINATE IMAGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	pm_DEC_mean = np.median(pm_DEC)

	DEC_rads = np.radians(DEC2_matched_array)
	correction = (np.cos(DEC_rads)) #CALCULATES CORRECTION FACTOR~~~~~~~~~~~~~~~~~~~~~~~~
	pm_RA_component = correction*(pm_RA - pm_RA_mean) * 1000 #in units of mas/yr
	pm_DEC_component = (pm_DEC - pm_DEC_mean) * 1000
	pm_mag = np.sqrt(pm_RA_component**2 + pm_DEC_component**2) #units is mas/yr

	pm_mag_std = np.std(pm_mag) #DEFINING VALUES TO REMOVE ANYTHING BIGGER THAN 1 STANDARD DEVIATION
	#THE OFFSET HAS REMOVED THE MEAN SO THIS SHOULDN'T BE TAKEN INTO ACCOUNT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	minus_sigma = 0 - pm_mag_std
	plus_sigma = 0 + pm_mag_std

	x = len(pm_mag)
	for i in range(0, x):

		if minus_sigma <= pm_mag[i] <= plus_sigma: #THIS REMOVES THE OBJECTS THAT ARE >1 OR <1 SIGMA FROM THE MEAN~~~~~~~~~~~~~~~~~~~~
			pm_mag_clipped.append(pm_mag[i])
			pm_RA_clipped.append(pm_RA_component[i])
			pm_DEC_clipped.append(pm_DEC_component[i])

	pm_mag_clipped_array = np.asfarray(pm_mag_clipped); pm_mag_mean = np.median(pm_mag_clipped_array)
	pm_RA_clipped_array = np.asfarray(pm_RA_clipped); pm_RA_mean = np.median(pm_RA_clipped_array)
	pm_DEC_clipped_array = np.asfarray(pm_DEC_clipped); pm_DEC_mean = np.median(pm_DEC_clipped_array)
	#THESE ARE THE median? FOR THE RESPECTIVE QUANTATIES~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	pm_mag_std = np.std(pm_mag_clipped_array) #CALCULATES THE STD FOR THE CLIPPED ARRAY. THIS WILL BECOME THE NEW ERROR FOR THE RESPECTIVE QUANTATIES~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	pm_RA_std = np.std(pm_RA_clipped_array)
	pm_DEC_std = np.std(pm_DEC_clipped_array)

	JAperMag2_matched_array = np.asfarray(JAperMag2_matched3); J_error_matched_array = np.asfarray(J_error_matched3)
	K1AperMag2_matched_array = np.asfarray(K1AperMag2_matched3); K_error_matched_array = np.asfarray(K_error_matched3)
	HAperMag2_matched_array = np.asfarray(HAperMag2_matched3); H_error_matched_array = np.asfarray(H_error_matched3)

	JminusK = JAperMag2_matched_array - K1AperMag2_matched_array #CALCULATES THE COLOUR INDICES OF THE OBJECTS~~~~~~~~~~~~
	JminusH = JAperMag2_matched_array - HAperMag2_matched_array
	HminusK = HAperMag2_matched_array - K1AperMag2_matched_array

	JminusK_error = np.sqrt((J_error_matched_array)**2 + (K_error_matched_array)**2) #CAUCLATES THE ERROR IN QUADRATURE IN THE COLOUR INDICES~~~~~~~~~~~
	JminusH_error = np.sqrt((J_error_matched_array)**2 + (H_error_matched_array)**2)
	HminusK_error = np.sqrt((H_error_matched_array)**2 + (K_error_matched_array)**2)

	JminusK = JminusK.tolist(); JminusK_error = JminusK_error.tolist() #turns this stuff into lists
	JminusH = JminusH.tolist(); JminusH_error = JminusH_error.tolist()
	HminusK = HminusK.tolist(); HminusK_error = HminusK_error.tolist()

	#DEFINE LINE Y=MX+C ON A GRAPH. X AXIS IS H-K AND Y AXIS IS J-H
	red_length = len(JminusK)
	for q in range(0, red_length): #THIS LOOP REMOVES COLOURS BELOW OR BELOW THE REDDENING LINE AS WELL AS ANYTHING THAT HAS AN OBSCENE PM~~~~~~~~~~~~~~~~~~~~~~~~~~

		if 0.0 <= HminusK[q] <= 0.7:
			reddening_line1 = 1.55 * HminusK[q] + 0.7; reddening_line2 = 1.55 * HminusK[q] + 0

			if pm_mag[q] < 200:

				if JminusH[q] > reddening_line2:

						if JminusH[q] < reddening_line1:
							index = JminusH.index(JminusH[q])
							RA_U.append(RA1_matched3[index])
							RA_F.append(RA2_matched3[index])
							DEC_U.append(DEC1_matched3[index])
							DEC_F.append(DEC2_matched3[index])
							pm_RA_component2.append(pm_RA_component[index])
							pm_DEC_component2.append(pm_DEC_component[index])
							pm_mag2.append(pm_mag[index])
							JminusK2.append(JminusK[index])
							JminusH2.append(JminusH[index])
							HminusK2.append(HminusK[index])
							JminusK_error2.append(JminusK_error[index])
							JminusH_error2.append(JminusH_error[index])
							HminusK_error2.append(HminusK_error[index])

	pm_mag2_std = np.std(pm_mag2) #TAKES STANDARD DEVITION OF NEW SAMPLE
	sigma_1 = plt.Circle((0,0), (pm_mag2_std), color = 'r', fill=False) #defines the circles to be over plotted for SD
	sigma_2 = plt.Circle((0,0), (2*pm_mag2_std), color = 'r', fill=False)
	sigma_3 = plt.Circle((0,0), (3*pm_mag2_std), color = 'r', fill=False)
	sigma_4 = plt.Circle((0,0), (4*pm_mag2_std), color = 'r', fill=False)
	sigma_5 = plt.Circle((0,0), (5*pm_mag2_std), color = 'r', fill=False)

	pm_mag_length = len(pm_mag2)
	pm_mag = pm_mag.tolist()
	sigma = 3*pm_mag2_std

	for z in range(0, pm_mag_length):

		if  pm_mag2[z] > sigma:
			index = pm_mag2.index(pm_mag2[z])
			RA_outlierU.append(RA_U[index])
			DEC_outlierU.append(DEC_U[index])
			RA_outlierF.append(RA_F[index])
			DEC_outlierF.append(DEC_U[index])
			pm_RA_outlier.append(pm_RA_component2[index])
			pm_RA_outlier_error.append(pm_RA_std)
			pm_DEC_outlier.append(pm_DEC_component2[index])
			pm_DEC_outlier_error.append(pm_DEC_std)
			pm_mag_outlier.append(pm_mag2[index])
			pm_mag_outlier_error.append(pm_mag_std)
			JminusK_outlier.append(JminusK2[index])
			JminusK_outlier_error.append(JminusK_error2[index])
			JminusH_outlier.append(JminusH2[index])
			JminusH_outlier_error.append(JminusH_error2[index])
			HminusK_outlier.append(HminusK2[index])
			HminusK_outlier_error.append(HminusK_error2[index])


	target_index = RA1_matched3.index(candidate) #index of the target star in the original pm_RA stuff

	fout = open("output.txt",'w')
	fout.write("Sigma 3 detections\nRA UWISH2, DEC UWISH2, RA UKIDSS, DEC UKIDSS, pm RA, pm RA error, pm DEC, pm DEC error, pm_mag, pm mag error, J-K, J-K error, J-H, J-H error, H-K, H-K error\n")
	fout.write("%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f\n\n" %(RA1_matched3[target_index], DEC1_matched3[target_index], RA2_matched3[target_index], DEC2_matched3[target_index], pm_RA_component[target_index], pm_RA_std, pm_DEC_component[target_index], pm_DEC_std, pm_mag[target_index], pm_mag_std, JminusK[target_index], JminusK_error[target_index], JminusH[target_index], JminusH_error[target_index], HminusK[target_index], HminusK_error[target_index]))

	fout_len = len(RA_outlierU)
	for out in range(0, fout_len):
		fout.write("%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f\n" %(RA_outlierU[out], DEC_outlierU[out], RA_outlierF[out], DEC_outlierF[out], pm_RA_outlier[out], pm_RA_outlier_error[out], pm_DEC_outlier[out], pm_DEC_outlier_error[out], pm_mag_outlier[out], pm_mag_outlier_error[out], JminusK_outlier[out], JminusK_outlier_error[out], JminusH_outlier[out], JminusH_outlier_error[out], HminusK_outlier[out], HminusK_outlier_error[out]))

	fout.close()

	length = len(pm_RA_component2)
	print "\nThere are now %d matched stars after applying colour reddening limitations." % length
	print "\nHPM candidate data:"
	print "\nPM RA component (mas/yr)"
	print "\n\t", pm_RA_component[target_index], u'\u00B1', pm_RA_std
	print "\nPM DEC component (mas/yr)"
	print "\n\t", pm_DEC_component[target_index], u'\u00B1', pm_DEC_std
	print "\nPM Magnitudes (mas/yr)"
	print "\n\t", pm_mag2[target_index], u'\u00B1', pm_mag_std

	fig = plt.figure() #PLOTS DATA TO A FIGURE~~~~~~~~~~~~
	ax = fig.add_subplot(1,1,1)
	plt.plot(pm_RA_component2, pm_DEC_component2, 'kx') #THIS IS THE DATA WITHOUT CRAP COLOURS~~~~~~~~~~~~~~~~~~~~
	plt.plot(pm_RA_component[target_index], pm_DEC_component[target_index], 'rx', linewidth=10) #THIS IS THE ORIGINAL TARGET STAR BEING OUTPUT AS A RED X~~~~~~~~~~~
	plt.ylabel("PM DEC (mas/yr)")
	plt.xlabel("PM RA * COS(DEC) (mas/yr)")
	plt.axis('equal')
	plt.grid(True)
	ax.add_patch(sigma_1) #adds sd overplots
	ax.add_patch(sigma_2)
	ax.add_patch(sigma_3)
	ax.add_patch(sigma_4)
	ax.add_patch(sigma_5)
	plt.show()

	restart(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate)

def matching(RA1_order, RA2_order, DEC1_order, DEC2_order, distance1_order, distance2_order, JAperMag2_order, HAperMag2_order, K1AperMag2_order, J_error_order, H_error_order, K_error_order, h2_mag_order):

	print "\n####################################################################"
	print "#                          Matching                                #"
	print "####################################################################"

	global RA1_matched3 #variables used in the next function
	global RA2_matched3
	global DEC1_matched3
	global DEC2_matched3
	global JAperMag2_matched3
	global HAperMag2_matched3
	global K1AperMag2_matched3
	global J_error_matched3
	global H_error_matched3
	global K_error_matched3
	global candidate

	RA1_order = RA1_order.tolist() #converting numpy arrays to lists
	RA2_order = RA2_order.tolist()
	DEC1_order = DEC1_order.tolist()
	DEC2_order = DEC2_order.tolist()
	JAperMag2_order = JAperMag2_order.tolist()
	HAperMag2_order = HAperMag2_order.tolist()
	K1AperMag2_order = K1AperMag2_order.tolist()
	J_error_order = J_error_order.tolist()
	H_error_order = H_error_order.tolist()
	K_error_order = K_error_order.tolist()

	RA1_matched = []; RA1_matched2 = []; RA1_matched3 = []
	RA2_matched = []; RA2_matched2 = []; RA2_matched3 = []
	DEC1_matched = []; DEC1_matched2 = []; DEC1_matched3 = []
	DEC2_matched = []; DEC2_matched2 = []; DEC2_matched3 = []
	JAperMag2_matched = []; JAperMag2_matched2 = []; JAperMag2_matched3 = []
	HAperMag2_matched = []; HAperMag2_matched2 = []; HAperMag2_matched3 = []
	K1AperMag2_matched = []; K1AperMag2_matched2 = []; K1AperMag2_matched3 = []
	J_error_matched = []; J_error_matched2 = []; J_error_matched3 = []
	H_error_matched = []; H_error_matched2 = []; H_error_matched3 = []
	K_error_matched = []; K_error_matched2 = []; K_error_matched3 = []

	a = len(distance1_order); b = len(distance2_order) #defines value of lists

	print "\nMatching stars by 1 arcsecond."
	#matching_criteria_input = float(raw_input(">")) #dunno if I want to put it as 1 or as a user input
	matching_criteria = 1.0 / 3600.0 #converts arcsecond input into degrees

	h2_mag_order = h2_mag_order.tolist()
	#print "\nMatching target star: The H2 UWISH2 mag should be similar to UKIDSS K-band magnitude."
	print "\nPlease input the exact RA co-ordinates of the candidate star in the UWISH2 data."
	candidate = float(raw_input("> "))
	target_h2 = h2_mag_order[RA1_order.index(candidate)]
	target_RA = candidate
	target_DEC = DEC1_order[RA1_order.index(candidate)]

	count = 0
	for cand1 in range(0,b): #this bit matches the target star
		kminus = K1AperMag2_order[cand1] - 0.5
		kplus = K1AperMag2_order[cand1] + 0.5
		W = RA2_order[cand1] - matching_criteria
		X = RA2_order[cand1] + matching_criteria
		Y = DEC2_order[cand1] - matching_criteria
		Z = DEC2_order[cand1] + matching_criteria

		if W <= target_RA <= X: #this will work better as I can always ID the target

			if Y <= target_DEC <= Z:

				if kminus <= target_h2 <= kplus:
					count = count + 1

	check = len(RA2_matched)
	if count == 1:
			print "\nTarget star matched."

	else:
			print "\nTarget star cannot be matched. Exiting program"
			exit()

	#this bit makes sure the matched RA's are in a sensible range to each other
	#a = UWISH2 data
	#b = UKIDSS data

	for c in  range(0, b): #starts at 1 to skip the HPM candidate -> some examples with double peaks would delete the HPM candidate in later code because it would be matched twice
		g = RA2_order[c] - matching_criteria
		h = RA2_order[c] + matching_criteria
		i = DEC2_order[c] - matching_criteria
		j = DEC2_order[c] + matching_criteria
		u = RA2_order[c]
		kminus = K1AperMag2_order[c] - 0.5
		kplus = K1AperMag2_order[c] + 0.5

		for d in range(0, a):  #second list THIS IF THE FIRST MATCHING

			if g <= RA1_order[d] <= h: #if RA1_order[d] is in the range...

				if i <= DEC1_order[d] <= j:#and if DEC1_order[d] is also in the range...

					if kminus <= h2_mag_order[d] <= kplus:
						v = RA1_order[d]
						index_RA1 = RA1_order.index(v)#finds index for v in the list
						index_RA2 = RA2_order.index(u) #;print(c, d, index_RA1, index_RA2)
						RA1_matched.append(v) #and then the lists get appended according to the indexes and matched things
						RA2_matched.append(u)
						DEC1_matched.append(DEC1_order[index_RA1])
						DEC2_matched.append(DEC2_order[index_RA2])
						JAperMag2_matched.append(JAperMag2_order[index_RA2])
						HAperMag2_matched.append(HAperMag2_order[index_RA2])
						K1AperMag2_matched.append(K1AperMag2_order[index_RA2])
						J_error_matched.append(J_error_order[index_RA2])
						H_error_matched.append(H_error_order[index_RA2])
						K_error_matched.append(K_error_order[index_RA2])

	#the thing below removes all duplicates in the RA list and preserves order THIS BIT REMOVES DUPLIACTE RA'S AND MATCHES OTHER DATA
	c = Counter(RA2_matched)
	RA2_matched2_1 = [x for x in RA2_matched if c[x] == 1]
	f = len(RA2_matched2_1)
	f1 = len(RA2_matched)

	for k in range(0,f): #this can start from 0 because the HPM candidate isn't added to the list outside of matching
		l = RA2_matched2_1[k]

		for m in range(0,f1): #lists now has to be re-written to remove the other data for the duplicates that have just been deleted
			n = RA2_matched[m]

			if l == n:
				index = RA2_matched.index(l)
				RA2_matched2.append(RA2_matched[index])
				DEC2_matched2.append(DEC2_matched[index])
				HAperMag2_matched2.append(HAperMag2_matched[index])
				JAperMag2_matched2.append(JAperMag2_matched[index])
				K1AperMag2_matched2.append(K1AperMag2_matched[index])
				RA1_matched2.append(RA1_matched[index])
				DEC1_matched2.append(DEC1_matched[index])
				J_error_matched2.append(J_error_order[index])
				H_error_matched2.append(H_error_order[index])
				K_error_matched2.append(K_error_order[index])


	#now do it for the DEC list, I'm not 100% if this bit is needed but I think it is just in case! THIS BIT REMOVE DUPLIACE DECS AND MATCHES OTHER DATA
	c = Counter(DEC2_matched2)
	DEC2_matched2_1 = [x for x in DEC2_matched2 if c[x] == 1] #removes duplicates RA and preserves order
	f = len(DEC2_matched2_1)
	f1 = len(DEC2_matched)

	for k in range(0,f):
		l = DEC2_matched2_1[k]

		for m in range(0,f1):
			n = DEC2_matched[m]

			if l == n:
				index = DEC2_matched2.index(l)
				RA2_matched3.append(RA2_matched2[index])
				DEC2_matched3.append(DEC2_matched2[index])
				HAperMag2_matched3.append(HAperMag2_matched2[index])
				JAperMag2_matched3.append(JAperMag2_matched2[index])
				K1AperMag2_matched3.append(K1AperMag2_matched2[index])
				RA1_matched3.append(RA1_matched2[index])
				DEC1_matched3.append(DEC1_matched2[index])
				J_error_matched3.append(J_error_order[index])
				H_error_matched3.append(H_error_order[index])
				K_error_matched3.append(K_error_order[index])

	indexxx = RA1_matched3.index(candidate)
	print "Co-ordinate in UKIDSS: ", RA2_matched3[indexxx]
	print "J magnitude: ", JAperMag2_matched3[indexxx]
	print "H magnitude: ", HAperMag2_matched3[indexxx]
	print "K magnitude: ", K1AperMag2_matched3[indexxx]

	proper_motions(RA1_matched3, RA2_matched3, DEC1_matched3, DEC2_matched3, JAperMag2_matched3, HAperMag2_matched3, K1AperMag2_matched3, J_error_matched3, H_error_matched3, K_error_matched3, candidate)

def numericalSort(value): #function used for sorting AUTHOR: Martijn Pieters

	#The numericalSort function splits out any digits in a filename, turns it into
	#an actual number, and returns the result for sorting.

	numbers = re.compile(r'(\d+)')
	parts = numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts

def Froebrich(file_address):

	print "\n####################################################################"
	print "#                          Froebrich                               #"
	print "####################################################################"
	print "\nEpoch read in from file, no input is needed."

	f = open(file_address, 'r') #opening CSV file
	read = csv.reader(f) #assigning variable read as the CSV read module

	global RA2_order  #defining these variables as global variables his allows
	global DEC2_order #them to be used in the global code
	global epoch2
	global distance2_order
	global JAperMag2_order
	global HAperMag2_order
	global K1AperMag2_order
	global J_error
	global H_error
	global K_error

	RA2 = [] #defining lists for RA, DEC and epoch data
	DEC2 = []
	distance2 = []
	distance2_order = []
	JAperMag2 = []
	HAperMag2 = []
	K1AperMag2 = []
	J_error = []
	H_error = []
	K_error = []

	#print "\nHow many lines do you want to skip?"
	print "18 lines skipped."
	#skip = int(raw_input("> ")) #I think 18 lines always needs to be skipped
	skip = 18

	#this bit is just reading in the data
	for row in islice(read, skip, None): #skipping first 18 lines of csv file

		for row in read:
			RA2.append(row[2]) #adding data to earlier defined lists
			DEC2.append(row[3])
			distance2.append(row[18])
			JAperMag2.append(row[6])
			J_error.append(row[7])
			HAperMag2.append(row[8])
			H_error.append(row[9])
			K1AperMag2.append(row[10])
			K_error.append(row[11])
			epoch = row[17]

		#this is where data gets sorted
		epoch2 = float(epoch)
		distance2_order = sorted(distance2, key=numericalSort)
		l = len(distance2) #finds length of the distance1 (unsorted) list

		RA2_order = np.zeros(l) #creates  array full of zeros of length t
		DEC2_order = np.zeros(l)
		JAperMag2_order = np.zeros(l)
		HAperMag2_order = np.zeros(l)
		K1AperMag2_order = np.zeros(l)
		J_error_order = np.zeros(l)
		H_error_order = np.zeros(l)
		K_error_order = np.zeros(l)

		for j in range(0, l):
			if distance2_order[j] == distance2[j]:
				RA2_order[j] = RA2[j]
				DEC2_order[j] = DEC2[j]
				JAperMag2_order[j] = JAperMag2[j] #magnitudes included now because we're interested in those from UWISH2 data
				HAperMag2_order[j] = HAperMag2[j] #this bit is commented in the U09B8 function!
				K1AperMag2_order[j] = K1AperMag2[j]
				J_error_order[j] = J_error[j]
				H_error_order[j] = H_error[j]
				K_error_order[j] = K_error[j]

			else:
				y = distance2.index(distance2_order[j])
				RA2_order[j] = RA2[y]
				DEC2_order[j] = DEC2[y]
				JAperMag2_order[j] = JAperMag2[y]
				HAperMag2_order[j] = HAperMag2[y]
				K1AperMag2_order[j] = K1AperMag2[y]
				J_error_order[j] = J_error[y]
				H_error_order[j] = H_error[y]
				K_error_order[j] = K_error[y]

	f.close()

	print "First element:", RA2[0] #prints first element

	matching(RA1_order, RA2_order, DEC1_order, DEC2_order, distance1_order, distance2_order, JAperMag2_order, HAperMag2_order, K1AperMag2_order, J_error_order, H_error_order, K_error_order, h2_mag_order)

def U09B8(file_address):

	print "\n####################################################################"
	print "#                            U09B8                                 #"
	print "####################################################################"

	global RA1_order  #defining these variables as global variables
	global DEC1_order #this allows them to be used in the code
	global epoch1     #I think I could use return but whatever, I don't know what I'm doing
	global distance1_order
	global h2_mag_order

	#U09B8 data has no epoch but it will be the epoch of the image on the UWISH2
	#website. The code will ask for that date and convert into the right format.
	print "\nWhat year was the image of the candidate taken in?"
	year = int(raw_input("> "))
	print "\nAnd what month? (number)"
	month = int(raw_input("> "))
	epoch1 = year + (month / 12.0) #converts date into something friendly

	f = open(file_address, 'r') #opening CSV file
	read = csv.reader(f) #assigning variable read as the CSV read module

	RA1 = [] #defining lists for RA, DEC and epoch data
	DEC1 = []
	distance1 = []
	h2_mag = []
	distance1_order = []

	#print "\nHow many lines do you want to skip?"
	print "First 14 lines of file skipped."
	#skip = int(raw_input("> ")) #generally tends to be 13 or 14 lines that needs to be skipped
	skip = 14 #I got tired of inputting this. The target star was never the missed one..

	#this bit is just reading in the data
	for row in islice(read, skip, None): #skipping first 14 lines of csv file

		for row in read:
			RA1.append(row[3]) #adding data to earlier defined lists
			DEC1.append(row[4])
			distance1.append(row[45])
			h2_mag.append(row[28])

		#this bit is all about sorting the data
		distance1_order = sorted(distance1, key=numericalSort) #sorts distance in ascending order
		t = len(distance1) #finds length of the distance1 (unsorted) list

		RA1_order = np.zeros(t) #creates array full of zeros of length t
		DEC1_order = np.zeros(t)
		h2_mag_order = np.zeros(t)

		for i in range(0, t): #loop from range of 0 - t

			if distance1_order[i] == distance1[i]:
				RA1_order[i] = RA1[i]
				#if distance1_order has same element as distance1
				#then the element of RA1 is the same in RA1_order
				DEC1_order[i] = DEC1[i]
				#if distance1_order has same element as distance1
				#then the element of DEC1 is the same in
				#DEC1_order
				h2_mag_order[i] = h2_mag[i]

			else:
				x = distance1.index(distance1_order[i]) #finds index
				RA1_order[i] = RA1[x] #the next element belongs in this index
				DEC1_order[i] = DEC1[x]
				h2_mag_order[i] = h2_mag[x]

	f.close() #it's only polite to close it after opening it

	print "First element:", RA1[0] #prints first element

	Froebrich(file_address_Froebrich)

def file_address(address):

	global file_address_U09B8
	global file_address_Froebrich

	print "\nFile address: ", address
	print "\nPlease input the remaining file address with file extension for UWISH2 data."
	file_name_U09B8 = raw_input("> ")
	print "\nPlease input the remaining file address with file extension for UKIDSS GPS data."
	file_name_Froebrich = raw_input("> ")

	file_address_U09B8 = address + file_name_U09B8
	file_address_Froebrich = address + file_name_Froebrich

	U09B8(file_address_U09B8)

def start():

	#this is the function that loads up first. It asks what OS I'm using
	#and then sets up the default file path

	global address

	print "\n####################################################################"
	print "# Annyeonghaseyo, it's time to calculate some proper motions oppa! #"
	print "####################################################################"
	print "\nWhich operating system are you using?"
	print "\n1. Windows"
	print "2. Linux"
	OS = int(raw_input("> ")) #reads in user input as an integer
	OS = 1

	if OS == 1:
		address = "/Users/saultyevil/Google Drive/"
		file_address(address)

	elif OS == 2:
		address = "/home/saultyevil/Downloads/PH600/data/"
		file_address(address)

	else:
		print "\nInvalid input. Restarting program."
		start() #with an invalid input the program loops back

#starts the program. Functions are done backwards so the program can loop using the restart function defined first.
start()
