#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from cmath import pi
from os import path

from acoular import L_p, integrate
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from numpy import save, load, array, minimum


def calcNegativeOfIntervall(intervalls):
    """
    :param intervalls: must have the form [ [StartA1, EndA1], [StartA2, EndA2]  ]
    :return: Inverted version of intervalls
    """
    intervallsToZero = []
    start, stop = 0, 0

    # merge intervalls:
    nrEntries = len(intervalls)
    mergedIntervalls = []
    n = 0
    while n < nrEntries - 1:
        if intervalls[n][1]+1 == intervalls[n+1][0]:
            mergedIntervalls.append([intervalls[n][0], intervalls[n+1][1]])
            n += 2
        else:
            mergedIntervalls.append(intervalls[n])
            n += 1

    # Invert intervalls
    for i in mergedIntervalls:
        if intervalls[0] == 0:  # Intervall beinhaltet das aller erste Sample
            start = i[1] + 1
        else:
            stop = i[0] - 1
            intervallsToZero.append([start, stop])
            start = i[1] + 1


    return intervallsToZero

def getIntervallsByDegreeOverlapping(angleRes, angletracker):
    """
    :param angleArray: From the trigger the angles for each time sample
    :param angleRes: Width of the intervals (in multiples of pi). E.g. res = pi/2 -> method returns 4 Intervall
    groups. [[g1],[g2],[g3],[g4] ; g1 = [[start1, stop1], [start2,stop2], ...] start1 = nr of time sample where degree
     = 0, stop1 = nr of time sample where degree = pi/2 first time, start2 = nr of time sample where degree = 0 for
     second time
    :return: intervall in form of: [ [ [StartA1, EndA1], [StartA2, EndA2]  ],   [  [StartB1, EndB1], [...]  ]  ]
    StartA ->
    Nr of the first sample which is in group A, StartA2 Nr of the sample where for the second time angle is in group A.
    group A e.g. 0°->10°, group B 10.01° -> 20°
    """
    angleArray = angletracker.angle
    intervals = []
    cuttedIntervals = []  #
    if 2 * pi % angleRes >= 0.05:
        raise Exception(f"2*pi is no multiple of resolution, res is {angleRes / pi} pi")

    for i in range(int(2 * pi / angleRes)):  # Number of intervalls
        intervals.append([])  # -> intervals= [[start],[]]
        cuttedIntervals.append([])
    n = len(intervals) - 1  # Höchster Index von intervals
    currN = int(angleArray[0] // angleRes)  # Gibt an, welches Winkelintervall vervollständigt wird. currN <= n?
    intervals[currN].append([0, None])

    lenAngleArray = len(angleArray)
    if angletracker.rot_direction == -1:  # Winkel werden kleiner mit Drehung
        for i in range(1, lenAngleArray - 1):
            # Is the current angle the last one for the current Intervall?
            if angleArray[i + 1] >= angleRes * (currN - 1) and angleArray[i + 1] < angleRes * currN:
                intervals[currN][-1][1] = i

                intervals[currN - 1].append([i + 1, None])
                currN -= 1
            elif currN == 0:
                if angleArray[i + 1] > angleRes:  # Check whether there is transition from e.g. 1° to 358°
                    # (#decrasing angles)
                    # [_2_,1, 359,358] 2° > 10° -> False, [2,_1_, 359,358],  359 > 10° -> True
                    intervals[currN][-1][1] = i
                    intervals[n].append([i + 1, None])
                    currN = n
    else:  # Angles are increasing
        for i in range(1, lenAngleArray - 1):
            if angleArray[i + 1] >= angleRes * (currN + 1):
                intervals[currN][-1][1] = i

                intervals[currN + 1].append([i + 1, None])
                currN += 1
            elif currN == n:
                if angleArray[i + 1] < angleRes:
                    intervals[currN][-1][1] = i
                    intervals[0].append([i + 1, None])
                    currN = 0

    # insert last Sample
    #     inervals[  angleIntervall(lastSample)] = Index of last Sample
    intervals[int(angleArray[-1] // angleRes)][-1][1] = lenAngleArray - 1

    return intervals

def averageAndMinimum(mapFolder, splittedMapName, highestMapNr, minSoundVolume, nrInts):
    """
    Calculates the average and minimum of all given maps and saves the maps minMap, average and standing
    :param mapFolder: Path to the folder, where the maps are saved
    :param splittedMapName: Name of the map. It MUST end with the nr of the map. e.g.: Blabla2, BlaBla3 ,...
    :param highestMapNr: Highest number at the end of the mapname
    :param minSoundVolume: Maps which highest volume level is below this wll be deleted
    :param nrInts: Nr of intervalls in which the revolution has been divided
    """
    print(f"Averaging maps (filename is {splittedMapName}Nr.npy)\n")
    maps = []
    for i in range(highestMapNr + 1):

        filename = splittedMapName + str(i) + '.npy'
        loadedMap = load(path.join(mapFolder, filename))

        if L_p(loadedMap).max() < minSoundVolume:
            print(
                f"\n\n \033[91mWARNING: MAP FOR Nr {i} maximum volume = {L_p(loadedMap).max()}dB "
                f"<{minSoundVolume}dB skipping this map for calculating minimum and average map\u001b[0m\n\n")
            continue

        maps.append(loadedMap)
        print(f"max: {format(L_p(loadedMap).max(), '.2f')}dB")

    minMap = maps[0]
    averageMap = maps[0]

    for i in range(1, len(maps)):
        minMap = minimum(maps[i], minMap)
        averageMap += maps[i]

    standingMap = averageMap - minMap*nrInts
    save(path.join(mapFolder, splittedMapName + "Average"), averageMap)
    save(path.join(mapFolder, splittedMapName + "Standing"), standingMap)
    save(path.join(mapFolder, splittedMapName + "Min"), minMap)
    print("Done Minimizing & Averaging")

def invertIntervalls(intervalls, nrOfSamples):
    """
    Inverts the given Intervalls. E.g. intervalls = [[5,10], [20,25]] -> returns: [[0,4], [11,19], [26,nrOfSamples+1]]
    :param intervalls: must have the form [ [StartA1, EndA1], [StartA2, EndA2]  ]. StartA1 = starting samples of
    intervall 1, which microphonedata shall not later get zeroed. EndA1 = end sample of intervall 1
    :return: Inverted intervalls in form [[StartToZeroA1, EndToZero1], [...], ...]
    """
    print("invertIntervalls")

    intervallsToZero = []

    start, stop = 0, 0

    # ----- merge intervalls:
    lenInts = len(intervalls)
    if lenInts >1:
        mergedIntervalls = []
        n = 0
        while n < lenInts - 1:
            if intervalls[n][1]+1 == intervalls[n+1][0]:    # e.g. [[1,10], [11,16]]
                mergedIntervalls.append([intervalls[n][0], intervalls[n+1][1]])
                n += 2  # skip one entry in intervalls, because 2 intervalls got merged
            else:
                mergedIntervalls.append(intervalls[n])
                n += 1
        try:
            if intervalls[-2][1] +1 == intervalls[-1][0]:
                mergedIntervalls.pop()
                mergedIntervalls.append([intervalls[-2][0], intervalls[-1][1]])
            else:
                mergedIntervalls.append(intervalls[-1])
        except:
            pass
    else:
        mergedIntervalls = intervalls

    # ----- Invert intervalls
    lenMerged = len(mergedIntervalls)
    for i in range(lenMerged - 1):
        x = 1
        if mergedIntervalls[i][0] == 0:   # This intervalls has sample nr 0 in it
            start = mergedIntervalls[i][1] + 1        # intervalls[0]= [0, 20] -> start = 21
        else:   # 0 isn't inside i
            stop = mergedIntervalls[i][0] - 1         # i= [15,25] -> stop current intervall at 14
            intervallsToZero.append([start, stop])
            start = mergedIntervalls[i][1] + 1        # start next intervalls with 26
        # If i[1] = last sample -> intervallsToZero won't get appended with [lastsample+1, ?]. Because after start get's
        # set to start = lastSample+1 intervallsToZero.append won't get triggered

    lastInt = mergedIntervalls[-1]
    if len(mergedIntervalls) != 1:
        if lastInt[1] == nrOfSamples - 1:  # the last Intervall which shall not get zeroed ends with the last sample
            intervallsToZero.append([start, lastInt[0]])
        else:   # intervallNotToZero = [lastStart, LastStop] -> intervallsToZero.append([lastStop, nrOfSamples])
            intervallsToZero.append([start, lastInt[0]-1])
            intervallsToZero.append([lastInt[1]+1, nrOfSamples - 1])
    else:   # Just one entry in merged intervalls
        if lastInt[0] == 0:
            if lastInt[1] == nrOfSamples-1:
                print("C_FilterFunctionality/InvertIntervalls will zero no Data? Check if the inputs are right!")
            else:
                intervallsToZero.append([lastInt[1]+1, nrOfSamples-1])
        else: # lastInt doesn't start with sample nr 0
            intervallsToZero.append([0, lastInt[0]-1])
            if lastInt[1] != nrOfSamples:
                intervallsToZero.append([lastInt[1]+1,nrOfSamples - 1])


    return intervallsToZero

def plotMaps(mapFolder, splittedMapName, grid, minFactorMultiply=1, sectors= [], freq=1500):
    """
    Plots three Average, Standing and min Map maps
    :param mapFolder: Path to the map folder
    :param splittedMapName: Name of the maps
    :param grid:
    :param minFactorMultiply: Factor with which the min map gets multiplied. Should be one or equal to number of angle
    intverals
    :param sectors: [[A],[B],...], where A and B are sectors (draws circles for these sectors)
    :param freq: frequency of the map. Just used for labeling the maps
    """
    averageMapDB = L_p(load(path.join(mapFolder, splittedMapName + "Average" + ".npy")))
    standingMapDB = L_p(load(path.join(mapFolder, splittedMapName + "Standing" + ".npy")))

    multiply = True if minFactorMultiply != 1 else False
    minMap = load(path.join(mapFolder, splittedMapName + "Min" + ".npy"))
    minMap = minMap*minFactorMultiply if multiply else minMap
    minMapDB = L_p(minMap)
    print(f"minMap max = {minMapDB.max()} dB")
    print(f"average Map max = {averageMapDB.max()}")

    f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex='all')
    f.canvas.set_window_title(splittedMapName)

    averageMapDBmx = averageMapDB.max()
    dbRange = 17

    ax1.imshow(standingMapDB, vmax=averageMapDBmx, vmin=averageMapDBmx - dbRange, interpolation='nearest',
               extent=grid.extend(),
               origin='lower')
    ax1.set_title(f"Gemittelte Map - min Map\n \n {freq} Hz")
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")


    ax2.imshow(averageMapDB, vmax=averageMapDBmx, vmin=averageMapDBmx - dbRange, interpolation='nearest',
               extent=grid.extend(),
               origin='lower')
    ax2.set_title(f"Gemittelte Map\n\n {freq} Hz")
    ax2.set_xlabel("x [m]")
    ax2.set_ylabel("y [m]")

    im = ax3.imshow(minMapDB, vmax=averageMapDBmx, vmin=averageMapDBmx - dbRange, interpolation='nearest',
               extent=grid.extend(),
               origin='lower')
    text = f"Werte multipliziert mit {minFactorMultiply}" if multiply else ""
    ax3.set_title(f"min Map\n {text}\n {freq} Hz")
    ax3.set_xlabel("x [m]")
    ax3.set_ylabel("y [m]")
    cax = f.add_axes([ax3.get_position().x0+0.25, ax3.get_position().y0-0.1, 0.02, ax3.get_position().height+0.2])
    f.colorbar(im,cax=cax, label="Schalldruckpegel [dB]", orientation='vertical')

    for a in [ax1, ax2, ax3]:
        for index, s in enumerate(sectors):
            circle = Circle((s[0], s[1]), radius=s[2], fill=False, linestyle='--', color='w', linewidth=2.5)
            a.add_patch(circle)
            a.text(s[0]-s[2]/2, s[1]+s[2]+0.05, f"S {index}", color="w", size=20.)

    plt.show()

def integrateSources(positions, mapFolder, splittedMapName, highestMapNr, grid, freq=1500):
    """
    Integrates  over the given position. And plots the sound volumes for each map in a scatter plot
    :param grid:
    :param positions: [[x0,y0,r0], [x1,y1,r1]] x,y = positions of the expected sound source, radius of the sector to
    integrate
    :param mapFolder: path as string to the npy map folder
    :param splittedMapName: name of the maps
    :param highestMapNr: highest nr of the map
    """
    integrated = []
    sectors = []

    for xyr in positions:
        sectors.append(array([xyr[1], xyr[0], xyr[2]]))

    for index, s in enumerate(sectors):
        integrated.append([])
        for i in range(highestMapNr +1):

            filename = splittedMapName + str(i) + '.npy'
            loadedMap = load(path.join(mapFolder, filename))
            integrated[index].append(L_p(integrate(loadedMap, grid, s)))

    degreeRes = 360/(highestMapNr+1)
    angles = [a * degreeRes for a in list(range(len(integrated[0])))]
    #angles = list(range(len(integrated[0])))       # Use this if you want to plot the numbers of maps instead
    for index,i in enumerate(integrated):
        plt.plot(angles, i, "-o", label=f"Sektor S {index}")
    plt.legend()
    plt.title(f"Integrierte Schalldruckpegel \n {freq} Hz Terzband")
    plt.xlabel("Startwinkel der Map")
    plt.ylabel("Schalldruckpegel [dB]")
    plt.show()

