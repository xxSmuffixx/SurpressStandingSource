#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# imports from other packages
"""
Inside this script are SpatialInterpolatorRotationZeroing and ZeroedMaskedTimeSamples defined, which derivated from
acoular classes, normally found in tprocess.py and sources.py
"""

from acoular import SpatialInterpolatorRotation, TimeSamples, SamplesGenerator, Calib, MaskedTimeSamples
from numpy import pi
from traits.api import List


# ========== tprocess
from traits.traits import Trait, Property


class SpatialInterpolatorRotationZeroing(SpatialInterpolatorRotation):
    """
    Spatial  Interpolation for rotating sources. Gets samples from :attr:`source`
    and angles from  :attr:`AngleTracker`.Generates output via the generator :meth:`result`
    Result function returns zeroed microphonedata, for samples inside of intsToZero
    """

    intsToZero = List()
    def result(self, num=128):
        """
        Python generator that yields the output block-wise.

        Parameters
        ----------
        num : integer
            This parameter defines the size of the blocks to be yielded
            (i.e. the number of samples per block).

        Returns
        -------
        Samples in blocks of shape (num, :attr:`numchannels`).
            The last block may be shorter than num.
        """
        period = 2 * pi                         # period for rotation
        angle = self.angle_source._get_angle()  # get angle
        count = 0       # counter to track angle position in time for each block
        for timeData in self.source.partiallyZeroedResult(num, self.intsToZero):
            phiDelay = angle[count:count + num]
            interpVal = self._result_core_func(timeData, phiDelay, period, self.Q, interp_at_zero=False)
            yield interpVal
            count += num


# ========== source


class ZeroedMaskedTimeSamples(MaskedTimeSamples):
    def partiallyZeroedResult(self, num, intsToZero):
        """
        Yields data blockwise. If a sample inside the borders of "intsToZero" all microphone data of that sample gets
        set to 0. Because it zeroes Data, there will be a small error in the sound evaluation if e.g. originally the data
        of a single microphone  looked like this: [4,5,6,5,4,3,2,1] and it gets set to [4,5,6,0,0,0] there is a jump
        of -6 instead of -1
        :param num: blocksize of the yielded data
        :param intsToZero: [[start0, end0], [s1,e1],...] samples with the number from start to end will be zeroed
        :return: Samples in blocks of shape (num, numchannels)
        """
        if self.numsamples == 0:
            raise IOError("no samples available")

        # loop for all samples
        n = 0  # intsToZero[n]= [start n, end n]

        if self.start is None :
            i = 0
        else:
            i = self.start
        if self.stop is None:
            stop = self.data.shape[0]
        else:
            stop = self.stop

        cuttedLastIntervall = False
        doneWithLastZero = False
        while i < stop:
            dataToYield = self.data[i:  min(i + num, stop)][:, self.channels] # SELF.CALIB NOT IMPLEMENTED!!
            intervallsToZeroNewIndices = []
            # loop to calculate intervallsToZeroNewIndices
            newIndexesCalculated = False
            while not (newIndexesCalculated or doneWithLastZero):

                a = intsToZero[n][0]
                e = intsToZero[n][1]
                A = i
                E = i + num
                # Following lines calculate the intersection of a [A,E] & [a,e]
                if a <= A and e < E:
                    # Case 1&2: zero is balancing on on left side of group to yield
                    intervallsToZeroNewIndices.append([0, e - A])
                elif a > A and e < E:
                    # Case 3
                    intervallsToZeroNewIndices.append([a - A, e - A])
                elif a > A and e == E:
                    # Case 4: zero inside yield, touching right border
                    intervallsToZeroNewIndices.append([a - i, num - 1])
                elif a > A and e > E and a < E:
                    # Case 5: zero balancing on right side of yield
                    intervallsToZeroNewIndices.append([a - A, num - 1])
                    cuttedLastIntervall = True
                elif a == E and e > E:
                    # case 10 / 5.2: zero starts exactly at end inside of yield, hangs over on right side
                    intervallsToZeroNewIndices.append([num - 1, num - 1])
                    cuttedLastIntervall = True
                elif a == E and e == E:
                    # case 11 / 5.3: zero has len 1 an sits at the end of yield
                    intervallsToZeroNewIndices.append([num - 1, num - 1])
                elif a < i and e > E or a == i and e > E:
                    # Case 6 and 9: zero hanging over on left and right side of yield
                    # OR ero starts at start of yield, but hangs over on the right side
                    intervallsToZeroNewIndices.append([0, num - 1])
                    cuttedLastIntervall = True
                elif a < A and e < A:
                    # Case 7: Zero left of yield. This should only happen if the last couple samples dont have to get
                    # zeroed
                    newIndexesCalculated = True  # is this right?
                elif a > E and e > E:
                    # Case 8: toZero right of Group to yield this shouldnt get triggered?
                    cuttedLastIntervall = True  # dont increase n by one, we went too far with this

                if not cuttedLastIntervall:
                    n += 1
                    if n >= len(intsToZero):
                        doneWithLastZero = True
                else:  # eliminate the effect of n+=1, so for next run start again with n, not n+1
                    cuttedLastIntervall = False
                    newIndexesCalculated = True  # if we had to cut zero, this means we are at end of yield

            # create dataToYield
            for z in intervallsToZeroNewIndices:
                dataToYield[z[0]:z[1], :] = 0

            yield dataToYield
            i += num
