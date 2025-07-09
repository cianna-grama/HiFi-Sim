# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# GenerateQscores.py

import random
import numpy as np

class GenerateQscores:
    '''
    Class to generate quality scores based on a given distribution.
    Distributions: random, normal, weighted.
    '''

    # initialize varibles for max/mins of hifi and phred scores
    def __init__(self):
        '''
        Initilizes variables: 
        - Maxium and minimum possible phred scores.
        - Maxium and minimum expected phred scores for HiFi seqeuencing.
        '''
        # possible phred scores
        self.minphred = 0
        self.maxphred = 93
        # expected scores 
        self.hifimin = 30
        self.hifimax = 40

    def checkphredscores(self, *scores):
        '''
        Checks that given scores are within the possible range of 
        phred scores. 
        Range: 0-93.
        '''
        for score in scores:
            if not (self.minphred <= score <= self.maxphred):
                raise ValueError(f"Phred score {score} must be between {self.minphred} and {self.maxphred}")

    def checklength(self, length):
        '''
        Checks that given length is not negative. 
        '''
        if length <= 0:
            raise ValueError("length must be positive")

    def randomdist(self, length, minphred=30, maxphred=40):
        '''
        Generates quality scores with a random distribution.
            length (int): Amount of quality scores to generate. 
            minphred (int): Minimum phred score to generate.
            maxphred (int): Maximum phred score to generate.
        '''

        # validify inputs
        self.checklength(length)
        self.checkphredscores(minphred, maxphred)
        if minphred > maxphred:
            raise ValueError(f"minphred {minphred} must be <= maxphred {maxphred}")

        # randomly generate phred scores in range and alter ASCII code to character 
        # then add to list of quality scores
        qualitychars = []
        for i in range(length):
            phredscore = random.randint(minphred, maxphred)
            qualitychar = chr(phredscore + 33)
            qualitychars.append(qualitychar)

        # return list of quality characters as a joined string
        return "".join(qualitychars)

    def normaldist(self, length, meanphred=35, stddev=2, minphred=30, maxphred=40):
        '''
        Generates quality scores with a normal distribution. 
            length (int): Amount of quality scores to generate. 
            meanphred (int): Average phred score of the normal distribution..
            stddev (int): Standard deviation of the normal distribution.
            minphred (int): Minimum phred score of the normal distribution.
            maxphred (int): Maximum phred score of the normal distribution.
        '''
        # check inputs 
        self.checklength(length)
        self.checkphredscores(meanphred, minphred, maxphred)
        
        if stddev <= 0:
            raise ValueError("standard deviation must be positive")
        
        qualitychars = []
        for i in range(length):
            # generate random score from normal distribution using numpy
            phredscore = np.random.normal(meanphred, stddev)
            # clip to specified range
            phredscore = max(minphred, min(maxphred, int(round(phredscore))))
            # transform to ASCII and add to list
            qualitychar = chr(phredscore + 33)
            qualitychars.append(qualitychar)

        # return list of characters as a joined string to be in hifi format
        return "".join(qualitychars)

    def weightedchoice(self, length, weights=None, phredrange=None):
        '''
        Generates quality scores with a weighted choice where 
        the inputted score range is made more likely that others. 
        Defaults to range 30-40 with weights on higher scores
        to mimic HiFi sequencing patters.
            length (int): Amount of quality scores to generate.
            weights (list): List of the weight values for each 
            score in the phredrange list. 
            phredrange (list): List of the possible scores in order. 
        '''
        self.checklength(length)
        # create score range list if none given 
        if phredrange is None:
            phredrange = list(range(30, 40))
    
        self.checkphredscores(*phredrange)
        # if no input for weights, set default
        # realistic hifi- weights favor higher quality scores
        if weights is None:
            weights = [1, 1, 2, 3, 5, 8, 10, 8, 5, 3]  # peaks around 36-37
        
        if len(weights) != len(phredrange):
            raise ValueError(f"weights list must have same number of elements as given phred range: {len(phredrange)} elements")

        # choose, transform, and return qualities as characters 
        qualitychars = []
        for i in range(length):
            phredscore = random.choices(phredrange, weights=weights)[0]
            qualitychar = chr(phredscore + 33)
            qualitychars.append(qualitychar)
        
        return "".join(qualitychars)
 
    def decodechars(self, qualitychars):
        '''
        Decode a string of quality scores from HiFi ASCII characters
        back to numeric phred scores.
        '''
        return [ord(char) - 33 for char in qualitychars]
    


def main():
    # testing class
    generate = GenerateQscores()
    length = 20

    print("=== Quality Score Generation Examples ===\n")

    print("1. Random (30-40):")
    randomdist = generate.randomdist(length)
    print(f"   Quality string: {randomdist}")
    print(f"   Phred scores: {generate.decodechars(randomdist)}\n")

    print("2. Normal Distribution (mean=35, std=2):")
    normaldist = generate.normaldist(length)
    print(f"   Quality string: {normaldist}")
    print(f"   Phred scores: {generate.decodechars(normaldist)}\n")

    print("3. Weighted Choice (favors higher scores):")
    weighted = generate.weightedchoice(length)
    print(f"   Quality string: {weighted}")
    print(f"   Phred scores: {generate.decodechars(weighted)}\n")

if __name__ == "__main__":
    main()
