#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 20:10:14 2021

@author: saimihirj
"""
# Bioinformatics - WEEK 2
# Where in the Genome Does Replication Begin?(Part 2)

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


# Reproduce the PatternCount function here.
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    skew_array = [0]
    for i in range(len(Genome)):
        if Genome[i] == 'C':
            skew_array.append(skew_array[i] - 1)
        elif Genome[i] == 'G':
            skew_array.append(skew_array[i] + 1)
        else:
            skew_array.append(skew_array[i])
    return skew_array

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    # generate an empty list positions
    # set a variable equal to SkewArray(Genome)
    # find the minimum value of all values in the skew array
    # range over the length of the skew array and add all positions achieving the min to positions
    positions = [] # output variable
    skew_array = SkewArray(Genome)
    count = 0
    min_value = min(skew_array)
    for i in skew_array:
        if i == min_value:
            positions.append(count)
        count+=1    
    return positions

# Input:  A DNA string Genome
# Output: A list containing all integers i maximizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MaximumSkew(Genome):
    # generate an empty list positions
    # set a variable equal to SkewArray(Genome)
    # find the maximum value of all values in the skew array
    # range over the length of the skew array and add all positions achieving the max to positions
    positions = [] # output variable
    skew_array = SkewArray(Genome)
    count = 0
    max_value = max(skew_array)
    for i in skew_array:
        if i == max_value:
            positions.append(count)
        count+=1    
    return positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    count = 0
    for i, j in zip(p, q):
        if i != j:
            count += 1
    return count

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count+=1
    return count