#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 12:20:38 2021

@author: saimihirj
"""

# Bioinformatics - WEEK 1
# Where in the Genome Does Replication Begin?(Part 1)


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count


# Input: A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text

def FrequncyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern]+=1
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequncyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            pattern = key
            words.append(pattern)
    return words


# Input: A DNA string Pattern
# Output: The reverse complement of Pattern

def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern

def Reverse(Pattern):
    return Pattern[::-1]

def Complement(Pattern):
    base_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    pattern = ""
    for char in Pattern:
        pattern += base_pairs.get(char)
    return pattern