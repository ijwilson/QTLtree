#!/usr/bin/python
""" Generate a data set to use in QTLtree

This is a very simple simpulation algorithm to produce data sets that can be read into
QTLtree, without sending vast files and using only python, which is fairly generally available

"""
import random










print(random.random())



nSNPS = 1000
nIND = 200

betaprior=[0.4,0.4]

freq = [random.betavariate(betaprior[0],betaprior[1]) for i in range(nSNPS)]






