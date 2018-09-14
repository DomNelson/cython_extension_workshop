from __future__ import division
import numpy as np
import sys, os
from collections import defaultdict, Counter


class Pedigree:
    def __init__(self, pedfile):
        """
        Methods for analyzing pedigree files, with the columns specified
        in the header as:

            Ind    Father    Mother

        which are the only columns loaded, regardless of others present.
        """
        self.pedfile = pedfile
        alldata = np.genfromtxt(pedfile, skip_header=1, usecols=(0, 1, 2),
                                dtype=int)

        self.inds = alldata[:, 0]
        self.fathers = alldata[:, 1]
        self.mothers = alldata[:, 2]
        indices = range(len(alldata[:,0]))
        self.ind_dict = dict(zip(alldata[:,0], indices))

        ## Build a dictionary containing the parents and offspring of
        ## each individual
        self.parent_dict = dict(zip(self.inds, alldata[:, 1:3]))

        self.offspring_dict = defaultdict(list)
        for ind, father, mother in alldata:
            self.offspring_dict[father].append(ind)
            self.offspring_dict[mother].append(ind)

        ## Probands are individuals who are neither mothers nor fathers
        probands = set(self.inds).difference(set(self.fathers))
        self.probands = probands.difference(set(self.mothers))


    def ordered_lineage(self, ind):
        """
        Returns all ancestors of the given ind, along with the lengths of
        each lineage connecting to them
        """
        ordered_lineage = defaultdict(list)
        indnum = self.ind_dict[ind]
        current_generation = [self.mothers[indnum], self.fathers[indnum]]

        if current_generation[0] == current_generation[1] == 0:
            ordered_lineage[0].append(0)
            return ordered_lineage

        i = 0

        while len(current_generation) > 0:
            i += 1
            next_generation = []
            for ind in current_generation:
                ## If there are multiple paths, we save both lengths in the list
                ordered_lineage[ind].append(i)

            for ancestor in current_generation:
                indnum = self.ind_dict[ancestor]
                if self.mothers[indnum] != 0:
                    next_generation.append(self.mothers[indnum])
                if self.fathers[indnum] != 0:
                    next_generation.append(self.fathers[indnum])

            current_generation = next_generation

        return ordered_lineage


    def getlineage(self, ind):
        """
        Returns all the ancestors of an individual, including themselves,
        in a single list
        """
        lineage = [ind]
        father, mother = self.parent_dict[ind]

        if father != 0:
            lineage.extend(self.getlineage(father))

        if mother != 0:
            lineage.extend(self.getlineage(mother))

        return lineage


    def ordered_descendants(self, ind):
        """
        Returns a dictionary of all descendants of 'ind', with a corresponding
        list of path lengths from 'ind' to each descendant, in number of
        generations, as:

            descendants[ind] = [path_length_1, path_length_2, ... ]
        """
        descendants = defaultdict(list)
        descendants[ind] = [0]
        currentgen = self.offspring_dict[ind]

        i = 0
        while len(currentgen) > 0:
            i += 1
            nextgen = []
            for ind in currentgen:
                ## If there are multiple paths, we save all lengths in the list
                descendants[ind].append(i)
            for offspring in currentgen:
                nextgen.extend(self.offspring_dict[offspring])

            currentgen = nextgen

        return dict(descendants)


    def allowedinds(self, indlist):
        """
        Returns:

            commonancestors - list of all common ancestors between inds in
                                indlist
            cone_inds - dict with format {anc: set(descendants of anc)}
            ind_cones - dict with format {ind: set(ancs ancestral to ind)}
        """
        ancestorsets = []
        for ind in indlist:
            ancestorsets.append(set(self.getlineage(ind)))

        commonancestors = list(set.intersection(*ancestorsets))

        if len(commonancestors) == 0:
            print "Warning: No common ancestors!"

        cone_inds = {}
        ind_cones = defaultdict(set)
        for anc in commonancestors:
            tmp_dict = self.ordered_descendants(anc)
            tmp = set(tmp_dict.keys())
            ## We only need to track the individuals who could possibly climb
            ## to some point on the tree. To find these individuals, we take
            ## the descendants of one common ancestor, and take the union of
            ## its intersection with the lineage of every affected proband.
            ## repeat for every common ancestor to get the descent cones.
            tmp_intersect = set()
            for ancestorset in ancestorsets:
                tmp_intersect.update(tmp.intersection(ancestorset))
            tmp = tmp_intersect
            cone_inds[anc] = tmp

            ## Build dict of cone memberships for each individual
            for ind in tmp:
                ind_cones[ind].add(anc)

            ## Nonexistant inds, denoted by '0', have no descendats
            ind_cones[0] = set()

        return commonancestors, cone_inds, ind_cones
