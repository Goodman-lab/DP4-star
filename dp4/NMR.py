# -*- coding: utf-8 -*-
"""
Created on Mon Jan  12 14:42:47 2015

@author: ke291

Takes care of all the NMR description interpretation, equivalent atom
averaging, Boltzmann averaging and DP4 input preparation and running DP4.py. Called by PyDP4.py

FUNCTIONS AFTER REWRITE:
Calculation of NMR shifts based on TMS reference
Equivalent atom averaging
NMR description parsing
NMR raw data interpretation top level organization
"""

import re
import os
import math
import pickle
import pkg_resources

from dp4.matching_assignment import find_assignment
from pathlib import Path

gasConstant = 8.3145
temperature = 298.15
hartreeEnergy = 2625.499629554010

# Data structure for loading and keeping all of experimental NMR data in one place.

class NMRData:
    def __init__(self,settings):

        self.cwd = Path(os.getcwd())
        self.InputPath = settings.NMRsource  # Initial structure input file
        self.Type = 'desc'          # desc or fid, depending on whether the description or raw data used
        self.Atoms = []             # Element labels
        self.Cshifts = []           # Experimental C NMR shifts
        self.Clabels = []           # Experimental C NMR labels, if any
        self.Hshifts = []           # Experimental H NMR shifts
        self.Hlabels = []           # Experimental H NMR labels, if any
        self.Equivalents = []       # Atoms assumed to be NMR equivalent in computational data
        self.Omits = []
        self.protondata = {}
        self.carbondata = {}

        print(self.InputPath)

        #print(self.InputPath.split('.'))

        #quit()

        if len(self.InputPath) == 0:
            print('No NMR Data Added, quitting...')
            quit()

        else:

            for ind1 , p in enumerate(self.InputPath):

                if p.exists():

                    if p.is_dir():
                        print('NMR data in description format required.')
                        quit()
                    else:
                        self.Type = 'desc'
                        self.ExpNMRFromDesc()

                else:
                    print('NMR data path does not exist, quitting...')
                    quit()

    def ExpNMRFromDesc(self):

        print('Loading NMR data from ' + str(self.InputPath))

        # Reads the experimental NMR data from the file
        ExpNMR_file = open(self.InputPath[0], 'r')
        Cexp = ExpNMR_file.readline()
        ExpNMR_file.readline()
        Hexp = ExpNMR_file.readline()

        # Check if exp NMR file contains info about equivalent atoms and read it
        # into an array
        # Also reads a list of atoms to omit from analysis

        equivalents = []
        omits = []

        ExpNMR_file.readline()
        for line in ExpNMR_file:
            if not 'OMIT' in line and len(line) > 1:
                equivalents.append(line[:-1].split(','))
            elif 'OMIT' in line:
                omits.extend(line[5:-1].split(','))

        ExpNMR_file.close()

        self.Clabels, self.Cshifts = self.ParseExp(Cexp)
        self.Hlabels, self.Hshifts = self.ParseExp(Hexp)
        self.Equivalents = equivalents
        self.Omits = omits

    def ParseExp(self, exp):

        if len(exp)>0:

            # Replace all 'or' and 'OR' with ',', remove all spaces and 'any'
            texp = re.sub(r"or|OR", ',', exp, flags=re.DOTALL)
            texp = re.sub(r" |any", '', texp, flags=re.DOTALL)

            # Get all assignments, split mulitassignments
            expLabels = re.findall(r"(?<=\().*?(?=\)|;)", texp, flags=re.DOTALL)
            expLabels = [x.split(',') for x in expLabels]

            # Remove assignments and get shifts


            ShiftData = (re.sub(r"\(.*?\)", "", exp.strip(), flags=re.DOTALL)).split(',')

            expShifts = [float(x) for x in ShiftData]

        else:

            expLabels = []
            expShifts=[]

        return expLabels, expShifts

def CalcBoltzmannWeightedShieldings(Isomers):

    energies = []

    for i, iso in enumerate(Isomers):


        # Calculate rel. energies in kJ/mol
        minE = min(iso.DFTEnergies)

        relEs = []

        for e in iso.DFTEnergies:
            relEs.append((e - minE) * hartreeEnergy)

        Isomers[i].Energies = relEs

        populations = []

        # Calculate Boltzmann populations
        for e in relEs:
            populations.append(math.exp(-e * 1000 / (gasConstant * temperature)))

        q = sum(populations)

        for p in range(0, len(populations)):
            populations[p] = populations[p] / q

        Isomers[i].Populations = populations

        # Calculate Boltzmann weighed shielding constants
        # by summing the shifts multiplied by the isomers population
        BoltzmannShieldings = []

        for atom in range(len(iso.Atoms)):
            shielding = 0

            c = 1

            for population, shieldings in zip(iso.Populations, iso.ConformerShieldings):

                c+=1

                shielding = shielding + shieldings[atom] * population

            BoltzmannShieldings.append(shielding)

        Isomers[i].BoltzmannShieldings = BoltzmannShieldings

    return Isomers


def GetTMSConstants(settings):
    TMSfile_path = pkg_resources.resource_filename(__name__, 'data/TMSdata')
    with open(TMSfile_path, 'r') as TMSfile:
        TMSdata = TMSfile.readlines()
    for i, line in enumerate(TMSdata):
        buf = line.split(' ')
        if len(buf) > 1:
            if settings.Solvent != '':
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == settings.Solvent.lower():
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    TMS_SC_C13 = float(buf[3])
                    TMS_SC_H1 = float(buf[4])
                    return TMS_SC_C13, TMS_SC_H1
            else:
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == 'none':
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    TMS_SC_C13 = float(buf[3])
                    TMS_SC_H1 = float(buf[4])
                    return TMS_SC_C13, TMS_SC_H1

    print("No TMS reference data found for these conditions, using defaults\n")
    print("Unscaled shifts might be inaccurate, use of unscaled models is" + \
          " not recommended.")

    return settings.TMS_SC_C13, settings.TMS_SC_H1


def NMRDataValid(Isomers):

    for isomer in Isomers:
        if (len(isomer.ConformerShieldings) == 0):
            return False

    return True


def CalcNMRShifts(Isomers, settings):

    print('WARNING: NMR shift calculation currently ignores the instruction to exclude atoms from analysis')
    for i, iso in enumerate(Isomers):

        BShieldings = iso.BoltzmannShieldings

        Cvalues = []
        Hvalues = []
        Clabels = []
        Hlabels = []

        for a, atom in enumerate(iso.Atoms):

            if atom == 'C':
                shift = (settings.TMS_SC_C13-BShieldings[a]) / (1-(settings.TMS_SC_C13/10**6))
                Cvalues.append(shift)
                Clabels.append('C' + str(a + 1))
        
            elif atom == 'H':
                shift = (settings.TMS_SC_H1-BShieldings[a]) / (1-(settings.TMS_SC_H1/10**6))
                Hvalues.append(shift)
                Hlabels.append('H' + str(a + 1))

        Isomers[i].Cshifts = Cvalues
        Isomers[i].Hshifts = Hvalues

        Isomers[i].Clabels = Clabels
        Isomers[i].Hlabels = Hlabels
    
        if settings.save_shift_data:
            shift_data = {
                'Cshifts': Cvalues,
                'Hshifts': Hvalues,
                'Clabels': Clabels,
                'Hlabels': Hlabels
            }
            with open(f'Isomer_{i}_shift_data.pkl', 'wb') as f:
                pickle.dump(shift_data, f)

        print('C shifts for isomer ' + str(i) + ": ")
        print(', '.join(['{0:.3f}'.format(x) for x in Isomers[i].Cshifts]))

        print('H shifts for isomer ' + str(i) + ": ")
        print(', '.join(['{0:.3f}'.format(x) for x in Isomers[i].Hshifts]))

        for conf in iso.ConformerShieldings:

            Cconfshifts = []
            Hconfshifts = []

            for a, atom in enumerate(iso.Atoms):

                if atom == 'C':

                    shift = (settings.TMS_SC_C13-conf[a]) / (1-(settings.TMS_SC_C13/10**6))
                    Cconfshifts.append(shift)

                if atom == 'H':
                    shift = (settings.TMS_SC_H1 - conf[a]) / (1 - (settings.TMS_SC_H1 / 10 ** 6))
                    Hconfshifts.append(shift)

            Isomers[i].ConformerCShifts.append(Cconfshifts)
            Isomers[i].ConformerHShifts.append(Hconfshifts)

    return Isomers


def PrintConformationData(AllSigConfs):
    # Make a list of populations and corresponding files for reporting
    # significant conformations
    """from operator import itemgetter
    ConfsPops = [list(x) for x in zip(args, populations)]
    ConfsPops.sort(key=itemgetter(1), reverse=True)
    totpop = 0
    i = 0
    while totpop < 0.8:
        totpop += ConfsPops[i][1]
        i += 1
    SigConfs = ConfsPops[:i]"""
    for Es, pops in zip(RelEs, populations):
        print('\nConformer relative energies (kJ/mol): ' + \
            ', '.join(["{:5.2f}".format(float(x)) for x in Es]))

        print('\nPopulations (%): ' + \
            ', '.join(["{:4.1f}".format(float(x)*100) for x in pops]))

    for i, SigConfs in enumerate(AllSigConfs):
        print("\nNumber of significant conformers for isomer "\
            + str(i+1) + ": " + str(len(SigConfs)) + "\n(pop, filename)")
        for conf in SigConfs:
            print("   " + format(conf[1]*100, "4.2f") + "%   " + conf[0])
        print('----------------')
        print("   " + format(100*sum([x[1] for x in SigConfs]), "4.2f") +\
            "%   in total")


def RemoveEquivalents(Noutp, equivs, OldCval, OldHval, OldClabels, OldHlabels):
    Cvalues = list(OldCval)
    Hvalues = list(OldHval)
    Clabels = list(OldClabels)
    Hlabels = list(OldHlabels)
    
    for eqAtoms in equivs:

        eqSums = [0.0]*Noutp
        eqAvgs = [0.0]*Noutp

        if eqAtoms[0][0] == 'H':
            #print eqAtoms, Hlabels
            for atom in eqAtoms:
                eqIndex = Hlabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Hvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Hlabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Hvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Hlabels.index(eqAtoms[atom])
                del Hlabels[del_index]
                for ds in range(0, Noutp):
                    del Hvalues[ds][del_index]

        if eqAtoms[0][0] == 'C':
            for atom in eqAtoms:
                eqIndex = Clabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Cvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Clabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Cvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Clabels.index(eqAtoms[atom])
                del Clabels[del_index]
                for ds in range(0, Noutp):
                    del Cvalues[ds][del_index]
                    
    return Cvalues, Hvalues, Clabels, Hlabels
    

def MAE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append(abs(L1[i]-L2[i]))
        return sum(L)/len(L)


def RMSE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append((L1[i]-L2[i])**2)
        return math.sqrt(sum(L)/len(L))


def ReadNMRShifts(Isomers, settings):
    for i, iso in enumerate(Isomers):
        with open(f'Isomer_{i}_shift_data.pkl', 'rb') as f:
            shift_data = pickle.load(f)
            Isomers[i].Cshifts = shift_data['Cshifts']
            Isomers[i].Hshifts = shift_data['Hshifts']

            Isomers[i].Clabels = shift_data['Clabels']
            Isomers[i].Hlabels = shift_data['Hlabels']
    return Isomers

def PairwiseAssignment(Isomers,NMRData,settings):

    # for each isomer sort the experimental and calculated shifts

    for i, iso in enumerate(Isomers):

        sortedCCalc = sorted(iso.Cshifts, reverse=True)

        sortedClabels =[  '' for i in iso.Clabels]

        for ind_1 ,  shift in enumerate(iso.Cshifts):

            ind_2 = sortedCCalc.index(shift)

            sortedClabels[ind_2] = iso.Clabels[ind_1]


        sortedHCalc = sorted(iso.Hshifts, reverse=True)

        sortedHlabels = ['' for i in iso.Hlabels]

        for ind_1, shift in enumerate(iso.Hshifts):

            ind_2 = sortedHCalc.index(shift)

            sortedHlabels[ind_2] = iso.Hlabels[ind_1]
        
        label_map = {}
        for shift, label in zip(sortedHCalc, sortedHlabels):
            label_map[shift] = label


        sortedCExp = sorted(NMRData.Cshifts, reverse=True)
        sortedHExp = sorted(NMRData.Hshifts, reverse=True)

        assignedCExp = [''] * len(sortedCCalc)
        assignedHExp = [''] * len(sortedHCalc)

        tempCCalcs = list(iso.Cshifts)
        tempHCalcs = list(iso.Hshifts)

        # do the assignment in order of chemical shift starting with the largest

        # Carbon

        for exp, shift ,label in zip(sortedCExp, sortedCCalc , sortedClabels):

            if label not in NMRData.Omits:

                ind = tempCCalcs.index(shift)

                assignedCExp[ind] = exp

                tempCCalcs[ind] = ''

        # Proton

        if settings.matching_assignment:
            assigned, exp_exclude, calc_exclude = find_assignment(sortedHExp, sortedHCalc, threshold=settings.threshold)
            iso.number_excluded = len(exp_exclude)
            print(f'\nIsomer {i}:')
            print(f'Experimental shifts excluded: {exp_exclude}')
            print(f'Calculated shifts excluded: {calc_exclude}')
            iso.Hshifts = []
            iso.Hexp = []
            iso.Hlabels = []
            for calc_shift, exp_shift in assigned.items():
                iso.Hshifts.append(calc_shift)
                iso.Hexp.append(assigned[calc_shift])
                iso.Hlabels.append(label_map[calc_shift])

        else:

            for exp,shift,label in zip(sortedHExp, sortedHCalc,sortedHlabels):

                if label not in NMRData.Omits:

                    ind = tempHCalcs.index(shift)

                    assignedHExp[ind] = exp

                    tempHCalcs[ind] = ''

            # update isomers class

            iso.Cexp = assignedCExp
            iso.Hexp = assignedHExp

    return Isomers

def number_exclusions(Isomers):
    print('\nNumber of exclusions:')
    for i, iso in enumerate(Isomers):
        print(f'Isomer {i + 1}: {iso.number_excluded}')
    print('\n')





