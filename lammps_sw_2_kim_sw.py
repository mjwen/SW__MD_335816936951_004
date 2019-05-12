#!/usr/bin/env python
'''
This program reads in a file in LAMMPS pair_style sw format and if possible
generates a parameter file for the KIM model driver for the Stillinger-Weber
potential (see https://openkim.org/MD_335816936951_004).

Note that the KIM model driver is less general than the LAMMPS format and
therefore not all conversions are possible.
'''

################################################################################
#
#  CDDL HEADER START
#
#  The contents of this file are subject to the terms of the Common Development
#  and Distribution License Version 1.0 (the "License").
#
#  You can obtain a copy of the license at
#  http:# www.opensource.org/licenses/CDDL-1.0.  See the License for the
#  specific language governing permissions and limitations under the License.
#
#  When distributing Covered Code, include this CDDL HEADER in each file and
#  include the License file in a prominent location with the name LICENSE.CDDL.
#  If applicable, add the following below this CDDL HEADER, with the fields
#  enclosed by brackets "[]" replaced with your own identifying information:
#
#  Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
#  CDDL HEADER END
#
#  Copyright (c) 2017, Regents of the University of Minnesota.
#  All rights reserved.
#
#  Contributor(s):
#     Ellad B. Tadmor
#
################################################################################

# Python 2-3 compatible code issues
from __future__ import print_function
try:
   input = raw_input
except NameError:
   pass

# Packages used
import math
import datetime
import sys

__version__ = "1.0.0"
__author__  = "Ellad Tadmor"
__date__    = "26-May-2018"

################################################################################
#
#   FUNCTIONS
#
################################################################################
def read_sw_triplets(swfile):
    '''
    Read in file in LAMMPS pair_style sw format.

    Line format:
    # comments
    spec1 spec2 spec3 epsilon sigma a lambda gamma  cos(theta) A B p q tol
    '''
    swspecs  = []
    swparams = []
    with open(swfile) as f:
        for line in f:
            linestrip = line.strip()
            if linestrip != '':
                if linestrip[0] != '#':
                    data = linestrip.split()
                    swspecs.append(data[0:3])
                    swparams.append([float(d) for d in data[3:14]])
    # Identify species in file
    flatswspec = [item for sublist in swspecs for item in sublist]
    speclist = []
    [speclist.append(item) for item in flatswspec if item not in speclist]
    print('Species in SW file: ',speclist)
    numspec = len(speclist)
    if numspec**3 != len(swspecs):
        raise Exception('Incorrect number of lines in input file.')
    return speclist, swspecs, swparams

def process_sw_triplets(speclist, swspecs, swparams):
    '''
    Process SW pair_style parameters in swspecs and swparams into
    SW KIM driver format.

    Line format:
    # comment line
    A B p q sigma lambda gamma costheta_0 cutoff

    where relative to the LAMMPS parameters:

    A := A*epsilon
    lambda := lambda*epsilon
    gamma  := gamma*sigma
    rcut   := a*sigma

    '''
    # Initialize parameter array
    numspecs = len(speclist)
    swparams_kim = [[None for x in range(numspecs)] for y in range(numspecs)]

    # Extract pairwise interactions from LAMMPS data
    for i in range(numspecs):
        for j in range(i,numspecs):
            # Find corresponding lines in SW file
            numfound = 0
            for n in range(len(swspecs)):
               if (swspecs[n][0]==speclist[i] and \
                   swspecs[n][1]==speclist[j] and \
                   swspecs[n][2]==speclist[j]) or \
                  (swspecs[n][0]==speclist[j] and \
                   swspecs[n][1]==speclist[i] and \
                   swspecs[n][2]==speclist[i]) :
                   kim_A        = swparams[n][6]*swparams[n][0]
                   kim_B        = swparams[n][7]
                   kim_p        = swparams[n][8]
                   kim_q        = swparams[n][9]
                   kim_sigma    = swparams[n][1]
                   kim_lambda   = swparams[n][3]*swparams[n][0]
                   kim_gamma    = swparams[n][4]*swparams[n][1]
                   kim_costheta = swparams[n][5]
                   kim_cutoff   = swparams[n][2]*swparams[n][1]
                   numfound += 1
                   if numfound == 1:
                       swparams_kim[i][j] = \
                           [ kim_A, kim_B, kim_p, kim_q, kim_sigma, kim_lambda,
                             kim_gamma, kim_costheta, kim_cutoff ]
                   elif numfound == 2:
                       if kim_A        != swparams_kim[i][j][0] or \
                          kim_B        != swparams_kim[i][j][1] or \
                          kim_p        != swparams_kim[i][j][2] or \
                          kim_q        != swparams_kim[i][j][3] or \
                          kim_sigma    != swparams_kim[i][j][4] or \
                          kim_lambda   != swparams_kim[i][j][5] or \
                          kim_gamma    != swparams_kim[i][j][6] or \
                          kim_costheta != swparams_kim[i][j][7] or \
                          kim_cutoff   != swparams_kim[i][j][8]:
                           print('Lines for ',speclist[i],' and ',speclist[j], \
                                 ' have in consistent parameters in LAMMPS SW file.')
                           raise Exception('Corrupt input file.')
            swparams_kim[j][i] = swparams_kim[i][j]
            if i == j and numfound != 1:
                print('Wrong number of lines containing just ',speclist[i],'. ' \
                      'There should be one such lines.')
                raise Exception('Corrupt input file.')
            if i != j and numfound != 2:
                print('Wrong number of lines containing ',speclist[i],' and ',speclist[j], \
                      ' with one of them twice at end. There should be two such lines.')
                raise Exception('Corrupt input file.')

    # Make sure three-body terms are consistent with the algorithm used by KIM driver:

    # Verify tolerance is zero
    for p in swparams:
        if p[10] != 0.0:
            print('Tolerance not set to zero in LAMMPS potential file.')
            raise Exception('LAMMPS file more general than supported by KIM driver.')

    # Verify costheta_0 is constant
    for n in range(1,len(swspecs)):
        if swparams[n][5] != swparams[0][5]:
            print('Parameter costheta_0 not constant in LAMMPS potential file.')
            raise Exception('LAMMPS file more general than supported by KIM driver.')

    # Verify that lambda_ijk = sqrt(lambda_ij)*sqrt(lambda_ik) is equalt to
    # epsilon_ijk*lambda_ijk in LAMMPS file
    differences_detected = False
    max_diff = 0.0
    for n in range(len(swspecs)):
        swspecnum = [ speclist.index(swspecs[n][isp]) for isp in range(3) ]
        i = swspecnum[0]
        j = swspecnum[1]
        k = swspecnum[2]
        kim_lambda_ijk = math.sqrt(swparams_kim[i][j][5]*swparams_kim[i][k][5])
        if abs(kim_lambda_ijk - swparams[n][0]*swparams[n][3]) > 1e-6:
            differences_detected = True
            print('WARNING: Three-body prefactor differs between KIM and LAMMPS file:')
            print('         Species = ',swspecs[n])
            print('         KIM                lambda_ijk = ',kim_lambda_ijk)
            print('         LAMMPS epsilon_ijk*lambda_ijk = ',swparams[n][0]*swparams[n][3])
            max_diff = max(max_diff,abs(kim_lambda_ijk - swparams[n][0]*swparams[n][3]))
    if differences_detected:
        print('')
        print('CAREFUL: Differences detected between pair_style parameterization')
        print('         and KIM parameterization.')
        print('')
        print('         Maximum difference between parameters = ',max_diff)

    return swparams_kim

def write_kim_paramfile(speclist, swparams_kim, swfile):
    '''
    Write SW parameter file conforming to the SW KIM driver format
    '''
    kimfile = 'SW_' + swfile.strip().rsplit( ".", 1 )[ 0 ] + '.params'
    numspecs = len(speclist)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(kimfile,'w') as f:
       f.write('# SW Parameter file for ' + ' '.join(speclist) + \
               ' generated from LAMMPS pair_style sw file' + \
               ' `' + swfile.strip() + '` ' + \
               ' by lammps_sw_2_kim_sw at ' + timestamp + '\n')
       f.write('\n')
       f.write('{}\n'.format(numspecs))
       f.write('\n')
       head = ['A(energy)', 'B', 'p',  'q',  'sigma(length)', 'lambda(energy)', \
               'gamma(length)', 'costheta_0', 'cutoff(length)' ]
       f.write('#')
       f.write('{:>23} '.format(head[0]))
       for k in range(1,9):
           f.write('{:>24} '.format(head[k]))
       f.write('\n')
       for i in range(numspecs):
           for j in range(i,numspecs):
               for k in range(9):
                   f.write('{:24.16e} '.format(swparams_kim[i][j][k]))
               f.write('\n')
       f.write('\n')
       f.write('#   First line, number of species\n')
       f.write('#\n')
       f.write('#   Each line lists the following 9 parameters for the interaction between two species:\n')
       f.write('#\n')
       f.write('#   A(energy)  B  p  q  sigma(length)  lambda(energy)  gamma(length)  costheta_0  cutoff(length)\n')
       f.write('#\n')
       f.write('#   In the following order:\n')
       for i in range(3):
           for j in range(i,numspecs):
               f.write('#   ({}, {})\n'.format(i+1,j+1))
       if numspecs>2:
           f.write('#    .\n')
           f.write('#    .\n')
           f.write('#    .\n')
       f.write('#   where ')
       if numspecs==1:
           f.write('{}={}.\n'.format(1,speclist[0]))
       else:
           for n in range(numspecs-1):
               f.write('{}={}, '.format(n+1,speclist[n]))
           f.write('and {}={}.\n'.format(numspecs,speclist[numspecs-1]))
       f.write('#\n')
       f.write('#   Note that the parameters in this file are not the standard Stillinger-Weber\n')
       f.write('#   potential parameters. Four of them are redefined according to:\n')
       f.write('#\n')
       f.write('#   A      := A*epsilon\n')
       f.write('#   lambda := lambda*epsilon\n')
       f.write('#   gamma  := gamma*sigma\n')
       f.write('#   cutoff := a*sigma\n')
       f.write('#\n')
       f.write('#   where the left-hand side are the parameters used here, and the right-hand\n')
       f.write('#   side are the standard Stillinger-Weber potential parameters.\n')

################################################################################
#
#   MAIN PROGRAM
#
###############################################################################
if __name__ == '__main__':
    # Read in SW triplet lines
    swfile = input("LAMMPS pair_style sw filename = ")
    try:
        speclist, swspecs, swparams = read_sw_triplets(swfile)
    except Exception as e:
        print(e)
        print('Translation failed.')
        sys.exit()

    # Translate LAMMPS pair_style SW to KIM SW driver format
    try:
        swparams_kim = process_sw_triplets(speclist, swspecs, swparams)
    except Exception as e:
        print(e)
        print('Translation failed.')
        sys.exit()

    # Write parameter file in SW KIM driver format
    write_kim_paramfile(speclist, swparams_kim, swfile)
