#  Copyright (C) 2011 Oliver Stegle, modified by Joerg Hagmann
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import scipy as SP
import sys
import scipy.stats as st
from shutil import copyfile
from builtins import range


def estimate_q_values(PV,m=None,pi=1.0):
    """estimate q vlaues from a list of Pvalues
    this algorihm is taken from Storey, significance testing for genomic ...
    m: number of tests, (if not len(PV)), pi: fraction of expected true null (1.0 is a conservative estimate)
    """
    if m is None:
        m = len(PV) * 1.0
    else:
        m*=1.0
    lPV = len(PV)
    
    #1. sort pvalues
    PV = PV.squeeze()
    IPV = PV.argsort()
    PV  = PV[IPV]

    #2. estimate lambda
    if pi is None:
        lrange = SP.linspace(0.05,0.95,max(lPV/100.0,10))
        pil    = SP.double((PV[:,SP.newaxis]>lrange).sum(axis=0))/lPV
        pilr   = pil/(1.0-lrange)
        #ok, I think for SNPs this is pretty useless, pi is close to 1!
        pi =1.0
        #if there is something useful in there use the something close to 1
        if pilr[-1]<1.0:
            pi = pilr[-1]
            
    #3. initialise q values
    QV_ = pi * m/lPV* PV
    QV_[-1] = min(QV_[-1],1.0)
    #4. update estimate
    for i in range(lPV-2,-1,-1):
        QV_[i] = min(pi*m*PV[i]/(i+1.0),QV_[i+1])
    #5. invert sorting
    QV = SP.zeros_like(PV)
    QV[IPV] = QV_
    return QV


if __name__ == '__main__':

    if len(sys.argv)<3:
        print("p2qv pv.csv qv.csv")
        sys.exit(1)

    pv_file = sys.argv[1]
    qv_file = sys.argv[2]

    M = SP.loadtxt(pv_file,dtype='str')
    if len(M.shape)<2:
        #M= M[:,SP.newaxis]
        copyfile(pv_file, qv_file)
        sys.exit()
    
    R = SP.empty([M.shape[0],M.shape[1]+1],dtype='object')
    R[:,0:-1] = M

    pv_str = M[:,-1]
    pv_str[pv_str=='NA'] = 'NAN'

    #ok convert to qv
    pv = SP.array(pv_str,dtype='float')
    pv[SP.isnan(pv)] = 1.0

    qv = estimate_q_values(pv)
    R[:,-1] = qv

    SP.savetxt(qv_file,R,fmt='%s',delimiter='\t')
