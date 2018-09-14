import sys, os
import numpy as np
import time

import ped
import pysignal


# pedfile = os.path.expanduser('~/project/anc_finder/data/BALasc_probands1930.txt')
# pedfile = os.path.expanduser('~/project/anc_finder/data/pedEx.txt')
pedfile = os.path.expanduser(
        '~/project/anc_finder/scripts/test/test_data/pedEx2.txt')
# pedfile = os.path.expanduser(
#         '~/project/anc_finder/scripts/test/test_data/pedEx3.txt')
P = ped.Pedigree(pedfile)

ped_arr = pysignal.sort_ped(P)
ninds = len(P.inds)
print len(P.probands), "probands"
sample_idx = [P.ind_dict[x] for x in list(P.probands)[:10]]
samples = np.array(sample_idx, dtype=np.int32)
print "Samples:", samples

cP = pysignal.cPed()
cP.load_ped(ped_arr, len(samples))
cP.load_samples(samples)
# cP.print_samples()
# for s in samples:
#     cP.update_ancestor_weights(s, 1)
# cP.print_nodes()
cP.init_sample_weights()

# cP.set_all_weights(1)
cP.print_nodes()
cP.climb_step()
cP.print_nodes()
cP.climb_step()
cP.print_nodes()
cP.climb_step()
cP.print_nodes()
cP.climb_step()
cP.print_nodes()

print("Success!")
