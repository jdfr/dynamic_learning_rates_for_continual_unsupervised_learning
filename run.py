import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.io import savemat, loadmat
import argparse
from sklearn.cluster import MiniBatchKMeans
#import code; code.interact(local=vars())


def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--num_clusters', type=int)
  parser.add_argument('--num_repeats', type=int)
  parser.add_argument('--num_meta_repeats', type=int)
  parser.add_argument('--input', type=str)
  parser.add_argument('--output', type=str)
  args = parser.parse_args()
  return args

if __name__ == "__main__":
  #args = Args()
  args = get_args()
  d = loadmat(args.input)
  batches_expanded = np.squeeze(d['batches_expanded'])
  num_batches = batches_expanded.size
  all_centroids = np.empty((num_batches,1), dtype=object)
  mbkm = MiniBatchKMeans(n_clusters=args.num_clusters)
  for k in range(num_batches):
    batches = batches_expanded[k].T
    num     = int(np.round(batches.shape[0]/args.num_repeats))
    low     = 0
    high    = num
    for k1 in range(args.num_meta_repeats):
      for k2 in range(args.num_repeats):
        batch = batches[low:high,:]
        #print(k, k1, k2, batch.shape)
        #import code; code.interact(local=vars())
        mbkm.partial_fit(batch)
        low  += num
        high += num
      if k1!=(args.num_meta_repeats-1):
        batches = batches[np.random.permutation(batches.shape[0]),:]
        low     = 0
        high    = num
    #import code; code.interact(local=vars())
    all_centroids[k,0] = mbkm.cluster_centers_
  savemat(args.output, {'all_centroids': all_centroids})
