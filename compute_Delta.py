from math import pi
import numpy as np
import matplotlib.pyplot as pl
import pdb
import pyfits as pf
import pickle as pk
from scipy import interpolate

def interp_avg_meanr(lam_list, z_list, lam_bins, z_bins, avg_meanr_mat):
    avg_meanr_list = np.zeros(len(lam_list))
    for ci in xrange(0,len(lam_list)):
        z = z_list[ci]
        lam = lam_list[ci]
        if ((lam > np.max(lam_bins)) or (lam < np.min(lam_bins))):
            print "Extrapolating in lambda!"
            pdb.set_trace()
        if ((z > np.max(z_bins)) or (z < np.min(z_bins))):
            print "Extrapolating in z!"
            pdb.set_trace()
        match_z_index = np.argmax(z_bins > z)-1
        match_lam_index = np.argmax(lam_bins > lam)-1
        avg_meanr_list[ci] = avg_meanr_mat[match_lam_index, match_z_index]
        if (np.isnan(avg_meanr_list[ci])):
            print "Mean is not a number!"
            pdb.set_trace()
    return avg_meanr_list


#Calculate meanr for all clusters in catalog file and output to new file
#also compute x = (meanr - <meanr>)/<meanr>

outdir = '/home/ebaxter/DES/redmapper_angular_clustering/DeltaCode/output/'
outfile = 'dr8_run_redmapper_v5.10_lgt5_catalog_meanr_withDeltagridspline_v01x.pk'

#Redmapper members
redmapmemdir = '/data/ebaxter/redmapper_data/v5.10/'
redmapmemfile = 'dr8_run_redmapper_v5.10_lgt5_catalog_members.fit'

#Redmapper data
redmapdir = '/data/ebaxter/redmapper_data/v5.10/'
redmapfile = 'dr8_run_redmapper_v5.10_lgt5_catalog.fit'

print "Reading redmapper catalog..."
# Data catalog
hdu = pf.open(redmapdir + redmapfile)
memmatchid_clusters = hdu[1].data.field('mem_match_id')
lam_clusters = hdu[1].data.field('lambda_chisq')
z_clusters = hdu[1].data.field('z_lambda')
ra_clusters = hdu[1].data.field('ra')
dec_clusters = hdu[1].data.field('dec')
hdu.close()
num_clusters = len(lam_clusters)

print "Reading redmapper member catalog..."
#Read in member catalog
hdu = pf.open(redmapmemdir + redmapmemfile)
memmatchid_members = hdu[1].data.field('mem_match_id')
memprob_members = hdu[1].data.field('pfree')*hdu[1].data.field('p')
r_members = hdu[1].data.field('r')
hdu.close()
num_members = len(r_members)

print "Calculating meanr..."
#Calculate meanr for every cluster
# This stores <R_mem> for every cluster
meanr_clusters = np.zeros(num_clusters)-1.
memmatchid_clusters_output = np.zeros(num_clusters)-1.
counter = 0
prev_index = 0
for mi in xrange(1,num_members):
    if (mi % 1000000 == 0):
        print "member = ", mi+1, " out of ", num_members+1
    if (memmatchid_members[mi] != memmatchid_members[prev_index]):
       meanr_clusters[counter] =  np.sum(memprob_members[prev_index:mi]*r_members[prev_index:mi])/np.sum(memprob_members[prev_index:mi])
       memmatchid_clusters_output[counter] = memmatchid_members[prev_index]
       counter += 1
       prev_index = mi
#Last cluster done separately
meanr_clusters[counter] =  np.sum(memprob_members[prev_index:]*r_members[prev_index:])/np.sum(memprob_members[prev_index:])
memmatchid_clusters_output[counter] = memmatchid_members[prev_index]

#Test for mismatch
bad = np.where(memmatchid_clusters_output != memmatchid_clusters)[0]
if (len(bad) > 0):
    print "Mismatching indices!"
    pdb.set_trace()

#Compute average meanr in bins of richness and redshift
num_lam_bins = 10
num_z_bins = 5
lam_bins_meanr = np.exp(np.linspace(np.log(0.999*np.min(lam_clusters)), np.log(1.001*np.max(lam_clusters)), num = num_lam_bins+1))
z_bins_meanr= np.linspace(0.999*np.min(z_clusters), 1.001*np.max(z_clusters), num = num_z_bins+1)
meanr_mat = np.zeros((num_lam_bins, num_z_bins))
print "Computing meanr across grid..."
for li in xrange(0,num_lam_bins):
    for zi in xrange(0,num_z_bins):
        in_bin = np.where((lam_clusters > lam_bins_meanr[li]) & (lam_clusters <= lam_bins_meanr[li+1]) & (z_clusters > z_bins_meanr[zi]) & (z_clusters <= z_bins_meanr[zi+1]))[0]
        meanr_mat[li,zi] = np.mean(meanr_clusters[in_bin])

print "Computing average value of meanr as a function of richness and z..."
#Get average value of meanr for clusters of that richness and redshift, evaluated at richness and redshift of all clusters

#Use simple interpolation of grid
avg_meanr_clusters_simple = interp_avg_meanr(lam_clusters, z_clusters, lam_bins_meanr, z_bins_meanr, meanr_mat)

#Get avg meanr using ungridded spline fit
spline = interpolate.SmoothBivariateSpline(lam_clusters, z_clusters, meanr_clusters)
avg_meanr_clusters_ungridspline = spline.ev(lam_clusters, z_clusters)

#Get avg meanr using gridded spline fit
meanr_mat = np.nan_to_num(meanr_mat)
rbs = interpolate.RectBivariateSpline(0.5*(lam_bins_meanr[1:] + lam_bins_meanr[:-1]), 0.5*(z_bins_meanr[1:] + z_bins_meanr[:-1]), meanr_mat)
avg_meanr_clusters_gridspline = rbs.ev(lam_clusters, z_clusters)

#############
#Some tests:
#Compare different methods
fig, ax = pl.subplots(1,1)
ax.plot(avg_meanr_clusters_simple, avg_meanr_clusters_ungridspline,'bs')
ax.plot(avg_meanr_clusters_simple, avg_meanr_clusters_gridspline,'rs')
ax.plot(np.array([0.,1.]), np.array([0.,1.]))

fig_lz, ax_lz = pl.subplots(1,1)
ax_lz.scatter(lam_clusters, z_clusters, c = avg_meanr_clusters_simple)


###########

pdb.set_trace()

avg_meanr_clusters = np.copy(avg_meanr_clusters_gridspline)
Delta_clusters = (meanr_clusters - avg_meanr_clusters)/avg_meanr_clusters

pdb.set_trace()

pickle_data = {"redmapfile":redmapfile, "redmapmemfile":redmapmemfile, "mem_matchid":all_memmatchid,'meanr':all_meanr, 'avg_meanr':all_avg_meanr, 'Delta':all_Delta, 'ra':orig_ra, 'dec':orig_dec, 'lambda':orig_lam}
pk.dump(pickle_data, open('meanr_data_v5.10_lamgt5_withcoords.pk','wb'))


pdb.set_trace()
