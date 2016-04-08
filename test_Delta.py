import matplotlib.pyplot as pl
import numpy as np
import pdb
import compute_Delta as cd

#Redmapper members
redmapmemdir = '/data/ebaxter/redmapper_data/v5.10/'
redmapmemfile = 'dr8_run_redmapper_v5.10_lgt5_catalog_members.fit'
#redmapmemdir = '/data/ebaxter/redmapper_data/v6.3.1/'
#redmapmemfile = 'dr8_run_redmapper_v6.3.1_lgt5_catalog_members.fit'

/data/ebaxter/redmapper_data/v5.10/dr8_run_redmapper_v5.10_lgt5_catalog_members.fit

dr8_run_redmapper_v5.10_lgt5_catalog.fit

#Redmapper data
redmapdir = '/data/ebaxter/redmapper_data/v5.10/'
redmapfile = 'dr8_run_redmapper_v5.10_lgt5_catalog.fit'
#redmapdir = '/data/ebaxter/redmapper_data/v6.3.1/'
#redmapfile = 'dr8_run_redmapper_v6.3.1_lgt5_catalog.fit'

#Compute relevant quantities
output_file = './output/test_Delta_out.fits'
num_lam_bins = 12
num_z_bins = 5
Delta_data = cd.compute_Delta(redmapdir + redmapfile, redmapmemdir + redmapmemfile, output_file, num_lam_bins, num_z_bins)
lam_clusters = Delta_data['lambda']
z_clusters = Delta_data['z']
meanr_clusters = Delta_data['meanr']
avg_meanr_clusters_simple = Delta_data['avg_meanr_simple']
avg_meanr_clusters_gridspline = Delta_data['avg_meanr_gridspline']
avg_meanr_clusters_ungridspline = Delta_data['avg_meanr_ungridspline']
#Compute different versions of Delta
Delta_simple = (meanr_clusters - avg_meanr_clusters_simple)/avg_meanr_clusters_simple
Delta_gridspline = (meanr_clusters - avg_meanr_clusters_gridspline)/avg_meanr_clusters_gridspline
Delta_ungridspline = (meanr_clusters - avg_meanr_clusters_ungridspline)/avg_meanr_clusters_ungridspline

# ############
# Some tests:
# Compare different methods
fig, ax = pl.subplots(1,1)
ax.plot(avg_meanr_clusters_simple, avg_meanr_clusters_ungridspline,'bs')
ax.plot(avg_meanr_clusters_simple, avg_meanr_clusters_gridspline,'rs')
ax.plot(np.array([0.,1.]), np.array([0.,1.]))

fig1, ax1 = pl.subplots(3,1)
ax1[0].scatter(lam_clusters, z_clusters, c = avg_meanr_clusters_simple, label = 'Simple grid', vmin = 0.22, vmax = 0.7)
ax1[1].scatter(lam_clusters, z_clusters, c = avg_meanr_clusters_gridspline, label = 'Spline on grid', vmin = 0.22, vmax = 0.7)
ax1[2].scatter(lam_clusters, z_clusters, c = avg_meanr_clusters_ungridspline, label = 'Ungridded spline', vmin = 0.22, vmax = 0.7)
ax1[0].legend()
ax1[1].legend()
ax1[2].legend()

ax1[0].set_xlabel(r'$\lambda$')
ax1[0].set_ylabel(r'$z$')
ax1[1].set_xlabel(r'$\lambda$')
ax1[1].set_ylabel(r'$z$')
ax1[2].set_xlabel(r'$\lambda$')
ax1[2].set_ylabel(r'$z$')
fig1.tight_layout()

fig2, ax2 = pl.subplots(3,1)
ax2[0].plot(lam_clusters, Delta_simple, 'bs', label = 'Simple grid')
ax2[1].plot(lam_clusters, Delta_gridspline, 'bs', label = 'Spline on grid')
ax2[2].plot(lam_clusters, Delta_ungridspline, 'bs', label = 'Ungridded spline')
ax2[0].legend()
ax2[1].legend()
ax2[2].legend()

ax2[0].set_xlabel(r'$\lambda$')
ax2[0].set_ylabel(r'$\Delta$')
ax2[1].set_xlabel(r'$\lambda$')
ax2[1].set_ylabel(r'$\Delta$')
ax2[2].set_xlabel(r'$\lambda$')
ax2[2].set_ylabel(r'$\Delta$')
fig2.tight_layout()

#fig1.savefig('z_vs_lambda.png')
#fig2.savefig('Delta_vs_lambda.png')

ra_clusters = Delta_data['ra']
dec_clusters = Delta_data['dec']

all_data = np.vstack((ra_clusters, dec_clusters, lam_clusters, z_clusters, Delta_gridspline))

pdb.set_trace()

fname = './output/dr8_redmapper_Delta_6.3.1_lamgt5_12lam_5z_ky2.txt'
np.savetxt(fname, all_data.transpose(), header = "#RA   DEC   Lambda    z    Delta")
test = np.loadtxt(fname)

