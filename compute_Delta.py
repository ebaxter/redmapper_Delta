import numpy as np
import sys
import pyfits as pf
from scipy import interpolate

def interp_avg_meanr(lam_list, z_list, lam_bins, z_bins, avg_meanr_mat):
    #Given matrix with average values of R_mem, return appropritate
    #value of mean R_mem for set of clusters
    avg_meanr_list = np.zeros(len(lam_list))
    for ci in xrange(0,len(lam_list)):
        z = z_list[ci]
        lam = lam_list[ci]
        if ((lam > np.max(lam_bins)) or (lam < np.min(lam_bins))):
            raise Exception("Extrapolating in lambda!")
        if ((z > np.max(z_bins)) or (z < np.min(z_bins))):
            raise Exception("Extrapolating in z!")
        match_z_index = np.argmax(z_bins > z)-1
        match_lam_index = np.argmax(lam_bins > lam)-1
        avg_meanr_list[ci] = avg_meanr_mat[match_lam_index, match_z_index]
        if (np.isnan(avg_meanr_list[ci])):
            raise Exception("Mean is not a number!")
    return avg_meanr_list

def compute_Delta(redmap_cluster_file, redmap_member_file, output_file, num_lam_bins, num_z_bins):
    #Compute Delta given redmap catalog and redmapmember catalog
    
    print "Reading redmapper catalog..."
    # Data catalog
    try:
        hdu = pf.open(redmap_cluster_file)
        memmatchid_clusters = hdu[1].data.field('mem_match_id')
        lam_clusters = hdu[1].data.field('lambda_chisq')
        z_clusters = hdu[1].data.field('z_lambda')
        ra_clusters = hdu[1].data.field('ra')
        dec_clusters = hdu[1].data.field('dec')
        hdu.close()
    except:
        print "Invalid redmapper cluster file!"
    num_clusters = len(lam_clusters)
    
    print "Reading redmapper member catalog..."
    # Read in member catalog
    try:
        hdu = pf.open(redmap_member_file)
        memmatchid_members = hdu[1].data.field('mem_match_id')
        pfree_members = hdu[1].data.field('pfree')
        p_members = hdu[1].data.field('p')
        r_members = hdu[1].data.field('r')
        theta_i_members = hdu[1].data.field('theta_i')
        theta_r_members = hdu[1].data.field('theta_r')
        hdu.close()
    except:
        print "Invalid redmapper member file!"
        
    num_members = len(r_members)
    #Compute membership probabilities
    memprob_members = pfree_members*p_members*theta_i_members*theta_r_members

    print "Calculating meanr..."
    # Calculate meanr for every cluster
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
    # Last cluster done separately
    meanr_clusters[counter] =  np.sum(memprob_members[prev_index:]*r_members[prev_index:])/np.sum(memprob_members[prev_index:])
    memmatchid_clusters_output[counter] = memmatchid_members[prev_index]

    # Test for mismatch
    bad = np.where(memmatchid_clusters_output != memmatchid_clusters)[0]
    if (len(bad) > 0):
        raise Exception("Mismatching indices!")

    # Compute average meanr in bins of richness and redshift
    lam_bins_meanr = np.exp(np.linspace(np.log(0.999*np.min(lam_clusters)), np.log(1.001*np.max(lam_clusters)), num = num_lam_bins+1))
    z_bins_meanr= np.linspace(0.999*np.min(z_clusters), 1.001*np.max(z_clusters), num = num_z_bins+1)
    meanr_mat = np.zeros((num_lam_bins, num_z_bins))
    print "Computing meanr across grid..."
    for li in xrange(0,num_lam_bins):
        for zi in xrange(0,num_z_bins):
            in_bin = np.where((lam_clusters > lam_bins_meanr[li]) & (lam_clusters <= lam_bins_meanr[li+1]) & (z_clusters > z_bins_meanr[zi]) & (z_clusters <= z_bins_meanr[zi+1]))[0]
            #Only compute mean if there are clusters in bin
            if (len(in_bin) > 0):
                meanr_mat[li,zi] = np.mean(meanr_clusters[in_bin])

    print "Computing average value of meanr as a function of richness and z..."
    # Get average value of meanr for clusters of that richness and redshift, evaluated at richness and redshift of all clusters

    # Use simple interpolation of grid
    avg_meanr_clusters_simple = interp_avg_meanr(lam_clusters, z_clusters, lam_bins_meanr, z_bins_meanr, meanr_mat)

    # Get avg meanr using ungridded spline fit
    spline = interpolate.SmoothBivariateSpline(lam_clusters, z_clusters, meanr_clusters)
    avg_meanr_clusters_ungridspline = spline.ev(lam_clusters, z_clusters)

    # Get avg meanr using gridded spline fit
    meanr_mat = np.nan_to_num(meanr_mat)
    rbs = interpolate.RectBivariateSpline(np.exp(0.5*(np.log(lam_bins_meanr[1:]) + np.log(lam_bins_meanr[:-1]))), 0.5*(z_bins_meanr[1:] + z_bins_meanr[:-1]), meanr_mat, ky = 2)
    avg_meanr_clusters_gridspline = rbs.ev(lam_clusters, z_clusters)

    #Compute Delta for all clusters
    avg_meanr_clusters = np.copy(avg_meanr_clusters_gridspline)
    Delta_clusters = (meanr_clusters - avg_meanr_clusters)/avg_meanr_clusters

    print "Outputting to fits file..."
    #Create columns for fits file
    ra_col = pf.Column(name='RA', format='E', array=ra_clusters)
    dec_col = pf.Column(name='DEC', format='E', array=dec_clusters)
    lam_col = pf.Column(name='lambda_chisq', format='E', array=lam_clusters)
    z_col = pf.Column(name='z_lambda', format='E', array = z_clusters)
    meanr_col = pf.Column(name='R_mem', format='E', array=meanr_clusters)
    avg_meanr_col = pf.Column(name='mean_R_mem', format='E', array=avg_meanr_clusters)
    Delta_col = pf.Column(name='Delta', format='E', array=Delta_clusters)
    #All columns
    cols = pf.ColDefs([ra_col, dec_col, lam_col, z_col, meanr_col, avg_meanr_col, Delta_col])
    #Write to fits
    tbhdu = pf.TableHDU.from_columns(cols)
    try:
        tbhdu.writeto(output_file)
    except:
        print "Unable to write to output file!"
        raise

if __name__ == "__main__":
    if (len(sys.argv) != 6):
        print "Proper usage: "
        print "python compute_Delta.py cluster_file member_file output_file num_lam_bins num_z_bins"
    else:
        num_lam_bins = 0
        num_z_bins = 0
        try:
            num_lam_bins = int(sys.argv[4])
            num_z_bins = int(sys.argv[5])
        except:
            print "Number of lambda, z bins not properly specified!"
        if (num_lam_bins*num_z_bins > 0):
            compute_Delta(sys.argv[1], sys.argv[2], sys.argv[3], num_lam_bins, num_z_bins)
