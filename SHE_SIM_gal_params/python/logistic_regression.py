""" This is David Harvey's code for predicting Sersic indices from stellar mass.

"""

import sklearn.linear_model as linear_model 
import numpy as np
import read_cosmos as rc
import random_forest as rf
import ipdb as pdb
from matplotlib import mlab as mlab
from matplotlib import pyplot as plt
from scipy.stats import norm
from matplotlib import gridspec


def sersic_index():
    '''
    Try and figure out functional forms using sersic_index

    '''
    
    cosmos_cat = rc.main()

    keyword='BULG'
    predictions = fit_logistic_reg( cosmos_cat, keyword, norder=1,
                                  nbins=40, nsamples=10000)


def physical_size():
    '''
    RADIUS_KPC

    '''

    cosmos_cat = rc.main()
    keyword='RADIUS_KPC'
    predictions = fit_logistic_reg( cosmos_cat, keyword )

def fit_logistic_reg(catalogue, keyword, nsamples=20000,
                     plot_name='linear_regression.pdf',
                     nbins=20, norder=2, features=None ):
    '''
    Using the sklearn logistical regression functions have a go at fitting
    functions


    '''
    

   
    training_features, training_answers, testing_samples, testing_answers, used_features \
        = rf.clean_sample( catalogue, keyword, nsamples=nsamples, features=features )
    
    #in order to do polynomial linear regression just use basis in orders, very easy

    new_training, new_features = PolynomialFeatures(used_features, training_features,\
                                                    update_names=True, norder=norder)

    new_testing, new_features = PolynomialFeatures(used_features, testing_samples,\
                                                   update_names=False, norder=norder)
          
    regression_obj = linear_model.LogisticRegression( )
    print 'Fitting'
    regression_fit = regression_obj.fit( new_training, training_answers)
    print 'Fitted'
    
    predictions = regression_fit.predict( new_testing )
    
    rf.plot_predictions( testing_answers, predictions,
                      new_features, regression_fit.coef_[0],
                      keyword, filename=plot_name, nbins=nbins)
    
    print 'With a score of ', regression_fit.score( new_testing, testing_answers)

    return regression_fit.intercept_, regression_fit.coef_, \
      new_testing, testing_answers

    


    
    
def PolynomialFeatures( names, sample, update_names=True, norder=2 ):
    '''
    Create more featues upto norder

    update_names : update the list of feeatures being used

    '''
    nNames = len(names)

    if not update_names:
        nNames /= norder
    
    new_sample = np.zeros( (sample.shape[0], sample.shape[1]*norder ))
    new_sample[ :, :nNames] = sample
    for iName in xrange(nNames):
        for iOrder in xrange( 2, norder+1 ):

            powered_sample = sample[:, iName]**(iOrder)
            
            new_sample[ :, nNames + iName - norder] =  powered_sample
            if update_names:
                names.append( names[iName]+'_'+str(iOrder) )

    return new_sample, names

def type_distributions(catalogue=None):
    '''
    
    Determine what the sersic index distributions are for
    each bulge type (which i can determine pretty well)

    This code will fit Gaussian distirbutions each distirbutions
    from the 
    

    '''
    if catalogue is None:
        catalogue = rc.main()
    nsamples= 2000
    keyword='MAG_AUTO'
    
    training_features, training_answers, testing_samples,\
        testing_answers, used_features \
        = rf.clean_sample( catalogue, keyword, nsamples=nsamples, features=None )


    nbins = 20

    bulge_values = np.array([0, 1, 2, 3, 9])
    nbulges= len(bulge_values)
    binned_dist = np.zeros( (nbulges, 2), float)
    sersic_col = np.arange(len(used_features))\
      [np.array(used_features) == 'SERSIC_N'][0]
    bulge_col = np.arange(len(used_features))\
      [np.array(used_features) == 'BULG'][0]
    for iBulge in xrange(nbulges):
        iBulge_val = bulge_values[iBulge]
        iX, iMean, iWidth, iPDF = \
            fit_gauss( \
                training_features[ training_features[:, bulge_col] == iBulge_val,
                                   sersic_col], nbins=nbins )
        binned_dist[iBulge, 0] = iMean
        binned_dist[iBulge, 1] = iWidth
        #plt.hist(  training_features[ training_features[:,4] == iBulge_val, 3],\
         #          bins=nbins)
        plt.plot( iX, iPDF,'-', label='Bulge value '+str(iBulge_val))
    plt.xlabel('Sersic index N')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig('../plots/sersic_distribiutions')   
    plt.show()
    return binned_dist



def fit_gauss( data, nbins=20 ):
    '''

    Fit a gaussian to data and return x and u

    '''

    bulge = np.histogram( data, bins=20)
    mean, width =  norm.fit( data )
    x = (bulge[1][:-1] + bulge[1][1:])/2.
    pdf = norm.pdf(x, mean, width)
    print mean, width
    return x, mean, width, pdf

def estimate_sersic():
    '''
    Using the logistic regerssion i will estimate the
    bulge fraction, then from this I will
    estimate the sersic index
    '''

    cosmos_cat = rc.main()

    keyword='BULG'
    intercept, coeffs, new_testing, answers = \
        fit_logistic_reg( cosmos_cat, keyword, norder=1,
                           nbins=40, nsamples=20000,
                           features=['MASS','REDSHIFT','BULG'])

    print coeffs, intercept
    sersic_dists = type_distributions( catalogue=cosmos_cat )
    
    my_predict = np.zeros(len(new_testing[:,0]))
    sersics = np.zeros(len(new_testing[:,0]))
    
    for iGal in  xrange(len(new_testing[:,0])):
        my_predict[iGal] = \
          determine_class( new_testing[iGal,:], \
                           coeffs, \
                           intercept,
                           answers)


    bulge_vals = np.unique( my_predict )
    for iBulge in xrange(len(bulge_vals)):

        iPredict = my_predict[ my_predict==bulge_vals[iBulge]]

        iSersics = np.random.normal( sersic_dists[iBulge,0], 
                                     sersic_dists[iBulge,1],
                                     len(iPredict) )

        sersics[ my_predict==bulge_vals[iBulge]] = iSersics

    bins=np.linspace(-5,15,60)

    gs = gridspec.GridSpec( 1,2)

    ax1 = plt.subplot(gs[0,0])
    correlation, x, y = np.histogram2d( sersics, new_testing[:,-1], bins=bins)
    ax1.imshow(correlation.T,origin='lower'\
               ,extent=[x[0],x[-1],y[0],y[-1]], \
               interpolation='none',vmax=5000)
    ax1.plot([0,10],[0,10],'k-')
    ax1.set_xlim(0,5)
    ax1.set_ylim(0,5)

    ax2 = plt.subplot(gs[0,1])
    
    ys, xs = np.histogram(sersics, bins=bins)
    xs_centre = (xs[:-1] + xs[1:])/2.
    ya, xa = np.histogram(new_testing[:,-1], bins=bins)
    xa_centre = (xa[:-1] + xa[1:])/2.

    
    ax2.plot(xs_centre, ys, '-', label='Predicted')
    ax2.plot(xa_centre, ya, '-', label='True')
    ax2.legend()
    plt.savefig('../plots/sersic_n.pdf', format='pdf')
    
    plt.show()


    return sersics, new_testing


#EVERYTHING BELOW HERE IS WHAT YOU NEED BRIAN


def euclid_estimate_sersic():
    '''
    So the estimate_sersic function will fit the logistic regression
    and then fit the gaussians and then estimate the sersic n
    which is not useful for euclid as you dont want to have the
    cosmos catalogue.

    So this is the code using the arleady fitted coeffocoemts
    and gaussians

    
    '''
    #First get some fake stellar masses and redshifts
    #and some answers to compare to
    #This is where you will input your own data
    #It will need to be in the form 2D array
    # data = Ngal x 2 where the columns go MASS, REDSHIFT
    data, answers = get_test_data()
    
    #Then get the fitted intercepts and coefficients
    intercept, coeffs = logistic_coefficients()
    
    #Set up the array which I will determine the bulge values
    #determine_class could be done much more efficiently
    #Than looping over each galaxy
    bulge_vals = np.zeros(len(answers))

    #What are the different class options
    classifications = np.unique( answers )
    
    for iGal in  xrange(len(answers)):
        bulge_vals[iGal] = \
          determine_class( data[iGal,:], \
                           coeffs, \
                           intercept,
                           classifications )

    #Taking the bulge values estimate the sersic indexes
    sersics = determine_sersics( bulge_vals )
    
    return sersics, answers


def determine_class( params, coeffs, intercept, classifications):
    '''
    Determine the classification from logistic regression
    params are the params of a galaxy

    The way this works is using the coefficients derived in the fit
    you work out the probability of each classification and then select the
    class which has the highest probability.

    The probability of a given classification is given by

    p(classification) = ( 1 + exp( -T ) )**-1
    where
    T = SUM(  COEFFS * PARAMS ) + INTERCEPT

    so here the input is the params of the galaxy e.g. MASS, REDSHIFT
    coeffs is a matrix of N_CLASSIFICATIONS x N_PARAMS

    This code loops through the number of choices and works out the prob and
    returns the classification....could be done better
    
    '''

    iProb = np.zeros(len(coeffs[:,0]))
    for i in xrange(len(coeffs[:,0])):
        
        iT =np.sum( params*coeffs[i,:] ) + intercept[i]
        iProb[i] = 1./(1+np.exp(-1.*iT))
        
    return classifications[np.argmax( iProb)]

def logistic_coefficients():
    '''
    Return the fitted coeffcients from the logistical regression
    This is based on estaimating the COSMOS bulge fraction
    using just stellar mass and redshift

    Also the y_interept of the equation

    Because there are five bulge classifications
    0, 1, 2, 3, 9
    And the parameters
    stellar mass and redshift (in that order)
    It will return an array 5x2
    '''
    coeffs = np.array([[ 1.49695144, -1.9082208 ],
                        [ 0.58817801, -0.87745313],
                        [-0.40911011, -0.36022461],
                        [-0.57818587,  0.78895037],
                        [ 0.30539142,  0.71603493]])
    y_intercept = np.array([-16.49570234,  -6.52302501,
                         3.32521522,   4.11680201,  -5.24419872])

    return y_intercept, coeffs

def bulge_sersic_distributions():
    '''
    So this gives the values for the fitted Gaussians
    to the distirbutions of sersics for each bulge
    ratio. Each bulge value correspsonds to a given
    Gaussian distributon of sersic indexes.

    Given that there are 5 bulge values, and each Gaussian
    is parameterised by mean and std I return an array
    of 5 x 2
    
    '''
    return np.array([[ 3.29313486,  1.23736668],
                    [ 1.90371314,  0.85641805],
                    [ 1.03420396,  0.4113925 ],
                    [ 0.64977613,  0.25539489],
                    [ 2.3155855 ,  1.84302273]])

def get_test_data():
    '''
    To test this return an array of stellar mass and redshift
    and the corresponding sersic indexes

    training_features will be an array of nsamples x 2
    training_answres will be an array of nsamples x 1
    '''

    cosmos_cat = rc.main()
    nsamples = 20000
    features = ['MASS', 'REDSHIFT']
    keyword = 'SERSIC_N'
    
    training_features, training_answers, testing_samples,\
        testing_answers, used_features \
        = rf.clean_sample( cosmos_cat, keyword, nsamples=nsamples,
                           features=features )

    return training_features, training_answers
        
    

    
    
def determine_sersics( bulge_vals):
    '''
    Determine the sersic indexes from the given bulge values
    It will get the gaussian distributions from the fitting early
    and draw randomly from these distributions.
    
    '''
    #Get the number of bulge classifications
    uniq_bulge_vals = np.unique( bulge_vals )
    sersic_dists =  bulge_sersic_distributions()
    sersics = np.zeros(len(bulge_vals))
    for iBulge in xrange(len(uniq_bulge_vals)):

        iPredict = bulge_vals[ bulge_vals==uniq_bulge_vals[iBulge]]

        iSersics = np.random.normal( sersic_dists[iBulge,0], 
                                     sersic_dists[iBulge,1],
                                     len(iPredict) )

        sersics[ bulge_vals==uniq_bulge_vals[iBulge]] = iSersics
        
    return sersics
