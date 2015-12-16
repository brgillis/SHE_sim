#..............................................................................
#.
#.                Makefile to link with the CCfits, cfitsio and gsl libraries
#.
#.      
#.     usage:    make image=
#.
#..............................................................................

#
# -- compiler options
#

CCFLAGS = -L/opt/local/lib -I/opt/local/include -lnlopt -lgsl -lgslcblas -lfftw3 ${INCLUDE} -lm ${LDLIBS}  
FITSIO_INCDIR =/Users/n_welikala/software/cfitsio/include

INCLUDE=-I${FITSIO_INCDIR} 
LDLIBS = -L/Users/n_welikala/software/cfitsio/lib -lcfitsio 


image=galsim3d_gen_models
#image=fit_goods_bvizyjh


exten   = c

# -- libraries
#LIBS_G = -L$(NR_DIR) \
	-lrecipes_c

# -- Actions
$(image)_always:

	gcc -o $(image) $(image).$(exten)  $(CCFLAGS) 

