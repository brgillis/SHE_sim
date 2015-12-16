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

CCFLAGS = -g3 -L/home/brg/lib -I/home/brg/include -lgsl -lgslcblas -lfftw3 ${INCLUDE} -lm ${LDLIBS}  
FITSIO_INCDIR =/home/brg/include/cfitsio

INCLUDE=-I${FITSIO_INCDIR} 
LDLIBS = -L/home/brg/lib/cfitsio -lcfitsio 


image=galsim3d_gen_models
#image=fit_goods_bvizyjh


exten   = c

# -- libraries
#LIBS_G = -L$(NR_DIR) \
	-lrecipes_c

# -- Actions
$(image)_always:

	gcc -o $(image) $(image).$(exten)  $(CCFLAGS) 

