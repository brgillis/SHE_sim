#3DCOLORFIT
#NIRAJ WELIKALA (UNIVERSITY OF OXFORD)
#7.12.2015
#Produce models of each type

#type r0_bulge r0_disk z0_disk bulge_compression Npoints nsersic_bulge"


#E type:bulge only
./galsim3d_gen_models_v2 bulge 1 -999 -999 0.64 10000 4
./galsim3d_gen_models_v2 disk -999 1 1 -999 10000 -999
cp 3dmodel_bulge.dat 2dmodel_bulge_e.dat
cp 3dmodel_disk.dat 3dmodel_disk_e.dat

#n=4 bulge + disk

./galsim3d_gen_models_v2 bulge 0.22 -999 -999 0.64 10000 4.00 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999
cp 3dmodel_bulge.dat 2dmodel_bulge_n4.00.dat
cp 3dmodel_disk.dat 3dmodel_disk_n4.00.dat

#n=3.5 bulge + disk

./galsim3d_gen_models_v2 bulge 0.22 -999 -999 0.64 10000 3.50 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999
cp 3dmodel_bulge.dat 2dmodel_bulge_n3.50.dat
cp 3dmodel_disk.dat 3dmodel_disk_n3.50.dat


#n=3 bulge + disk

./galsim3d_gen_models_v2 bulge 0.22 -999 -999 0.64 10000 3.00 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999
cp 3dmodel_bulge.dat 2dmodel_bulge_n3.00.dat
cp 3dmodel_disk.dat 3dmodel_disk_n3.00.dat


#S0 type
./galsim3d_gen_models_v2 bulge 0.22 -999 -999 0.64 10000 2.71 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999
cp 3dmodel_bulge.dat 2dmodel_bulge_n2.71.dat
cp 3dmodel_disk.dat 3dmodel_disk_n2.71.dat


#Sa type
./galsim3d_gen_models_v2 bulge 0.25 -999 -999 0.64 10000 2.56  
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999 
cp 3dmodel_bulge.dat 2dmodel_bulge_n2.56.dat
cp 3dmodel_disk.dat 3dmodel_disk_n2.56.dat

#Sb type
./galsim3d_gen_models_v2 bulge 0.21 -999 -999 0.54 10000 2.00 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999 
cp 3dmodel_bulge.dat 2dmodel_bulge_n2.00.dat
cp 3dmodel_disk.dat 3dmodel_disk_n2.00.dat

#Sc type
./galsim3d_gen_models_v2 bulge 0.22 -999 -999 0.54 10000 1.78  
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999  
cp 3dmodel_bulge.dat 2dmodel_bulge_n1.78.dat
cp 3dmodel_disk.dat 3dmodel_disk_n1.78.dat

#Sd type
./galsim3d_gen_models_v2 bulge 0.24 -999 -999 0.54 10000 1.80 
./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999 
cp 3dmodel_bulge.dat 2dmodel_bulge_n1.80.dat
cp 3dmodel_disk.dat 3dmodel_disk_n1.80.dat



#------------------------
#n=1 pure disk type (this is redundant because the fitting would find this if we had a pure disk with no bulge)

#./galsim3d_gen_models_v2 bulge 0.24 -999 -999 0.54 10000 1.0 
#./galsim3d_gen_models_v2 disk -999 1 0.1667 -999 10000 -999 
#cp 3dmodel_bulge.dat 3dmodel_bulge_n4.dat
#cp 3dmodel_disk.dat 3dmodel_disk_n4.dat
