


"""
Instructions:


# create the dataset file

python3 -m ImageD11.sinograms.dataset                      \
     /data/visitor/ma5415/id11/20221027                    \
     /data/visitor/ma5415/id11/20221027/analysis/test      \
     NSCOPE_SiMo1000_6                                     \
     DT600nm_Z50                                           \
     dstest_NSCOPE_SiMo1000_6_DT600nm_Z50.h5


# do some segmentation
python3 -m ImageD11.sinograms.lima_segmenter setup dataset.h5
sbatch /data/visitor/ma5415/id11/20221027/analysis/test/NSCOPE_SiMo1000_6/NSCOPE_SiMo1000_6_DT600nm_Z50/slurm/lima_segmenter_slurm.sh

# next step(s)
label the pixels (save that or not?)  
   cost = np.cumsum( dataset.nnz )

save the peaks from each frame?
    s_raw f_raw S_I     float32, float32, float32 == 4 * 3 = 12 bytes  (number of pixels?)
    sc    fc            1 pixel is (uint16, uint16, uint16) = (2 * 3) = 6 bytes
    

determine dataset.npks_per_frame  # keeping all? Keeping Some?

make the graph of all the p



"""

