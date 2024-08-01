if __name__ == "__main__":

    import sys
    import ImageD11.sinograms.point_by_point
    import ImageD11.sinograms.dataset
    from ImageD11 import cImageD11
    import ast

    dsfile = sys.argv[1]
    hkltol = float(sys.argv[2])
    fpks = float(sys.argv[3])
    dstol = float(sys.argv[4])
    etacut = float(sys.argv[5])
    ifrac = float(sys.argv[6])
    costol = float(sys.argv[7])
    y0 = float(sys.argv[8])
    symmetry = sys.argv[9]
    foridx = [int(x) for x in ast.literal_eval(sys.argv[10])]
    forgen = [int(x) for x in ast.literal_eval(sys.argv[11])]
    uniqcut = float(sys.argv[12])
    minpkint = int(sys.argv[13])

    print('Loading dset')
    ds = ImageD11.sinograms.dataset.load(dsfile)

    print('Loading peaks')
    ImageD11.cImageD11.cimaged11_omp_set_num_threads(ImageD11.cImageD11.cores_available())
    cf_2d = ds.get_cf_2d()
    ImageD11.cImageD11.cimaged11_omp_set_num_threads(1)
    cf_2d.filter(cf_2d.Number_of_pixels > minpkint)

    print('Making pbp object')
    pbp_object = ImageD11.sinograms.point_by_point.PBP(ds.parfile,
                                                       ds,
                                                       hkl_tol=hkltol,
                                                       fpks=fpks,
                                                       ds_tol=dstol,
                                                       etacut=etacut,
                                                       ifrac=ifrac,
                                                       cosine_tol=costol,
                                                       y0=y0,
                                                       symmetry=symmetry,
                                                       foridx=foridx,
                                                       forgen=forgen,
                                                       uniqcut=uniqcut)

    pbp_object.setpeaks(cf_2d)

    print('Go for pbp')
    pbp_object.point_by_point(ds.pbpfile, loglevel=3)
