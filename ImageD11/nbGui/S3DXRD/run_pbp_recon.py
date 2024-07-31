if __name__ == "__main__":
    # horrible workaround to include id11 code path
    import sys
    
    id11_code_path = sys.argv[1]
    
    sys.path.insert(0, id11_code_path)

    import ImageD11.sinograms.point_by_point
    import ImageD11.sinograms.dataset
    from ImageD11 import cImageD11
    import ast

    dsfile = sys.argv[2]
    hkltol = float(sys.argv[3])
    fpks = float(sys.argv[4])
    dstol = float(sys.argv[5])
    etacut = float(sys.argv[6])
    ifrac = float(sys.argv[7])
    costol = float(sys.argv[8])
    y0 = float(sys.argv[9])
    symmetry = sys.argv[10]
    foridx = [int(x) for x in ast.literal_eval(sys.argv[11])]
    forgen = [int(x) for x in ast.literal_eval(sys.argv[12])]
    uniqcut = float(sys.argv[13])
    minpkint = int(sys.argv[14])

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