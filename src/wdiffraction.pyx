
import numpy as np
cimport numpy as np

# get this from np.intc.__name__ but I don't understand how to
# to it "live"
ctypedef np.int32_t CINT


cdef extern void c_choosegrains(
        float *,
        float *,
        float, float, float,
        float *,
        float *,
        float *,
        int, int,
        CINT *, 
        float *
        )

def choosegrains(
    np.ndarray[np.float32_t, ndim=2, mode="c"] XL,
    np.ndarray[np.float32_t, ndim=1, mode="c"] omega,
    float wedge, float chi, float wvln,
    np.ndarray[np.float32_t, ndim=3, mode="c"] UBs,
    np.ndarray[np.float32_t, ndim=3, mode="c"] UBIs,
    np.ndarray[np.float32_t, ndim=2, mode="c"] Ts,
    np.ndarray[CINT, ndim=1, mode="c"]    labels,
    np.ndarray[np.float32_t, ndim=2, mode="c"]    hkls
    ):
    cdef int ngrains
    cdef int npeaks
    # Check array shapes
    assert XL.shape[1] == 3, "XL must be shape (npks,3)"
    npeaks = XL.shape[0]
    assert omega.shape[0] == npeaks, "omega must be shape (npks,)"
    assert UBs.shape[1] == UBs.shape[2] == 3 , "UBs must be (ngrains,3,3)"
    assert UBIs.shape[1] == UBIs.shape[2] == 3 , "UBs must be (ngrains,3,3)"
    ngrains = UBs.shape[0]
    assert UBIs.shape[0] == ngrains, "UBIs, UBs shapes must agree"
    assert Ts.shape[0] ==  ngrains , "Ts shape must be (ngrains, 3)"
    assert Ts.shape[1] ==  3 , "Ts shape must be (ngrains, 3)"
    assert labels.shape[0] ==  npeaks,"Labels.shape must be (npeaks,)"
    assert hkls.shape[0] == npeaks, "hkls.shape must be (npeaks,3 )"
    assert hkls.shape[1] == npeaks, "hkls.shape must be (npeaks,3 )"

    # Call the external

    c_choosegrains( &XL[0,0], &omega[0],
            wedge, chi, wvln,
            &UBs[0,0,0], &UBIs[0,0,0], &Ts[0,0],
            ngrains, npeaks,
            &labels[0],
            &hkls[0,0] )







