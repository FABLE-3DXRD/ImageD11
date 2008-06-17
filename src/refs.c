


void get_key( 
        int h,
        int k,
        int l,
        float *sym,
        int nsym,
        int *kh,
        int *kk,
        int *kl){
    /* Based on George Sheldrick's fortran code from the Sienna computing
     * school
     *
     * h, k, l are the integer hkl indices of a diffraction spot
     *
     * sym[] holds 12 floating point numbers per symmetry matrix
     * nysm is the number of symmetry matrices
     *
     * kh, kk, kl will hold the hkl indices transformed into a unique a.s.u
     *
     * All credit to George, all bugs from me.
     */
    int o, s, a, b, c ;
    /* Initialise as current hkl */
    *kh = h*1.0f;    *kk = k*1.0f;     *kl = l*1.0f;
    for( s=-1; s < 2 ; s = s + 2 ) { /* sign in -1, 1 */
        for( o = 0 ; o < nsym ; o++){ /* Foreach symmop */
            a = s * (sym[12*o+0]*h + sym[12*o+3]*k + sym[12*o+6]*l);
            b = s * (sym[12*o+1]*h + sym[12*o+4]*k + sym[12*o+7]*l);
            c = s * (sym[12*o+2]*h + sym[12*o+5]*k + sym[12*o+8]*l);
            if(c < *kl) continue;
            if(c > *kl) { *kh=a; *kk=b; *kl=c; continue; }
            if(b < *kk) continue;
            if(b > *kk) { *kh=a; *kk=b; *kl=c; continue; }
            if(a <= *kh) continue;
            if(a > *kh) { *kh=a; *kk=b; *kl=c; continue; }
        }
    }
    /* At this point we have checked all symmetry equivalent hkls and found
     * and stored the one with max hkl in the return args
     */
}
        
void get_many_keys( int * hkl, int * keys, int nref, float * sym, int nsym){
    /* Apply the get_key code to a (flat) array of hkls
     */
    int i;
    for(i=0; i<nref; i++){
        get_key(  hkl[3*i],   hkl[3*i+1],   hkl[3*i+2],
                 sym, nsym,
                 &keys[3*i], &keys[3*i+1], &keys[3*i+2]);
                }
}


