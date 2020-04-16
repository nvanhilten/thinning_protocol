/*
This code was written and implemented by Niek van Hilten and Herre Jelger Risselada for GROMACS 2019.3.
When using this code, please cite:
-----------------
Van Hilten, N., Stroh, K.S. and Risselada, H.J., (2020). Membrane Thinning Induces Sorting of Lipids and the Amphipathic Lipid Packing Sensor (ALPS) Protein Motif. Front. Physiol. 11:250. doi:10.3389/fphys.2020.00250 -----------------



For the original GROMACS paper, please see:
-----------------
Abraham, M. J., Murtola, T., Schulz, R., Pall, S., Smith, J. C., Hess, B., et al. (2015).
GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. Soft. X 1-2, 19–25
-----------------



This protocol can be used to induce thinning of lipid membranes in MD simulations. It should be called from the .mdp file, as demonstrated in the THINNING_EXAMPLE.mdp file.
Also, it is required to define a index group "THINNING" which contains the atoms on which the inward directed force will act. 
We recommend adding all tail and glycerol atoms of all lipid molecules to this "THINNING" index group.
*/



/*insert the following lines at line 147 in /src/gromacs/mdrun/md.cpp*/
void do_thinning(tensor vir, rvec *x, rvec *f, int START, int HOMENR, t_commrec *cr, matrix box, unsigned short *cU1, real thickness, real kforce, real xzone, real lshift)
{
    int n,d,j,m;
    rvec CENT,dr,df;
    real lengthx,scale,hthick,dz=0;

    /* half the box*/
    for(d=0;d<DIM;d++)
        CENT[d]=0.5*box[d][d];

    /*loop over atoms on each node, note that atoms are split over the nodes*/

    hthick = 0.5*thickness;

    for(n=START;n<START+HOMENR;n++)
    {
        if(cU1[n] == 0)
        {
            clear_rvec(df);
            clear_rvec(dr);
            rvec_sub(CENT,x[n],dr);

            lengthx = 2.0*fabs(CENT[XX] - x[n][XX]);

            /* define zones in x dimension*/

            /*thin zone, buffer zone, and normal zone*/
            if(lengthx <= xzone)
                scale = 1.0;
            else if(lengthx xzone + lshift)
                scale = ((xzone + lshift) - lengthx)/lshift;
            else 
                scale=0;

            dz = fabs(dr[ZZ]) - hthick;

            if(dz > 0) 
            {
                if(dr[ZZ] > 0)
                    df[ZZ] = scale*kforce*dz;

                if(dr[ZZ] < 0)
                    df[ZZ] = -1.0*scale*kforce*dz;

                f[n][ZZ] += df[ZZ];
                        
                /* virial update... only for pressure calculation*/
                vir[ZZ][ZZ] -= -0.5*df[ZZ]*dz;
            }
        }
    }
}





/*The function is called by inserting the following lines at line 938 in /src/gromacs/mdrun/md.cpp*/
            if(ir->userint1 != 0)
            {
                do_thinning(force_vir, 
                            as_rvec_array(state->x.arrayRefWithPadding().unpaddedArrayRef().data()), 
                            as_rvec_array(f.arrayRefWithPadding().unpaddedArrayRef().data()), 
                            0, mdatoms->homenr, cr, state->box, mdatoms->cU1,
                            ir->userreal1, ir->userreal2, ir->userreal3, ir->userreal4);
            }

