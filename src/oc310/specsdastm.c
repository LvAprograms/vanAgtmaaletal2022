/* === Output data every time step during sequenatial data assimilation === */
// Copied and cleaned from postproc()
void outputsda_lab()
{
	FILE *fldG;
	FILE *fldG2;
	FILE *fldM;
	FILE *fldM2;
	FILE *fldD;

	if (printmod) fprintf(fp_log,"Saving data each time step to visualization during data assimilation ...\n"); fflush(fp_log);

	count1++;

	if (count1 == 1)
	{
		if (restart==0)
		{       
			// Need like this, since in we are going to append data in loop  below
			fldM = fopen(fileTxtMarkerOutput,"wt");
			for (int im=0;im<nr_physmarkers;im++)
				{fprintf(fldM,"%d ", mphystrack[im]);}
			fprintf(fldM,"\n");
			fclose(fldM);
						
			fldG = fopen(fileTxtOutputGPS,"wt");
			for (int im=0;im<nr_gpsmarkers;im++)
				{fprintf(fldG,"%d ", mgpstrack[im]);}
			fprintf(fldG,"\n");
			fclose(fldG);
		}
	}
			
	// --- Marker track data ---
	fldM2 = fopen(fileTxtMarkerOutput,"a");
	for (int imt=0; imt<nr_physmarkers; imt++)
		{ fprintf(fldM2,"%e %.13e %.13e %.8e %.8e %e %e %e %e %e %e \n", timesum, markx[mphystrack[imt]], marky[mphystrack[imt]], markvx[mphystrack[imt]], markvy[mphystrack[imt]], markp[mphystrack[imt]], msbrit[mphystrack[imt]], markxx[mphystrack[imt]],markxy[mphystrack[imt]],markexx[mphystrack[imt]],markexy[mphystrack[imt]]); }
	fclose(fldM2);

	// --- GPS marker track data ---
	// Added vx and vy here, so could also potentially use as data input
	// Add also when reading in MATLAB!
	fldG2 = fopen(fileTxtOutputGPS,"a");
	for (int imt=0; imt<nr_gpsmarkers; imt++)
		{ fprintf(fldG2,"%e %.13e %.13e %e %e %e %e %e %e %e \n", timesum, markx[mgpstrack[imt]], marky[mgpstrack[imt]], markvx[mgpstrack[imt]], markvy[mgpstrack[imt]], markp[mgpstrack[imt]], markxx[mgpstrack[imt]],markxy[mgpstrack[imt]],markexx[mgpstrack[imt]],markexy[mgpstrack[imt]]); }
	fclose(fldG2);
	
	// --- Data in analogue borehole nodes ---
	// Location: same as nodes loaded in matlab from hdf5
	// ATAT Add command line passing or smt, so can not forget to update
	int nx_borehole = 269;
	int ny_borehole = 70;
	int nn_borehole = (nx_borehole-1)*ynumy+ny_borehole;
	
	// Print to file
	fldD = fopen(fileTxtOutputData,"a");
	fprintf(fldD,"%e %e %e %e %e %e \n", timesum, vx[nn_borehole], vy[nn_borehole], sxx[nn_borehole], sxy[nn_borehole], pr[nn_borehole]);
	fclose(fldD);
	
}


/* === Load changes from sequenatial data assimilation algorithm and interpolate them to markers === */
void inputfromsda_lab(int assimstep)
{
	char fileSV[100];
	
	// --- Load nodal point values from .txt file ---

	sprintf(fileSV,"posterior_statevector_%s_step%03d.txt",exp_name,assimstep);
	fl = fopen(fileSV,"rt");
	
	// Get info for storage of state variables
	ffscanf(); int nrStates = atoi(sa);
	ffscanf(); int typeofstates =atoi(sa);
	
	// Load each state variable from a column that has same sequence as in here (for each x, for each y)
	int isv;
	// Load matlab node (1) into C-array element (0), as 1 difference at array start
	for (isv=0;isv<nrStates/typeofstates;isv++)
		{ffscanf(); vx[isv]=atof(sa);}
	for (isv=0;isv<nrStates/typeofstates;isv++)
		{ffscanf(); vy[isv]=atof(sa);}
	// sxx and sxy are updated stresses from last time step, which will be used as old stresses in viscalc latter. sxx/ye here are only to interpolate only stress change from nodes to markers, but not done now, since do not need to do it often
	for (isv=0;isv<nrStates/typeofstates;isv++)
		{ffscanf(); sxx[isv]=atof(sa);}
	for (isv=0;isv<nrStates/typeofstates;isv++)
		{ffscanf(); sxy[isv]=atof(sa);} 
	for (isv=0;isv<nrStates/typeofstates;isv++)
		{ffscanf(); pr[isv]=atof(sa);}
	
	// Check if read all correctly: one more variable printed as 60.6060606060
	ffscanf(); 
	double okread = atof(sa);
	if(okread == 60.6060606060)
		{fprintf(fp_log,"New state variables from data assimilation update read correctly!\n");}
	else
		{
			fprintf(fp_log,"ERROR: Problem with reading of new state variables from data assimilation! Since okread = %f\n",okread); 
			exit(1);
		}
	
	fclose(fl);
	
	// ATAT SDA Check if no shift in location. Should not be I guess. But good to be 100% sure.
	
	// Calculate strain rates from velocities on nodes, since these used for slip rate calculation below
	// As in sxxcalc() and sxycalc() in move.c
	long int v[4];
	for (int m1=1;m1<xnumx;m1++)
	{
		for (int m2=1;m2<ynumy;m2++)
		{
			// Node number, i.e., position in vx[], etc. 
			long int mcmax1=m1*ynumy+m2;
				
			// Get neighbouring nodal numbers
			v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
			v[2]=v[0]+ynumy;v[3]=v[2]+1;
		
			// ATAT SDA Go over again why calculate like this. See klad 28-07-2015. Took from move.c straight... Deviatoric version???
			// Calculate Exx=1/2(dVx/dX-dVy/dY)
			exx[mcmax1]=0.5*((vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])-(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]));
		
			// Calculate Exy=1/2(dVx/dY+dVy/dX)=0? = p59 eq, w no rigid body rotation
			exy[mcmax1]=((vx[v[3]]-vx[v[2]])/(gy[m2+1]-gy[m2-1])) + ((vy[v[3]]-vy[v[1]])/(gx[m1+1]-gx[m1-1]));
			// ATAT SDA Why not still * 0.5? 
		}
	}
		
	// --- Interpolate changes from nodes to markers for viscosity calc. ---
	// So can calculate viscosity on each marker in step 3, before Stokes solve in step 4
	// During Navier-Stokes solve will obtain conservation ok again
	// But what impact of velocity update? strain rate for stress prediction, so right stress_old and viscosity ..
	// Viscosity affects new velocity calculation, but direct influence enough??
	double EXX,EXY,SXX,SXY,PR;
	
	// From move.c l. 737-849:
	for (long int m6=0;m6<=marknum;m6++) 
	{
		// Check that exclude non-rock markers and those out of grid
		if (markx[m6]>0 && marky[m6]>0 && markx[m6]<xsize && marky[m6]<ysize)
		{
			// Interpolate Stresses 
			int m10=m1serch(markx[m6]);
			int m20=m2serch(marky[m6]);

			// marker elements are passed by reference & to update value directly
			//allintersdaomp(markx[m6],marky[m6],m10,m20,&markexy[m6],&markxy[m6],&markexx[m6],&markxx[m6],&markp[m6]);
			allintersdaomp(markx[m6],marky[m6],m10,m20,&EXY,&SXY,&EXX,&SXX,&PR);

			// Stress diffusion here, if necessary - see move.c l. 737-849
			// Not now at first test. Not done so often yet, so initially likely fine.
			
			// Update marker variables. Now directly from state variables out of data assimilation
			// Before in forward model only update change w.r.t. previous time step to reduce interpolation artifacts, e.g., markxx[m6]+=(SXX-SXXE) - alternatively output only change in parameters....
			// For now keep simple and interpolate straight from nodes to markers = what get out of assimilation, i.e., what need
			markexx[m6] = EXX;
			markexy[m6] = EXY;
			markxx[m6]  = SXX;
			markxy[m6]  = SXY;
			markp[m6] 	= PR;
							
			// Should also take care of delta vx/vy for inertial term, which is advected via markers?
			// Would create "waves" due to assimilation changes in vx/vy; not desired...
			// --> run without inertia during assimilation time step
			//markvx[m6]+=(VX-MVX)*timestep/timestepe0;
			//markvy[m6]+=(VY-MVY)*timestep/timestepe0;
		}
	}

	// Test/Think after first tries about:
	// X Also miss stress rotation... important to still do? have vorticity from new vel.. or will be done afterwards??: stress rotation in our type of small time step model is small
	// Also miss marker advection following velocity update, but if not assimilate often guess is ok: ok due to small time step

	
}


/* Calculation of state variables for data assimilation by interpolation */
void allintersdaomp(double x, double y,long int m10, long int m20, double *EXY,double *SXY,double *EXX,double *SXX,double *PR)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3;
	/* en-NormalisedDistance */
	double ival,e,n,xrat;

	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;

	/* Store Up Left Node X,Y Num for later re-usage */
	m1=m10;
	m2=m20;

	/* Check weighting for interpolation */
	xrat=2.0/3.0;
	if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
	if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
	if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
	if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;

	/* SIGxy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXY=*SXY=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);

	/* EPSxy Interpolate after interpolation weights */
	m3=m10*ynumy+m20;
	*EXY=(1.0-e)*(1.0-n)*exy[m3]+(1.0-e)*n*exy[m3+1]+e*(1.0-n)*exy[m3+ynumy]+e*n*exy[m3+ynumy+1];
	*SXY=(1.0-e)*(1.0-n)*sxy[m3]+(1.0-e)*n*sxy[m3+1]+e*(1.0-n)*sxy[m3+ynumy]+e*n*sxy[m3+ynumy+1];
	/* End SIGxy interpolation ------------------------ */

	/* SIGxx,SIGyy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXX=*SXX=*PR=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* P interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*EXX =(1.0-e)*(1.0-n)*exx[m3]+(1.0-e)*n*exx[m3+1]+e*(1.0-n)*exx[m3+ynumy]+e*n*exx[m3+ynumy+1];
	*SXX =(1.0-e)*(1.0-n)*sxx[m3]+(1.0-e)*n*sxx[m3+1]+e*(1.0-n)*sxx[m3+ynumy]+e*n*sxx[m3+ynumy+1];
	*PR  =(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];

	//	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
}
/* Calculation of state variables for data assimilation by interpolation */
