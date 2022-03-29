/* === Information needed on the fly: parameters & markers to be tracked === */
void flyinput_lab()
{
	int mm1,im,imt;
	double xloc_vb_trench, xloc_pick;	
	FILE *fileexists;
	FILE *prnexists;
	char filePRN[50];
	
	// --- Overwrite or add input ---
	// Write data here so in large-scale model can restart from prior prn while updating values here as they will overwrite initial once
	// Load additional input data written by python for potential on the fly changes
	load_flyinput_lab();

	// Settings that may vary through time (in the large-scale model)
	veldepfric = 1;		// Switch for STM part (off in long-term geodynamic part)
	inertyn = 1;		// Inertia Yes (1) or No (0)
	lambfld = 1; 		// pore fluid pressure ratio (so not 1-this as in older version!); not used in lab as no fluid markers

	// --- Calculate from input: ---

	// Location seismogenic zone
	start_updip     = (start_sez+before_trench)-half_range;
	end_updip       = (start_sez+before_trench)+half_range;
	start_downdip   = (end_sez+before_trench)-half_range;
	end_downdip     = (end_sez+before_trench)+half_range;
	
	// 1 cm is input from the lab; where they measured
	// -1, to correct for going from 0 to ny-1 in vx-array
	above_fric = (int)(0.01/(gy[n_glayer]-gy[n_glayer-1]));
	vb_layer = n_glayer-1-above_fric;

	slab_dip        = (slab_dip_deg*3.14159)/180;

	// *** Select or read markers to be tracked ***
	// number of markers to be tracked is set statically in head.c
	
	// If not at start of simulation (t=0 s) than reload markers to be followed from previous txt file 
#if (setup == 1)	
	sprintf(filePRN,"%s1006.prn",exp_name);
	prnexists = fopen(filePRN,"r");
	if (timesum>=1e-15 && prnexists!=NULL)
#elif
	if (timesum>=1e-15)
#endif
	{
		if (printmod) fprintf(fp_log,"This is a restart. Collecting markers to be tracked and appending info to previous files. \n");
		restart = 1;
	
		// Read them from dt = 0 file  
		// GPS on surface 
		fileexists = fopen(fileTxtOutputGPS,"r");
		if (fileexists)
		{
			if (printmod) fprintf(fp_log,"PAY ATTENTION : LOADING GPS MARKERS TO FOLLOW FROM %s = RIGHT FILE !? \n",fileTxtOutputGPS); fflush(fp_log);	
			
			// Note needs to be 'fl' for usage of ffscanf	
			fl = fopen(fileTxtOutputGPS,"rt");
			if (fl==NULL)
				{ fprintf(fp_log,"ERROR: Can not open gps marker file : %s \n ",fileTxtOutputGPS); fflush(fp_log); exit(0); }
			for (im=0;im<nr_gpsmarkers;im++)
				{ffscanf(); mgpstrack[im] = atoi(sa);}
			fclose(fl);
		}
		else
			{fprintf(fp_log,"ERROR: Can not open gps marker file %s, but it is a restart \n",fileTxtOutputGPS); fflush(fp_log);  exit(0);}
	
		// Physics on fault
		fileexists = fopen(fileTxtMarkerOutput,"r");
		if (fileexists)
		{
			if (printmod) fprintf(fp_log,"PAY ATTENTION : LOADING PHYSICS MARKERS TO FOLLOW FROM %s = RIGHT FILE !? \n",fileTxtMarkerOutput); fflush(fp_log);	
			
			// Note needs to be 'fl' for usage with ffscanf	                
			fl = fopen(fileTxtMarkerOutput,"rt");
			if (fl==NULL)
				{ fprintf(fp_log,"ERROR: Can not open physics marker file : %s \n ",fileTxtMarkerOutput); fflush(fp_log); exit(0); }
			for (im=0;im<nr_physmarkers;im++)
				{ ffscanf(); mphystrack[im] = atoi(sa); }
			fclose(fl);	
		}
		else
			{fprintf(fp_log,"ERROR: Can not open physics marker file %s, but it is a restart \n",fileTxtMarkerOutput); fflush(fp_log);  exit(0);}
	}
	else
	{
		// --- Select new surface (GPS) markers ---
		fprintf(fp_log,"Selecting markers to be followed... \n"); fflush(fp_log); 

		// Get locations for GPS surface markers assuming a line at t=0
		for (im=0;im<nr_gpsmarkers;im++)
		{
			xgps[im] = gelx0 + startgps + dxgps*im;
			ygps[im] = gely0 - (xgps[im]-gelx0)*tan(slab_dip);
			
			fprintf(fp_log,"GPS markers for surface tracking will be located near: im = %d : y = %f cm, x = %f \n",im,ygps[im],xgps[im]);
		}
		// Add marker for data assimilation data point. Number comes from enkf_auto.m.
		// Small, fast setup: x=0.269-0.0025 , y=0.070+0.0025
		// Regular setup: x=0.3760-0.0025, y=0.0760+0.0025; calc from surface + depth equation
#if (setup == 1)
		{
			xgps[nr_gpsmarkers-1]=0.3760-0.0025;
			ygps[nr_gpsmarkers-1]=0.0760+0.0025;
			
			fprintf(fp_log,"GPS markers for surface tracking will be located near: im = %d : y = %f cm, x = %f \n",nr_gpsmarkers-1,ygps[nr_gpsmarkers-1],xgps[nr_gpsmarkers-1]);
		}
#endif


		// Select specific marker numbers
		for (mm1=0;mm1<marknum;mm1++)
		{
			if (markt[mm1]==2) 	
			{	
				for (im=0;im<nr_gpsmarkers;im++)
				{
					if (markx[mm1]>xgps[im] && markx[mm1]<(xgps[im]+0.005) && marky[mm1]>(ygps[im]-0.005) && marky[mm1]<ygps[im]) 
					{
						// To homogenize distance to free surface select marker nearest to the surface/selected point
						// do nearest since main interest in displacement, not stresses etc that might possibly be affected through interpolation with air
						if (mgpstrack[im]==0 || marky[mm1] < marky[mgpstrack[im]])
						{
							mgpstrack[im]=mm1;
							if (printmod){printf(" pick gps markers for im = %d : y = %f cm, x = %f, markernum %d \n",im,marky[mm1],markx[mm1], mgpstrack[im]);}
							continue;
						}
					}	
				}
			}
		}
			
		// --- Select new fault (physics) markers ---
		for (mm1=0;mm1<marknum;mm1++)
		{
			// - Pick markers to track/follow -
			// Markers within gel (1-7)
			// All at same height above interface as horizontal velocity, which retrieved before in postproc ( 1 cm above interface for now )
			if (marky[mm1]>gy[vb_layer] && marky[mm1]<gy[vb_layer+1] && markt[mm1]==2) 	
			{	
				// m1) +2 cm of virtual trench				
				xloc_vb_trench = before_trench + (0.01 / (tan(slab_dip)));
				xloc_pick = xloc_vb_trench + 0.02; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
				{
					mphystrack[0]=mm1;
				}	
			
				// m2) -1 cm updip sez 				
				xloc_pick = before_trench + start_sez - 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[1]=mm1; }	
			
				// m3) +1 cm updip sez 				
				xloc_pick = before_trench + start_sez + 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[2]=mm1; }	
			
				// m4)  				
				xloc_pick = before_trench + start_sez + 0.5*(end_sez-start_sez); 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[3]=mm1; }	
			
				// m5) -1 cm downdip sez 				
				xloc_pick = before_trench + end_sez - 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[4]=mm1; }	
			
				// m6) +1 cm downdip sez 				
				xloc_pick = before_trench + end_sez + 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[5]=mm1; }	
			
				// m7)  15 cm before backstop		
				xloc_pick = before_trench + ((w_height/sin(slab_dip))-d_bstop);
				if (setup==1){xloc_pick = before_trench + end_sez + 0.1; } 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[6]=mm1; }	
			}

			// Markers within fbl (8-14)
			// All at same height  just inside the fbl (n=126) 
			if (marky[mm1]>gy[n_glayer] && marky[mm1]<gy[n_glayer+1] && markt[mm1]==5) 	
			{	
				// m8) +2 cm of virtual trench				
				xloc_vb_trench = before_trench + (0.01 / (tan(slab_dip)));
				xloc_pick = xloc_vb_trench + 0.02; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[7]=mm1; }	
			
				// m9) -1 cm updip sez 				
				xloc_pick = before_trench + start_sez - 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[8]=mm1; }	
			
				// m10) +1 cm updip sez 				
				xloc_pick = before_trench + start_sez + 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[9]=mm1; }	
			
				// m11)  				
				xloc_pick = before_trench + start_sez + 0.5*(end_sez-start_sez); 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[10]=mm1; }	
			
				// m12) -1 cm downdip sez 				
				xloc_pick = before_trench + end_sez - 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[11]=mm1; }	
			
				// m13) +1 cm downdip sez 				
				xloc_pick = before_trench + end_sez + 0.01 ; 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[12]=mm1; }	
			
				// m14)  15 cm before backstop		
				xloc_pick = before_trench + ((w_height/sin(slab_dip))-d_bstop);
				if (setup==1){xloc_pick = before_trench + end_sez + 0.1;} 
				if (markx[mm1]>xloc_pick && markx[mm1]<(xloc_pick+0.005)) 
					{ mphystrack[13]=mm1; }	
			}
		}
	}
	
    	// Check if all markers filled, else exit, so can ensure right output
	int nr_gpsmarkers_missing = 0;	
	for (imt=0;imt<nr_gpsmarkers;imt++)
	{
		if (mgpstrack[imt]==0)
		{
			nr_gpsmarkers_missing++;	
			printf("  WARNING: The following GPS marker number to be tracked is not assigned: %d \n",imt); 
			
			if (nr_gpsmarkers_missing>=3)
			{
				printf("  ERROR: At least three GPS marker numbers to be tracked are not assigned: fix output!\n"); 
				exit(0);
			}
		}
	}
	
	int nr_physmarkers_missing = 0;	
	for (imt=0;imt<nr_physmarkers;imt++)
	{
		if (mphystrack[imt]==0)
		{
			nr_physmarkers_missing++;	
			printf("  WARNING: The following physics marker number to be tracked is not assigned: %d \n",imt); 
			
			if (nr_physmarkers_missing>=3)
			{
				printf("  ERROR: At least three physics marker numbers to be tracked are not assigned: fix output!\n"); 
				exit(0);
			}
		}
	}
}


/* === Post-processing of output at each timestep for printing to text files === */
void postproc_lab()
{
	int pos,il,im,imt,ilx,ily,iln,event_nr_print;
	double timesum05,sxy_tmp,sxx_tmp,vis_s_tmp, sxy_atsxx, exy_atexx;
	FILE *fld;
	FILE *fldQ;
	FILE *fldQ2;
	FILE *fld2;
	FILE *fldG;
	FILE *fldG2;
	FILE *fldM;
	FILE *fldM2;
	FILE *fldP;
	FILE *fldP2;
	FILE *fldS;
	FILE *fldS2;
	FILE *fldE;
	FILE *fldE2;
	FILE *fldDRg;
	FILE *fldDRg2;
	FILE *fldDRfd;
	FILE *fldDRfd2;
	FILE *fldDRfs;
	FILE *fldDRfs2;

	if (printmod) fprintf(fp_log,"Saving data each time step to *txt ...\n"); fflush(fp_log);

	count1++;

	if (count1 == 1)
	{
		if (restart==0)
		{
			// Write initial data to file
			fld = fopen(fileTxtOutput,"w");
			fprintf(fld," Location each node, total: %ld \n",xnumx);
			for(il=0; il<xnumx; il++)
				{fprintf(fld,"%f ", gx[il]);}
			fprintf(fld," \n");
 	
			fprintf(fld," Horizontal velocity at y= %f cm, node %d \n",gy[vb_layer],vb_layer);
        
			fprintf(fld," #   t(s)  v_bottom_for_each_y_node(cm/s) \n");
			fclose(fld);
       
			// Need like this, since in we are going to append data in loop  below
			fldM = fopen(fileTxtMarkerOutput,"wt");
			for (im=0;im<nr_physmarkers;im++)
				{fprintf(fldM,"%d ", mphystrack[im]);}
			fprintf(fldM,"\n");
			fclose(fldM);
						
			fldG = fopen(fileTxtOutputGPS,"wt");
			for (im=0;im<nr_gpsmarkers;im++)
				{fprintf(fldG,"%d ", mgpstrack[im]);}
			fprintf(fldG,"\n");
			fclose(fldG);
			
			fldS = fopen(fileTxtOutputS,"w");
			fprintf(fldS," t(s)  sbrit_ave_sez_top   sii_for_each_y_node(Pa)  \n");
			fclose(fldS);
	
			fldP = fopen(fileTxtOutputP,"w");
			fprintf(fldP," t(s)  pressure_for_each_y_node(Pa) \n");
			fclose(fldP);
        
			fldE = fopen(fileTxtOutputE,"w");
			fprintf(fldE," t(s)  eii_for_each_y_node(Pa) \n");
			fclose(fldE);


			fldQ = fopen(fileTxtOutputQ,"w");
			fprintf(fldQ," t(s)  yieldstress_for_each_y_node(Pa) \n");
			fclose(fldQ);
		}

		// not for restart; had frame start nr in it anyway, so only overwrite if start at same frame (=ok) and event_nr would not add then
		fldDRg = fopen(fileTxtOutputDRg,"w");
		fprintf(fldDRg," event_nr time vx_allnodes acc_slip_allnodes stressii_allnodes vis_allnodes \n");
		fclose(fldDRg);

		fldDRfd = fopen(fileTxtOutputDRfd,"w");
		fprintf(fldDRfd," event_nr time vx_allnodes acc_slip_allnodes stressii_allnodes vis_allnodes \n");
		fclose(fldDRfd);
	
		fldDRfs = fopen(fileTxtOutputDRfs,"w");
		fprintf(fldDRfs," event_nr time vx_allnodes acc_slip_allnodes stressii_allnodes  vis_allnodes \n");
		fclose(fldDRfs);
	
		acc_slip_gel = malloc( sizeof(double) * xnumx);
		if(acc_slip_gel == NULL)
		{
			fprintf(stderr, "Out of memory. \n");
			exit(EXIT_FAILURE); 
		}

		acc_slip_fbls = malloc( sizeof(double) * xnumx);
		if(acc_slip_fbls == NULL) 
		{
			fprintf(stderr, "Out of memory. \n");
			exit(EXIT_FAILURE); 
		}

		acc_slip_fbld = malloc( sizeof(double) * xnumx);
		if(acc_slip_fbld == NULL) 
		{
			fprintf(stderr, "Out of memory. \n");
			exit(EXIT_FAILURE); 
		}
	}

	timesum05 = timesum + 0.5*timestep;


	// --- Horizontal velocity data ---
	fld2 = fopen(fileTxtOutput,"a");
	fprintf(fld2," %d %e    ",count1,timesum05);
	for(il=0; il<xnumx; il++)
		{fprintf(fld2,"%f ", vx[il*ynumy+vb_layer]*100);}
	fprintf(fld2," \n ");
	fclose(fld2);

	// Next lines are all taken at interface top! Not 1 cm above, as velocity above.
	// Values taken at staggered nodes in center of cells
	
	// --- Stress data ---
	// Calculate averages
	sbrit_ave = sbrit_ave / count_sezm;

	fldS2 = fopen(fileTxtOutputS,"a");
	fprintf(fldS2," %e %f  ",timesum05, sbrit_ave);
	for(il=0; il<xnumx; il++)
	{
		// Interpolate sxy, so at same location as sxx
		sxy_atsxx = 0.25*(sxy[il*ynumy+(n_glayer+1)]+sxy[il*ynumy+(n_glayer+1)+1]+sxy[il*ynumy+(n_glayer+1)+ynumy]+sxy[il*ynumy+(n_glayer+1)+ynumy+1]);
		
		fprintf(fldS2,"%f ", pow(sxx[il*ynumy+(n_glayer+1)]*sxx[il*ynumy+(n_glayer+1)]+sxy_atsxx*sxy_atsxx,0.5));
	}
	fprintf(fldS2," \n ");
	fclose(fldS2);

	// output for interpolated yield stresses
	fldQ2 = fopen(fileTxtOutputQ,"a");
	fprintf(fldQ2," %e  ",timesum05);
	for(il=0; il<xnumx; il++)
		{ fprintf(fldQ2,"%f ", sbritn[il*ynumy+(n_glayer+1)] ); }
	fprintf(fldQ2," \n ");
	fclose(fldQ2);


	// --- Pressure data ---
	fldP2 = fopen(fileTxtOutputP,"a");
	fprintf(fldP2," %e  ",timesum05);
	for(il=0; il<xnumx; il++)
		{fprintf(fldP2,"%f ", pr[il*ynumy+(n_glayer+1)] );}
	fprintf(fldP2," \n ");
	fclose(fldP2);

	// --- Strain rate data ---
	fldE2 = fopen(fileTxtOutputE,"a");
	fprintf(fldE2," %e  ",timesum05);
	for(il=0; il<xnumx; il++)
	{
		// Interpolate exy, so at same location as exx
		exy_atexx = 0.25*(exy[il*ynumy+(n_glayer+1)]+exy[il*ynumy+(n_glayer+1)+1]+exy[il*ynumy+(n_glayer+1)+ynumy]+exy[il*ynumy+(n_glayer+1)+ynumy+1]);
	
		fprintf(fldE2,"%e ", pow(exx[il*ynumy+(n_glayer+1)]*exx[il*ynumy+(n_glayer+1)]+exy_atexx*exy_atexx,0.5));
	}
	fprintf(fldE2," \n ");
	fclose(fldE2);

	// --- Marker track data ---
	fldM2 = fopen(fileTxtMarkerOutput,"a");
	for (imt=0; imt<nr_physmarkers; imt++)
		{ fprintf(fldM2,"%e %.13e %.13e %e %e %e %e %e %e \n", timesum05, markx[mphystrack[imt]], marky[mphystrack[imt]], markp[mphystrack[imt]], msbrit[mphystrack[imt]], markxx[mphystrack[imt]],markxy[mphystrack[imt]],markexx[mphystrack[imt]],markexy[mphystrack[imt]]); }
	fclose(fldM2);

	// --- GPS marker track data ---
	fldG2 = fopen(fileTxtOutputGPS,"a");
	for (imt=0; imt<nr_gpsmarkers; imt++)
		{ fprintf(fldG2,"%e %.13e %.13e %e %e %e %e %e \n", timesum05, markx[mgpstrack[imt]], marky[mgpstrack[imt]], markp[mgpstrack[imt]], markxx[mgpstrack[imt]],markxy[mgpstrack[imt]],markexx[mgpstrack[imt]],markexy[mgpstrack[imt]]); }
	fclose(fldG2);


	// --- Determine start of an event using velocity threshold ---
	// At same depth as in post-processing
	if (event == 0)
	{
		// Set slip to zero always and when start of event 
		memset( acc_slip_gel, 0, sizeof(double) * xnumx);
		memset( acc_slip_fbls, 0, sizeof(double) * xnumx);
		memset( acc_slip_fbld, 0, sizeof(double) * xnumx);

		// ATAT wrong event nr; printed wring (format?) and what have wrong as never Inside event ...
		for(il=nstart_gel; il<nend_gel; il++)
		{
			if (vx[il*ynumy+vb_layer] < vtresh)
			{
				event = 1;
				event_nr = event_nr + 1;
				printf(" Start event %d ! \n", event_nr);

				break; 
			}
		}
	}
	else if (event == 1)
	{
		for(il=nstart_gel; il<nend_gel; il++)
		{
			// Continue event
			if (vx[il*ynumy+vb_layer] < vtresh)
			{ 
				event = 1; 
				printf(" Inside event ... %d \n", event_nr);
				break; 
			}
			// Terminate event
			// Will break at 1 when still is going on, so end result will be 1 if still a part with velocity above threshold
			else
				{ event = 0; }
		}
	}

	// --- Dynamic Rupture-like data ---
	// Once for whole run; filename incl event_nr ; open and first write somewhere in here ; acc_slip=0 here (these 2 with counter)

	// If inside event accumulate slip (cm)
	// vx is one-sides displacement, vpush is constant (and +), so get slip (=difference) like this;
	// vx is negative during event, so correct for positive slip with *-1
	// start count only when start event ... ass before is negligible ; no
	if (event == 1)
	{
		for(il=nstart_gel; il<nend_gel; il++)
		{
			// Only add if qualified as in event, as derived on 1 cm top
			if (vx[il*ynumy+vb_layer] < vtresh)
				// already secured add slip if in same event above with memset=0 if new event 
			{ 	
				acc_slip_gel[il] = acc_slip_gel[il] + -1*(vx[il*ynumy+vb_layer]-vpush)*timestep; 
				acc_slip_fbls[il] = acc_slip_fbls[il] + -1*(vx[il*ynumy+n_glayer]-vpush)*timestep; 
				acc_slip_fbld[il] = acc_slip_fbld[il] + -1*(vx[il*ynumy+n_glayer+1]-vpush)*timestep; 
			}
		}
	}

	// --- gelatine layer ---
	// Always write data to also see 
	// - 1 cm above interface; inside gelatine, for comparison with lab -
	// all located at line of staggered nodes, though vx in center and stresses and pressures one horizontal grid cell shifted
	fldDRg2 = fopen(fileTxtOutputDRg,"a");

	if (event == 0)
		{ event_nr_print = 0;}
	fprintf(fldDRg2," %d %e \n", event_nr_print, timesum05);
	
	// vx - max particle velocity (cm/s)
	for(il=nstart_gel; il<nend_gel; il++)
		{fprintf(fldDRg2,"%f ", vx[il*ynumy+vb_layer]*100);}
	fprintf(fldDRg2," \n ");

	// accumulated slip (cm)
	for(il=nstart_gel; il<nend_gel; il++)
		{ fprintf(fldDRg2,"%e ", acc_slip_gel[il]*100); }
	fprintf(fldDRg2," \n ");

	// sii (Pa)
	for(il=nstart_gel; il<nend_gel; il++)
		// Make sure at same location for calculation invariants, so linearly interpolate shear properties from normal nodes to staggered nodes where normal properties are located
		// Source is output here of sxy to prn
	{
		sxy_tmp = 0.25*(sxy[il*ynumy+vb_layer] + sxy[il*ynumy+vb_layer+1] + sxy[il*ynumy+vb_layer+ynumy] + sxy[il*ynumy+vb_layer+ynumy+1]);
		fprintf(fldDRg2,"%f ", pow(sxx[il*ynumy+vb_layer]*sxx[il*ynumy+vb_layer]+sxy_tmp*sxy_tmp,0.5));
	}
	fprintf(fldDRg2," \n ");
	
	// viscosity (Pa s)
	// need to average for viscosity for shear stress in stokes and normal deviatoric stress in stokes as have interface with mainly shear
	// Make sure at same location for proper averaging, so linearly interpolate shear properties from normal nodes to staggered nodes where normal properties are located
	// though really smoothens distribution, since a difference of one vertical grid cell makes a huge difference !
	// does not matter as will always be the viscosity of gelatine
	for(il=nstart_gel; il<nend_gel; il++)
		{ fprintf(fldDRg2,"%e ", nd[il*ynumy+vb_layer]); }
	fprintf(fldDRg2," \n ");

	fclose(fldDRg2);


	// - at top nodes fbl for what do in DR models -
	// really at top fbl, where still moving to trench all (= positive side there fault)
	// all located at line of staggered nodes
	fldDRfs2 = fopen(fileTxtOutputDRfs,"a");
	if (fldDRfs2==NULL) 
		{ printf("ERROR: Can not open DR file \n "); exit(0);}

	if (event == 0)
		{ event_nr_print = 0;}
	fprintf(fldDRfs2," %d %e \n", event_nr_print, timesum05);
	
	// vx - max particle velocity (cm/s)
	for(il=nstart_gel; il<nend_gel; il++)
		{fprintf(fldDRfs2,"%f ", vx[il*ynumy+n_glayer]*100);}
	fprintf(fldDRfs2," \n ");

	for(il=nstart_gel; il<nend_gel; il++)
		{ fprintf(fldDRfs2,"%e ", acc_slip_fbls[il]*100); }
	fprintf(fldDRfs2," \n ");

	// sii (Pa)
	// Make sure at same location for calculation invariants, so linearly interpolate shear properties from normal nodes to staggered nodes where normal properties are located
	// Source is output here of sxy to prn
	// Make sure is located exactly at top of frictional boundary layer and normal nodes (not staggered that ussualy do)
	for(il=nstart_gel; il<nend_gel; il++)
	{
		sxx_tmp = 0.25*(sxx[il*ynumy+n_glayer] + sxx[il*ynumy+n_glayer+1] + sxx[il*ynumy+n_glayer-ynumy] + sxx[il*ynumy+n_glayer+1-ynumy]);
		fprintf(fldDRfs2,"%f ", pow(sxy[il*ynumy+n_glayer]*sxy[il*ynumy+n_glayer]+sxx_tmp*sxx_tmp,0.5));
	}
	fprintf(fldDRfs2," \n ");

	// viscosity (Pa s)
	// need to average for viscosity for shear stress in stokes and normal deviatoric stress in stokes as have interface with mainly shear
	// Make sure at same location for proper averaging, so linearly interpolate shear properties from normal nodes to staggered nodes where normal properties are located
	// though really smoothens distribution, since a difference of one vertical grid cell makes a huge difference !
	// Make sure is located exactly at top of frictional boundary layer and normal nodes (not staggered that ussualy do)
	for(il=nstart_gel; il<nend_gel; il++)
	{ 
		// only at n_glayer+1 are in a zone where viscosity is really reduced
		// use only nd as log averaging is not good, and anyway is viscosity at that exact location, just used for solving exx part of stokes
		fprintf(fldDRfs2,"%e ", (nd[il*ynumy+n_glayer]));
	}
	fprintf(fldDRfs2," \n ");

	fclose(fldDRfs2);

	
	// - at top nodes fbl for what do in DR models -
	// one cell lower, so that in main slip zone, so see large viscosity drop
	// all located at line of staggered nodes
	fldDRfd2 = fopen(fileTxtOutputDRfd,"a");
	if (fldDRfd2==NULL)
		{printf("ERROR: Can not open DR file \n "); exit(0);}

	if (event == 0)
		{event_nr_print = 0;}
	fprintf(fldDRfd2," %d %e \n", event_nr_print, timesum05);
	
	// vx - max particle velocity (cm/s)
	for(il=nstart_gel; il<nend_gel; il++)
		{fprintf(fldDRfd2,"%f ", vx[il*ynumy+n_glayer+1]*100);}
	fprintf(fldDRfd2," \n ");

	// accumulated slip (cm)	
	for(il=nstart_gel; il<nend_gel; il++)
		{ fprintf(fldDRfd2,"%e ", acc_slip_fbld[il]*100); }
	fprintf(fldDRfd2," \n ");

	// sii (Pa)
	for(il=nstart_gel; il<nend_gel; il++)
	{
		sxx_tmp = 0.25*(sxx[il*ynumy+n_glayer+1] + sxx[il*ynumy+n_glayer+2] + sxx[il*ynumy+n_glayer+1-ynumy] + sxx[il*ynumy+n_glayer+2-ynumy]);
		fprintf(fldDRfd2,"%f ", pow(sxy[il*ynumy+n_glayer+1]*sxy[il*ynumy+n_glayer+1]+sxx_tmp*sxx_tmp,0.5));
	}
	fprintf(fldDRfd2," \n ");

	// viscosity (Pa s)
	for(il=nstart_gel; il<nend_gel; il++)
		{fprintf(fldDRfd2,"%e ", (nd[il*ynumy+n_glayer+1]));}
	fprintf(fldDRfd2," \n ");

	fclose(fldDRfd2);
}
