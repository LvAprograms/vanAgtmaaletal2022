void load_flyinput_lab()
{
       // Velocity-dependent frictional parameters
       // gamma = 1 - mu_d / mu_s , but also (beta-alpha)/mu_s in eq 3 in Ampuero & Ben-Zion (2008) -> mu_d = ( 1 - gamma ) * mu_s

       // Amount of Velocity Weakening / Strengthening (gamma) 
	    mgamma_vs = -77.5;
	    mgamma_vw = 0.825;

       // Characteristic velocity, slip rate at which friction is increased by a half
	    mvc_vs = 3.9e-05;
	    mvc_vw = 0.0002;

       // Static friction coefficient. mus_vw is set on the fly depending on rock type. 
	    mus_vs = 0.002;

        // Location seismogenic zone
	    before_trench   = 0.0827; 		// before trench
	    start_sez       = 0.13830;             // updip limit of seismogenic zone
	    end_sez         = 0.29830;            // downdip limit of seismogenic zone
	    half_range      = 0.002;

        // Wedge geometry parameters for selecting (surface) markers
	    gelx0           = 0.0847;       // x-coordinate of point 0 of the gel-wedge as defined in the init.t3c file
	    gely0           = 0.1193;       // y-coordinate of point 0 of the gel-wedge as defined in the init.t3c file

	    w_height        = 0.11;              // height of the wedge
	    d_bstop         = 0.15;    		// distance from the backstop for marker m7 (track markers:figure 7 in vanDinther(2012))
	    slab_dip_deg    = 10;          	// slab dip in degrees
	    vpush           = 3.900e-05;   		// m/s
	    vtresh          = -7.5e-05;    		// m/s, official: 1.5e-04, but use this to output more so see more of event

        // Location values for GPS markers
	    startgps        = 0.0268;		// relative to gelx0
	    dxgps           = 0.05;		// distance between gps markers

        // Get nodal point number for analysis layer
	    n_glayer 	    = 124;                // y-node of the top of the friction boundary layer

	    nstart_gel	    = 83; 
	    nend_gel 	    = 691; 
}
