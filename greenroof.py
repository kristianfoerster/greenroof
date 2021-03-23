# greenroof.py provides a generic class for 1D flow in green roofs with
# arbitrary dimensions
#
# Authors: Kristian Foerster, Philipp Kraft

import cmf
import pandas as pd
import numpy as np

class GreenRoof(cmf.project):
    """
    The GreenRoof class is a CMF implementation to compute Darcy and Richards
    flow in a green roof.
    """
    def __init__(self, init_state=None,dt_min=0.25, width=1., length=20, nl=5, 
                 dz=0.02, rheight=0.,depth=0.08, ksat = 1500, 
                 porosity = 0.5, b=5, mannings_n=0.08, duration=180, 
                 d_init_pot=0.3):
        """
        Constructor method of the GreenRoof class, which creates a model,
        based on a number of paramters.
        
        Parameters
        ----------
        init_state : List, optional
            List of two arrays (1st: x positioning of observations, 1nd: water level observed at this position).
            The default is None.
        dt_min : float, optional
            Model time step [min]. The default is 0.25.
        width : float, optional
            Width of the model [m]. The default is 1..
        length : float, optional
            Length of the model [m]. The default is 20.
        nl : float, optional
            Number of nodes / cells per meter. The default is 5.
        dz : float, optional
            Thickness of numerical cells [m], which together with depth determines the vertical resolution. The default is 0.02.
        rheight : float, optional
            Height of right side (left is 0), which is length*slope. The default is 0..
        depth : float, optional
            Vertical extent of the green roof [m]. The default is 0.08.
        ksat : float, optional
            Hydraulic conductivity :math:`K_{s}` [m/d]. The default is 1500.
        porosity : float, optional
            Porosity of the substrate [-]. The default is 0.5.
        b : float, optional
            Parameter b (Brooks-Corey retention curve) [-]. The default is 5.
        mannings_n : float, optional
            Manning's roughness for surface flow [-]. The default is 0.08.
        duration : int, optional
            Duration of the simulation period [min]. The default is 180.
        d_init_pot : float, optional
            Delta in matric potential subtracted from the initial value (0.3 would subtract 0.3 from the standard value in each layer) [m]. The default is 0.3.

        Returns
        -------
        None.

        """
                
        # call cmf project constructor
        super().__init__()

        # model setup
        self._set_time(dt_min,duration)
        self._set_dimensions(width, length, nl, dz, rheight, depth)
        self._define_geometry(init_state)
        self._set_parameters(ksat,porosity,b,mannings_n)
        self._create_cells(d_init_pot)
        self._define_fluxes()
                
    def _set_time(self,dt_min,duration):
        """
        Sets the timing of the model (temporal resolution and duration)

        Parameters
        ----------
        dt_min : float
            Time step [min].
        duration : int
            Duration of the simulation [m].

        Returns
        -------
        None.

        """
        self.tstart = cmf.Time(1, 7, 2019)
        self.dt = cmf.min * dt_min
        self.dt_min = dt_min
        self.duration = cmf.min*duration

    def _set_dimensions(self,width,length,nl,dz,rheight,depth):
        """
        Sets the dimensions of the model (horizontal and vertical)

        Parameters
        ----------
        width : float, optional
            Width of the model [m]. The default is 1..
        length : float, optional
            Length of the model [m]. The default is 20.
        nl : float, optional
            Number of nodes / cells per meter. The default is 5.
        dz : float, optional
            Thickness of numerical cells [m], which together with depth determines the vertical resolution. The default is 0.02.
        rheight : float, optional
            Height of right side (left is 0), which is length*slope. The default is 0..

        Returns
        -------
        None.

        """
        self.width = width
        self.length = length
        self.ncell = int(length * nl) # nl ... # of nodes per meter length
        self.rheight = rheight
        self.depth = depth
        self.slope = self.rheight/self.length
        self.slope_per_cell = self.rheight/self.ncell
        self.gradient = self.rheight/self.length
        self.dz = dz

    def _define_geometry(self,init_state):
        """
        This function computes the numerical grid for which flow is solved.

        Parameters
        ----------
        init_state : List, optional
            List of two arrays (1st: x positioning of observations, 1nd: eater level observed at this position).

        Returns
        -------
        None.

        """
        self.x_dim = np.arange(0, self.length, self.length / self.ncell)
        self.bottom = np.linspace(0,self.rheight, self.ncell, endpoint=True)
        
        if init_state is not None:
            xdim_c = init_state[0] / 20. * self.length
            self.init_z_interpol = np.interp(self.x_dim, xdim_c,init_state[1])
        else:
            self.init_z_interpol = np.zeros(len(self.x_dim))
        
        self.init_z_interpol+=self.bottom        
        
    def _set_parameters(self, ksat = 1500, porosity = 0.68,
                 b=5,mannings_n=0.08):
        """
        This function sets all model parameters required for simulations.

        Parameters
        ----------
        ksat : float, optional
            Hydraulic conductivity :math:`K_{s}` [m/d]. The default is 1500.
        porosity : float, optional
            Porosity of the substrate [-]. The default is 0.5.
        b : float, optional
            Parameter b (Brooks-Corey retention curve) [-]. The default is 5.
        mannings_n : float, optional
            Manning's roughness for surface flow [-]. The default is 0.08.

        Returns
        -------
        None.

        """
        self.Ksat = ksat
        self.porosity = porosity 
        self.mannings_n = mannings_n 
        self.b = b

    def _create_cells(self,d_init_pot=0.3):
        """
        Creates an obejct for each cell which consists of vertical layers

        Parameters
        ----------
        d_init_pot : float, optional
            Delta in matric potential subtracted from the initial value (0.3 would subtract 0.3 from the standard value in each layer) [m]. The default is 0.3.

        Returns
        -------
        None.

        """
        for i in range(self.ncell):
            c: cmf.Cell = self.NewCell(i * self.length / self.ncell, 0, 
                                       self.depth + (i-0) * self.slope_per_cell,
                                       self.length / self.ncell * self.width, True)                        
            if i:
                c.topology.AddNeighbor(self[i - 1], self.width)
            self.retention_curve = cmf.BrooksCoreyRetentionCurve(ksat=self.Ksat,
                                                                 porosity=self.porosity,
                                                                 _b=self.b,
                                                                 theta_x=0.20,
                                                                 psi_x=cmf.pF_to_waterhead(2.5))

            list_z = list()
            for zi in np.arange(0,self.depth, self.dz):
                c.add_layer(zi+self.dz, self.retention_curve)
                list_z.append(zi+self.dz)
            c.install_connection(cmf.Richards)
            c.saturated_depth=self.depth-(self.init_z_interpol-self.bottom)[i]
            for ii,li in enumerate(c.layers):
                if list_z[ii] <= self.depth-(self.init_z_interpol-self.bottom)[i]:
                    li.potential -= d_init_pot
            del list_z
            c.surfacewater.puddledepth = 0.001
            c.surfacewater.nManning = self.mannings_n        

    def _define_fluxes(self):
        """
        Defines process description for fluses, meteorological forcing and boundary conditions (outlets).

        Returns
        -------
        None.

        """
        # meteorology
        self.rainfall_stations.add('', 0.0, (0, 0, 0))
        self.use_nearest_rainfall()
        self.ET = []

        # cell to cell flux
        cmf.connect_cells_with_flux(self, cmf.Richards_lateral)
        cmf.connect_cells_with_flux(self, cmf.DiffusiveSurfaceRunoff)
        
        # define outlets
        # outlet Darcy flow
        self.outlet = self.NewOutlet('out', -self.length / self.ncell, 0, 0-self.gradient*self.length/self.ncell)
        # outflow surface runoff
        self.outlet_s = self.NewOutlet('out_s', -self.length / self.ncell, 0, 0+self.depth-self.gradient*self.length/self.ncell)

        for li in self[0].layers:
            cmf.Darcy(li, self.outlet, self.width)
        cmf.DiffusiveSurfaceRunoff(self[0].surfacewater, self.outlet_s, self.width)

    def set_design_rain(self, rain_duration=15, rain_amount=27):
        """
        Defines a uniform rainfall input.

        Parameters
        ----------
        rain_duration : int, optional
            Rainfall duration [min]. The default is 15.
        rain_amount : float, optional
            Rainfall total [mm]. The default is 27.

        Returns
        -------
        None.

        """
        data = cmf.timeseries(self.tstart, self.dt)
        while data.end < self.tstart + rain_duration*cmf.min:
            data.add(rain_amount * cmf.day/(rain_duration* cmf.min))
        data.add(0)
        self.rainfall_stations[0].data = data

    def set_rain(self, fromarray):
        """
        Sets rainfall from an array of values (temporal resolution 1 minute).

        Parameters
        ----------
        fromarray : np.array
            Rainfall intensities [mm/min].

        Returns
        -------
        None.

        """
        data = cmf.timeseries.from_array(self.tstart,cmf.min,fromarray*1440)
        self.rainfall_stations[0].data = data

    def _to_liters_per_timestep(self,value):
        """
        Transforms the standard flux values in m^3/d (as computed by CMF) to liters per time step

        Parameters
        ----------
        value : float
            Flux [m^3/d].

        Returns
        -------
        float
            Flux [L/time step].

        """
        return value * 1000 / 1440  * self.dt_min
        
    def run(self, verbose=False):
        """
        Runs the model.

        Parameters
        ----------
        verbose : bool, optional
            Whether progress is printed. The default is False.

        Returns
        -------
        resdf : pd.DataFrame
            Model results provided as time series (DatetimeIndex)
            Column 0 ('Darcy flow'): Modelled outflow (Darcy flow) [mm/min].
            Column 1 ('Surface runoff'): Modelled outlow (Surface runoff) [mm/min].
            Column 2 ('Water level'): Time series of lists including the water level along horizontal axis [m].
            Column 3 ('Ponding'): Time series of lists including the ponding on the surface along horizontal axis [m].

        """
        solver = cmf.CVodeIntegrator(self, 1e-9)
        list_t = list()
        list_q = list()
        list_qs = list()
        list_wlevel = list()
        list_sd = list()
        for t in solver.run(self.tstart, self.tstart + self.duration, self.dt):
            if verbose: print(t)
            list_q.append(self._to_liters_per_timestep(self.outlet(t)))
            list_qs.append(self._to_liters_per_timestep(self.outlet_s(t)))
            list_t.append(t.as_datetime())
            wlevel = np.zeros(self.ncell)
            sd = np.zeros(self.ncell)
            for ii in range(len(self.cells)):
                wlevel[ii] = self.bottom[ii]+max(0,self.depth-self.cells[ii].layers[0].get_saturated_depth())
                sd[ii]     = self.cells[ii].surfacewater.depth
            list_wlevel.append(wlevel)
            list_sd.append(sd)
        resdf = pd.DataFrame(index=list_t, data={'Darcy flow':list_q,
                                                'Surface runoff':list_qs,
                                                'Water level':list_wlevel,
                                                'Ponding':list_sd})
        return resdf

    @property
    def potential(self):
        """
        Returns the potential in each cell.

        Returns
        -------
        List
            Potential along x axis [m].

        """
        return [self.outlet.potential] + [c.layers[0].potential for c in self]

    @property
    def x(self):
        """
        Returns the x coordinate for each cell (as distance from the outlet)

        Returns
        -------
        List
            Distance along x axis [m].

        """
        return [self.outlet.position.x] + [c.x for c in self]

        