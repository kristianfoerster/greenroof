# A GreenRoof model utilising the Catchment modelling Framework (CMF)

It is published as free software under [GPLv3](http://www.gnu.org/licenses/gpl.html). (c) 2021 by Kristian Förster and Philipp Kraft.

This repository includes a Python implementation for a numerical flow model representing green roofs (Förster and Westerholt et al., 2021). The class is inherited from a Catchment Modelling Framework project ([CMF](https://philippkraft.github.io/cmf/)), which is a framework to develop hydrological models of arbitrary complexity (Kraft et al., 2011). The GreenRoof class represents Darcy and Richards flow and surface runoff (diffusive wave approximation). See [requirements.txt](requirements.txt) for a list of required Python packages. A [Jupyter notebook](Greenroof_test.ipynb) demonstrates how the GreenRoof class works in principle. The measurements were kindly provided by Daniel Westerholt (Westerholt and Lösken, 2021).

## References
Förster, K., Westerholt, D., Kraft, P., Lösken, G. Unprecedented retention capabilities of extensive green roofs – New design approaches and an open-source model, *Frontiers in Water*, 3, 689679, [https://doi.org/10.3389/frwa.2021.689679](https://doi.org/10.3389/frwa.2021.689679), 2021.

Kraft, P., Vaché, K. B., Frede, H.-G. and Breuer, L.: CMF: A Hydrological Programming Language Extension For Integrated Catchment Models, *Environmental Modelling & Software*, 26(6), 828–830, [https://doi.org/10.1016/j.envsoft.2010.12.009](https://doi.org/10.1016/j.envsoft.2010.12.009), 2011.

Westerholt, D., Lösken, G. Artificial rainfall experiments for green roofs depending on the flow length at 0 % and 2 % slope and different rain intensities (Version 1) [Data set]. *Zenodo*. [https://doi.org/10.5281/zenodo.4651339](https://doi.org/10.5281/zenodo.4651339), 2021.