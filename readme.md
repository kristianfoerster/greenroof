# A GreenRoof model utilising the Catchment modelling Framework (CMF)

This repository includes a Python implementation for a numerical flow model representing green roofs (Förster et al, in prep.). The class is inherited from a Catchment Modelling Framework project (CMF), which is a framework to develop hydrological models of arbitrary complexity (Kraft et al., 2011). The GreenRoof class represents Darcy and Richards flow and surface runoff (diffusive wave approximation). See [requirements.txt](requirements.txt) for a list of required Python packages. A [Jupyter notebook](Greenroof_test.ipynb) demonstrates how the GreenRoof class works in principle.

## References
Förster, K., Westerholt, D., Kraft, P., Lösken, G. (in prep.). Unprecedented retention capabilities of extensive green roofs – New design approaches and an open-source model.

Kraft, P., Vaché, K. B., Frede, H.-G. and Breuer, L.: CMF: A Hydrological Programming Language Extension For Integrated Catchment Models, *Environmental Modelling & Software*, 26(6), 828–830, https://doi.org/10.1016/j.envsoft.2010.12.009, 2011.
