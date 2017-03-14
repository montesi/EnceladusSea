# EnceladusSea
Calculation and analysis of the stress field expected in an ice shell pressurized by a regional sea confined at a pole or a global internal ocean with thin ice region at one or both poles. See Johnston and Montési, Journal of Geophysical Research, 2017.

* **EnceladusBatch**: run *EnceladusBatch* to automatically setup all models and analyse the stress at the surface of the model. This calls enceladusBuild
* **EnceladusBuild**: run *[out, Label] = EnceladusBuild(ThicknessShell,SeaThick,ModelType,nA)* to generate a series of models with given shell thickness *ThicknessShell*, and south polar sea/indendation thickness *Seathick* for a parametric sweep of *nA* values of sea/indentation angle equally spaced from 0° to 90°. *ModelType* is a flag for the desired model configuration
  * *ModelType=1*: no slip at the base of the ice shell (fixed)
  * *ModelType=2*: free slip at the base of the ice shell (roller)
  * *ModelType=3*: floating ice shell with an indentation only at the South Pole (ocean)
  * *ModelType=4*: floating ice shell with indentations at both Pole (north)
* **makefigureEnceladus**: Run *makefigureEnceladus(model,Label)* to visualize the stress field in the ice shell for a specific model; *model* is a comsol structure (figure 4 in Johnston and Montési 2017).
* **stressProfileEnceladus**: run *[Tectonics, MidFail]=stressProfileEnceladus(model,Yshear,Ycrack,Label)* to extract the stress profile at the surface and generate most figures in the paper. *model* is a comsol structure, *Yshear* and *Ycrack* the von Mises strenght and tensile strength.
* **classifyTectonics**: run *[Tectonics,Asouth,Anorth,MidFail]=classifyTectonics(Crack, Fmode, Regime, Alind)* to analyse the surface stress profile and classify the tectonic regime according the codes defined in Johnston and Montési (2017), figure 7.
* **basalStressEnceladus**: run this script to extract the profile of stress at the base of the ice shell and report the stress at a short distance of the sea edge as a function of sea angle.
* **NorthRelate**; run this script to calculate the north polar indentation angle as a function of the south polar indentation angle, ice shell thickness, and south polar indentatino thickness. Does not use COMSOL Multiphysics.
