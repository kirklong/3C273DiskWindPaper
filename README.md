# 3C273DiskWindPaper
Code and related files for Long+2023

Note: There is an error on line 172 of `functions.jl` where the size is fit in terms of $r_g$ instead of $r_s$ (units published in paper). The function `get_rMinMax()` is also incorrect, as there was an extra factor of $r$ that was omitted in deriving it by hand. An erratum has been submitted further detailing and correcting these errors, but the main takeaway is that this combined these errors lead to the $\bar{r}$ in the fits is underestimated by a factor of $\sim 1.3$ for a characteristic $r_\rm{fac} \sim 45$. We have uploaded a [corrected corner plot](cornerUpdated.png) showcasing the true values for $\bar{r}$ in the fits.

The code here is left as is to preserve the scientific record, but note an updated version of the `DiskWind` module with these errors corrected is available at https://github.com/kirklong/DiskWind

  
