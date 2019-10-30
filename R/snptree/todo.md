
Improvements
============

* Integration with `fastPHASE` and `shapeit`. 
    * The way to go about this might be to use `rehh`, either borrowing code or importing the library
* Integration with` haplostrips` so that I can read `vcf` files and do the same with that
    * This is proving troublesome as the python dependencies 
    are a little troublesome.  I could try to get a docker to make it work but that seems even more effort
* Integration with `ms` for simulations
* Integration with my own simulations
* Get plotting working in both directions
* Can I integrate the PAC idea into tree building a tree viewing

## List

1. See how haplostrips deals with vcf files
2. Get my own simulation system working with `bifurcation_split`
3. join left and right trees in a single tree removing the centre
     * write a function for left and right trees together
     * write a function that splits left and right trees by the centre position
4. Let's see if I can't rework the plot.phylo code to just do one thing.


## Talk
In my the usual way, the way to do this is to get a talk prepared on the material.

1. Genetic information can be represented by trees

1. How do we represent the information in haplotypes with:
	* extra information on origin?
	* case control status?
	* other covariates? 
	
2.  Some examples
	* bifurcation diagram
	* haplostrips
	
3. Can we do better with models?