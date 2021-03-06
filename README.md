# PedyR

R package to deal with pedigrees.

*upcoming updates:*
  * pedStats() - A function for more in depth description of a pedigree.
  * getBdate() - Estimates a birth date for animals without that record in the pedigree.
  * getFounder() - Calculates the effective number of founders for a pedigree.
  
## General Functions (version 0.1.4 - Aug 8, 2016)

For now, PedyR package has the following functions:

### checkPed()

Checks for the (most common) errors in pedigrees, such as duplicated ids, sire and dams with the same ids, etc...


### getF()

Calculates the inbreeding coefficient [[1](http://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-24-4-305)] for each individual in the pedigree.


### getGen()

Calculates the generation of each individual in the pedigree.


### getOrdPed()

A function to order a pedigree using the Kahn's topological sorting algorithm [[2](http://dl.acm.org/citation.cfm?id=369025)]. 


### effSize()

Calculates the effective population size for a given pedigree using methods discribed in [[3](http://pendientedemigracion.ucm.es/info/prodanim/html/JP_Web_archivos/jbg810.pdf)]. 
## Example:

```
id = c(1,2,3,4,5,6,7,8,9,10)
sire = c(0,0,1,1,3,5,5,5,7,7)
dam = c(0,0,2,2,2,4,4,6,6,8)
ped = as.data.frame(cbind(id,sire,dam))
checkPed(ped)
ped$F = getF(ped)
ped$Gen = getGen(ped)
```

## Authors

* **Bruno C. Perez** - [GitHub](https://github.com/BrnCPrz)
* **Gerson A. Oliveira Junior** - [GitHub](https://github.com/gersonjr)


## Contributing

Please feel free to submit new ideas, code corrections or pull requests to us.


## License

This project is licensed under the MIT License.

## References

 * [[1](http://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-24-4-305)] Meuwissen THE and Luo Z, 1992. Computing inbreeding coefficients in large populations. Genet. Sel. Evol. 24, 305-313.
 * [[2](http://dl.acm.org/citation.cfm?id=369025)] Kahn AB, 1962. "Topological sorting of large networks", Communications of the ACM.
 * [[3](http://pendientedemigracion.ucm.es/info/prodanim/html/JP_Web_archivos/jbg810.pdf)] Gutiérrez, JP, Cervantes I, Goyache F, (2009). "Improving the estimation of realized effective populationsizes in farm animals", J. Anim. Breed. Genet., 126, 327-332.
