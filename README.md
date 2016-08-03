# PedyR

R package to deal with genealogical structured data (pedigrees)

*next updates:*
  * **getOrdPed()** - A function to order a pedigree using the Kahn's topological sorting algorithm 

## General Functions (version 0.1.3)

For now, PedyR package has the following functions:

### checkPed()

Checks for the (most common) errors in pedigrees, such as duplicated ids, sire and dams with the same ids, etc...



### getF()

Calculates the inbreeding coefficient (see *Meuwissen and Luo, 1992*) for each individual in the pedigree.



### getGen()

Calculates the generation of each individual in the pedigree.



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

THE Meuwissen and Z Luo, 1992. Computing inbreeding coefficients in large populations. Genet. Sel. Evol. 24, 305-313.
