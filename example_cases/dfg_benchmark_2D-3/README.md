For an intro to this benchmark setup, please refer the [DFG benchmark 2D-3](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark3_re100.html) page.

The drag and lift coefficients are defined as follows :

- Cd = 2*assemble(dot(traction, self.nx)*ds(5))
- CL = 2*assemble(dot(traction, self.ny)*ds(5))
