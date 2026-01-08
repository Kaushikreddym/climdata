## Creation of your own CPP/C Library : 
<table>
    <tr>
        <td>requirements windows</td>
        <td>Armadillo, MinGW</td>
    </tr>
    <tr>
        <td>requirements linux/mac</td>
        <td>build-essential, g++, libarmadillo-dev</td>
    </tr>
</table>


### Generation

Navigate to the C++ directory:
```sh
$ cd ./wrapper/AlgoCollection
```
Generate your shared object file for ImputeGAP:
```sh
$ make libCDREC.so
```

### List of Available Algorithms:
- `libSTMVL.so`
- `libIterativeSVD.so`
- `libGROUSE.so`
- `libDynaMMo.so`
- `libROSL.so`
- `libSoftImpute.so`
- `libSPIRIT.so`
- `libSVT.so`
- `libTKCM.so`
- `libCDREC.so`
