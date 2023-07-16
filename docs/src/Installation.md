# Installation

We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:

```bash
git pull https://github.com/singular-gpispace/GromovWitten.git
```

In the same folder execute the following command:

```bash
julia --project
```

This will activate the environment for our package. In Julia install missing packages:

```bash
import Pkg; Pkg.instantiate()
```

and load our package. On the first run this may take some time.

```bash
using GromovWitten  
```
